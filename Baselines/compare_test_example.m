function compare_test_example()
    % Step 0: Environment Setup
    %initCobraToolbox(false);  
    %startcna;
    %setJavaHeapSize(8); 
    %parpool('local', 12);

    % Load model and initialize data
    load('e_coli_core.mat');
    ori_model = e_coli_core;
    growth_idx = find(ori_model.c);
    initial_remaining_gene_pool = readtable('e_coli_core_ex.csv', 'Delimiter', ',').Remaining_gene.';
    map_gene_list = readtable('e_coli_core_RID.csv', 'Delimiter', ',');
    D_list = map_gene_list{:, 1};      
    R_names = map_gene_list{:, 2}; 
    test_target = 'succ_e';
    
    % Create or open the output text file
    output_file = 'compare_test_example_results.txt';
    fid = fopen(output_file, 'w');  % Open for writing (creates the file if it doesn't exist)
    
    if fid == -1
        error('Failed to open file for writing.');
    end

    % Step 1: DBgDel
    overallTimer = tic;
    [gvalue, gr, pr, it, success] = DBgDel(ori_model, test_target, 10, 0.001, 0.001, initial_remaining_gene_pool);
    overallElapsedTime = toc(overallTimer);
    
    % Run the GRPRchecker function
    [GR_1, PR_1] = GRPRchecker(ori_model, test_target, gvalue);
    
    % Step 2: gMCSE
    [model, targetRID, extype] = modelSetting(ori_model, test_target);
    cnap = CNAcobra2cna(model);
    [cnap, enzymes, genes, gpr_rules] = CNAgenerateGPRrules(cnap);
    [T, t] = initializeTargetRegion(cnap.numr, growth_idx, targetRID, 0.001, 0.001);
    [D, d] = initializeDesiredRegion(cnap.numr, growth_idx, 0.001);
    options.milp_time_limit = 600;
    % Start timer for gMCSE
    tic;
    % Call the CNAgeneMCSEnumerator2 function
    [rmcs, full_mcs, full_cnap, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, status, obj] = ...
        CNAgeneMCSEnumerator2(cnap, T, t, D, d, [], [], ...
                              1, [], [], [], gpr_rules, options, []);
    % Initialize deleted_gene_str as an empty string by default
    deleted_gene_str = '';
    % Extract deleted genes **only if status is 0 or 3**
    if status == 0 || status == 3
        if size(full_mcs, 2) >= 1  % Check if full_mcs has at least one column
            [rowIdx, ~, val] = find(full_mcs(:,1)); % Find nonzero elements in column 1
            list_n = rowIdx(val == -1); % Extract indices where value is -1
            [~, idx] = ismember(list_n, D_list); % Find indices in D_list that match list_n
            deleted_gene = R_names(idx(idx > 0)); % Extract names, ignoring unmatched values

            % Convert deleted_gene to a single comma-separated string
            deleted_gene = strtrim(deleted_gene);
            deleted_gene_str = strjoin(deleted_gene, ', ');
        end
    end
    % End timer for gMCSE
    time_cost_gmcse = toc;
    % Get deleted genes
    geneD =deleted_gene_str; 
    % Store gene names and determine if they are in the deleted genes list
    for j = 1:length(ori_model.genes)
        gvalue{j, 1} = ori_model.genes{j};  % Store gene name
        if ismember(ori_model.genes{j}, strsplit(deleted_gene_str, ', '))
            gvalue{j, 2} = 0;  % Mark as 0 if gene is in deleted genes list
        else
            gvalue{j, 2} = 1;  % Mark as 1 otherwise
        end
    end
    % Run the GRPRchecker function
    [GR_2, PR_2] = GRPRchecker(ori_model, test_target, gvalue);

    % Save results to the text file
    fprintf(fid, '--- DBgDel Results ---\n');
    fprintf(fid, 'Target: %s\n', test_target);
    fprintf(fid, 'Time cost: %.4f seconds\n', overallElapsedTime);
    fprintf(fid, 'GR: %.4f\n', GR_1);
    fprintf(fid, 'PR: %.4f\n', PR_1);
    
    fprintf(fid, '\n--- gMCSE Results ---\n');
    fprintf(fid, 'Target: %s\n', test_target);
    fprintf(fid, 'Time cost: %.4f seconds\n', time_cost_gmcse);
    fprintf(fid, 'GR: %.4f\n', GR_2);
    fprintf(fid, 'PR: %.4f\n', PR_2);

    % Close the file after writing
    fclose(fid);

    % Display confirmation message
    disp(['Results have been saved to ', output_file]);
end