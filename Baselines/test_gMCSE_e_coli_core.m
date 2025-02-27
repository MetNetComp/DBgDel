%Step1: gMCSE

% CobraToolbox and CellNetAnalyzer
initCobraToolbox(false);  
startcna;
% Set Java heap memory for EFMTool before computation
% setJavaHeapSize(32); 
% Adjust MATLAB parpool based on available memory if necessary
% parpool('local', 12);

load('e_coli_core.mat');
ori_model = e_coli_core;
% Read the target CSV
Targetlist = readtable('possible_target_e_coli_core_no_e.csv', 'Delimiter', ',');
% Extract the first column as a cell array (useful for strings)
productionList = Targetlist{:, 1}; 
% Read the e_coli_core_RID file (gene-pseudoreaction-ID)
map_gene_list = readtable('e_coli_core_RID.csv', 'Delimiter', ',');
D_list = map_gene_list{:, 1};      
R_names = map_gene_list{:, 2}; 

growth_idx = find(ori_model.c); %Growth reaction ID
GR_threshold = 0.001; % Minimum growth threshold
PR_threshold = 0.001; % Maximum production threshold
maxSolutions = 1; % # Solution needed for records
rkoCost = [];
rkiCost= [];
maxCost = [];
gkoCost = [];
gkiCost = []; 
options.milp_time_limit = 300;
verbose = [];

% Initialize an empty cell array to store the results
results = {};
% Initialize variables to accumulate the total time cost and success count
total_time = 0;
success_count = 0;

for i = 1:length(productionList)
    % Display the current target
    targetMet = char(productionList(i)); 
    disp([num2str(i), ' Target: ', targetMet]);

    [model, targetRID, extype] = modelSetting(ori_model, targetMet);
    cnap = CNAcobra2cna(model);
    [cnap, enzymes, genes, gpr_rules] = CNAgenerateGPRrules(cnap);
    [T, t] = initializeTargetRegion(cnap.numr, growth_idx, targetRID, GR_threshold, PR_threshold);
    [D, d] = initializeDesiredRegion(cnap.numr, growth_idx, GR_threshold);

    % Start timer
    tic;

    % Call the CNAgeneMCSEnumerator2 function
    [rmcs, full_mcs, full_cnap, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, status, obj] = ...
        CNAgeneMCSEnumerator2(cnap, T, t, D, d, rkoCost, rkiCost, ...
                              maxSolutions, maxCost, gkoCost, gkiCost, gpr_rules, options, verbose);

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
        % Increase success count
        success_count = success_count + 1;
    end

    % End timer
    time_cost = toc;

    % Accumulate the total time cost
    total_time = total_time + time_cost;

    % Store results
    results{end+1, 1} = targetMet;      % Target metabolite
    results{end, 2} = time_cost;        % Time cost
    results{end, 3} = status;           % Status
    results{end, 4} = deleted_gene_str; % Deleted genes as string (or empty if status ≠ 0/3)

    % Display results in the loop
    disp(['Time cost for target ', targetMet, ': ', num2str(time_cost), ' seconds']);
    disp(['Status: ', num2str(status)]);
    disp(['Deleted genes: ', deleted_gene_str]);
end

% Calculate the average time
avg_time = total_time / length(productionList);

% Calculate the success rate
success_rate = success_count / length(productionList);

% Display the average time and success rate
disp(['Average time per iteration: ', num2str(avg_time), ' seconds']);
disp(['Success rate: ', num2str(success_rate * 100), '%']);

% Add the average time and success rate as the last row in the results
results{end+1, 1} = 'Average Time';  % Label for the average time row
results{end, 2} = avg_time;  % Add average time
results{end, 3} = success_rate;  % Add success rate 

% Convert the results to a table for easy saving to CSV
results_table = cell2table(results, 'VariableNames', {'Target', 'TimeCost', 'Status', 'DeletedGenes'});

% Write the table to a CSV file
writetable(results_table, 'e_coli_core_results.csv');

%Step2: Check the PR/GR

%Get deleted genes
geneD = results_table.DeletedGenes; 

% Initialize the output cell arrays
GR_all = [];  % Store GR values for each target
PR_all = [];  % Store PR values for each target
success_all = [];  % Store Success (1 or 0) for each target

% Loop through each target in the CSV
for i = 1:height(results_table)-1
    targetMet = results_table.Target{i};  % Get target name
    deleted_genes = strsplit(geneD{i}, ', ');  % Split the deleted genes into a cell array
    
    % Store gene names and determine if they are in the deleted genes list
    for j = 1:length(ori_model.genes)
        gvalue{j, 1} = ori_model.genes{j};  % Store gene name
        if ismember(ori_model.genes{j}, deleted_genes)
            gvalue{j, 2} = 0;  % Mark as 0 if gene is in deleted genes list
        else
            gvalue{j, 2} = 1;  % Mark as 1 otherwise
        end
    end
    
    % Run the GRPRchecker function
    [GR, PR] = GRPRchecker(ori_model, results_table.Target{i}, gvalue);
    
    % Store the results for each target
    GR_all = [GR_all; GR];
    PR_all = [PR_all; PR];
    
    % Calculate success: 1 if both GR and PR >= threshold, else 0
    if GR >= GR_threshold && PR >= PR_threshold
        success_all = [success_all; 1];
    else
        success_all = [success_all; 0];
    end
end

% Add GR, PR, and Success columns to the table
results_table.GR = [GR_all; NaN];  % Append empty row for GR
results_table.PR = [PR_all; NaN];  % Append empty row for PR
results_table.Success = [success_all; NaN];  % Append empty row for Success

% Calculate and append the success rate as the summary row
success_rate = sum(success_all) / length(success_all);
results_table{end, end} = success_rate;
% Save the updated table to a new CSV
disp(['Success rate: ', num2str(success_rate)]);
writetable(results_table, 'e_coli_core_results_check.csv'); 