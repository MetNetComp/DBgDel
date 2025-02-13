function RC_gene_selector(model, RC_size)
    % Function to randomly select a specified number of genes from a metabolic model
    % Inputs:
    %   - model: path to the model file (e.g., 'e_coli_core.mat')
    %   - RC_size: number of genes to randomly select

    % Randomly select genes
    rng('shuffle');  % Shuffle the random seed
    random_genes = model.genes(randperm(numel(model.genes), RC_size));

    % Display selected genes
    disp('Randomly selected genes:');
    disp(random_genes);

    % Convert to a table and save as CSV
    random_genes_table = table(random_genes, 'VariableNames', {'Remaining_gene'});
    output_file = fullfile('Ablation_Test_and_Baselines/', ['RC_genes_' inputname(1) '.csv']);
    writetable(random_genes_table, output_file);

    % Confirm save
    disp(['RC genes saved to ', output_file]);
end