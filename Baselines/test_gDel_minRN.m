%CobraToolbox&Solver
initCobraToolbox(false)
changeCobraSolver('ibm_cplex', 'all', 1);

%load a test model, e_coli_core/iMM904/iML1515
load('e_coli_core.mat');
model = e_coli_core;

%test on all possible target in above model
target_names = readtable('possible_target_e_coli_core_no_e.csv', 'Delimiter', ',');
% Start the overall timer for the loop
overallTimer = tic;

% Initialize the cell array to accumulate all results
allResults = cell(numel(target_names), 6); 

% Loop through each target
for i = 1:numel(target_names)
    target = target_names{i,:}{1};
    
    % Start the timer for gDel_minRN
    tic;
    [gvalue, gr, pr, it, success] = gDel_minRN(model, target, 10, 0.1, 0.1);
    % Stop the timer for gDel_minRN and display the elapsed time
    elapsedTime = toc;
    
    % Save the gDel_minRN results into the allResults array
    allResults{i, 1} = target;
    allResults{i, 2} = elapsedTime;         % Elapsed Time 
    allResults{i, 3} = gr;                  % gr 
    allResults{i, 4} = pr;                  % pr 
    allResults{i, 5} = it;                  % it 
    allResults{i, 6} = success;             % success=1
end

% Stop the overall timer for the loop and display the elapsed time
overallElapsedTime = toc(overallTimer);
disp(['Overall Execution Time: ', num2str(overallElapsedTime), ' seconds']);

% Convert the allResults cell array into a table for final CSV export
finalResultsTable = cell2table(allResults, 'VariableNames', {'Target', 'Elapsed Time', 'gr', 'pr', 'it', 'success'});

% Calculate the success rate (percentage of successful computations)
successRate = sum([allResults{:, 6}] == 1) / numel(allResults) * 100;
disp(['Overall Success Rate: ', num2str(successRate), '%']);

% Save the final results into a single CSV file
writetable(finalResultsTable, 'result_gDelminRN_e_coli_core.csv', 'WriteRowNames', false);