function [outputArg1,outputArg2] = example1()
%example1 calculates the gene deletion strategy for growth coupling
%for pantothenate in iMM904.

load('CBM_model/iMM904.mat');
model=iMM904;
initial_remaining_gene_pool = readtable('initial_remaining_gene/iMM904_ex.csv', 'Delimiter', ',').Remaining_gene.';
% Start the timer
tic;
[gvalue,gr,pr,it,success]=DBgDel(model,'pnto__R_c',10,0.1,0.1,initial_remaining_gene_pool);

% Stop the timer and display the elapsed time
elapsedTime = toc;
fprintf('Elapsed Time for example1: %.2f seconds\n', elapsedTime);
[GR,PR]=GRPRchecker(model,'pnto__R_c',gvalue)

save('example1.mat');
end