function [outputArg1,outputArg2] = RC_example1()
%RC_example1 calculates the gene deletion strategy for growth coupling for
%pantothenate in iMM904, using RC genes as the initial remaining gene set.

load('CBM_model/iMM904.mat');
model=iMM904;
%The RC gene set size is defaulted to be 195 in iMM904.
RC_gene_selector(iMM904,195);

initial_remaining_gene_pool = readtable('initial_remaining_gene/RC_genes_iMM904.csv', 'Delimiter', ',').Remaining_gene.';
% Start the timer
tic;
[gvalue,gr,pr,it,success]=DBgDel(model,'pnto__R_c',50,0.1,0.1,initial_remaining_gene_pool);

% Stop the timer and display the elapsed time
elapsedTime = toc;
fprintf('Elapsed Time for RC_example1: %.2f seconds\n', elapsedTime);
[GR,PR]=GRPRchecker(model,'pnto__R_c',gvalue)

save('RC_example1.mat');
end