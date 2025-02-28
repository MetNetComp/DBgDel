function [outputArg1,outputArg2] = example4()
%example4 calculates the gene deletion strategy for growth coupling
%for riboflavin in iML1515, using RC genes as the initial remaining gene set.

load('CBM_model/iML1515.mat');
model=iML1515;
%The RC gene set size is defaulted to be 267 in iML1515.
RC_gene_selector(iML1515,267);

initial_remaining_gene_pool = readtable('initial_remaining_gene/RC_genes_iML1515.csv', 'Delimiter', ',').Remaining_gene.';
% Start the timer
tic;
[gvalue,gr,pr,it,success]=DBgDel(model,'ribflv_c',50,0.1,0.1,initial_remaining_gene_pool);

% Stop the timer and display the elapsed time
elapsedTime = toc;
fprintf('Elapsed Time for RC_example4: %.2f seconds\n', elapsedTime);
[GR,PR]=GRPRchecker(model,'ribflv_c',gvalue)

save('RC_example4.mat');
end