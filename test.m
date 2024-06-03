function [outputArg1,outputArg2] = test()
%test example calculates the gene deletion strategy for growth coupling
%for succ_e in e_coli_core.
%
%
initCobraToolbox(false);
load('CBM_model/e_coli_core.mat');
initial_remaining_gene_pool = readtable('initial_remaining_gene/e_coli_core_ex.csv', 'Delimiter', ',').Remaining_gene.';
model=e_coli_core;
[gvalue,gr,pr,it,success]=DBgDel(model,'succ_e',10,0.001,0.001,initial_remaining_gene_pool)
[GR,PR]=GRPRchecker(model,'succ_e',gvalue)
save('test.mat');
end