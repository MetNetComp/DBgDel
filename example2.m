function [outputArg1,outputArg2] = example2()
%example2 calculates the gene deletion strategy for growth coupling
%for riboflavin in iML1515.
%

load('CBM_model/iML1515.mat');
model=iML1515;
initial_remaining_gene_pool = readtable('initial_remaining_gene/iML1515_ex.csv', 'Delimiter', ',').Remaining_gene.';
% Start the timer
tic;
[gvalue,gr,pr,it,success]=DBgDel(model,'ribflv_c',10,0.1,0.1,initial_remaining_gene_pool);

% Stop the timer and display the elapsed time
elapsedTime = toc;
fprintf('Elapsed Time for example2: %.2f seconds\n', elapsedTime);
[GR,PR]=GRPRchecker(model,'ribflv_c',gvalue)

save('example2.mat');
end

