function [GR ,PR] = GRPRchecker(model,targetMet,givenGvalue)

options=cplexoptimset('cplex');
options.mip.tolerances.integrality=10^(-12);

[model,targetRID,extype] = modelSetting(model,targetMet)

m=size(model.mets,1);
n=size(model.rxns,1);
g=size(model.genes,1);
gid=find(model.c);
pid=targetRID;


model2=model;
[grRules0] = calculateGR(model,givenGvalue);
lb2=model.lb;
ub2=model.ub;

for i=1:n
    if grRules0{i,4}==0
        lb2(i)=0;
        ub2(i)=0;
    end
end
[opt0.x, opt0.f, opt0.stat, opt0.output] = ...
    cplexlp(-model.c, [],[], model.S, zeros(m,1),lb2, ub2);
[opt0.x(gid) opt0.x(pid)]

GR0=-opt0.f;
lb2(gid)=GR0;
ub2(gid)=GR0;
model2.c(gid)=0;
model2.c(pid)=1;
[opt1.x, opt1.f, opt1.stat, opt1.output] = ...
    cplexlp(model2.c, [],[], model.S, zeros(m,1),lb2, ub2);
GR=GR0
PR=opt1.x(pid)

save('GRPRchecker.mat');
return;
end

