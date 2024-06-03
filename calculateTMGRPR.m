function [] = calculateTMGRPR(model2)

options=cplexoptimset('cplex');
options.mip.tolerances.integrality=10^(-12);

m0=size(model2.mets,1);
for i=1:m0
    i
    model=model2;
    targetMet=model.mets{i}
    [model,targetRID,extype] = modelSetting(model,targetMet)
    m=size(model.mets,1);
    gid=find(model.c);
    pid=targetRID;
    
    [opt1.x, opt1.f, opt1.stat, opt1.output] = ...
        cplexlp(-model.c, [],[], model.S, zeros(m,1),model.lb, model.ub);
    list0(i,1)=-opt1.f;
    
    model.c(gid)=0;
    model.c(pid)=1;
    [opt2.x, opt2.f, opt2.stat, opt2.output] = ...
        cplexlp(-model.c, [],[], model.S, zeros(m,1),model.lb, model.ub);
    if opt2.stat>0
        list0(i,2)=-opt2.f;
    else
        list0(i,2)=0;
    end
    model.lb(gid)=-opt1.f;
    model.ub(gid)=-opt1.f;
    
    [opt3.x, opt3.f, opt3.stat, opt3.output] = ...
        cplexlp(model.c, [],[], model.S, zeros(m,1),model.lb, model.ub);
    
    
    if opt3.stat>0
        list0(i,3)=opt3.f;
    else
        list0(i,3)=0;
    end
  
end

save('calculateTMGRPR.mat');
return;
end

