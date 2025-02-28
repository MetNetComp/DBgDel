function [gvalue gr pr it success] = DBgDel(model,targetMet,maxLoop,PRLB,GRLB,initial_remaining_gene_pool)
% DBgDel is a framework that determines gene deletion strategies 
%by integrating the extracted initial remaining gene pool from the database with 
%downstream algorithms to narrow the algorithmic search space.
%
%function [gvalue gr pr it success]  
%   = DBgDel(model,targetMet,maxLoop,PRLB,GRLB,initial_remaining_gene_pool)
%
%INPUTS
% model     COBRA model structure containing the following required fields to perform DBgDel.
%   rxns                    Rxns in the model
%   mets                    Metabolites in the model
%   genes               Genes in the model
%   grRules            Gene-protein-reaction relations in the model
%   S                       Stoichiometric matrix (sparse)
%   b                       RHS of Sv = b (usually zeros)
%   c                       Objective coefficients
%   lb                      Lower bounds for fluxes
%   ub                      Upper bounds for fluxes
%   rev                     Reversibility of fluxes
%
% targetMet   target metabolites
%             (e.g.,  'btn_c')
% maxLoop   the maximum number of iterations in DBgDel
% PRLB           the minimum required production rates of the target metabolites
%                     when DBgDel searches the gene deletion
%                     strategy candidates. 
%                     (But it is not ensured to achieve this minimum required value
%                      when GR is maximized withoug PRLB.)
% GRLB         the minimum required growth rate 
%                     when DBgDel searches the gene deletion
%                     strategy candidates. 
%
%OUTPUTS
% gvalue     The first column is the list of genes in the original model.
%            The second column contains a 0/1 vector indicating which genes should be deleted.
%             0 indicates genes to be deleted.
%             1 indecates genes to be remained.
% gr         the growth rate obained when the gene deletion strategy is
%              applied and the growth rate is maximized.
% pr         the target metabolite production rate obained 
%              when the gene deletion strategy is applied and the growth rate is maximized.
% it       indicates how many iterations were necessary to obtain the
%              solution.
% success  indicates whether DBgDel obained an appropriate gene
%               deletion strategy. (1:success, 0:failure)
%
%   May. 21, 2024   Ziwei YANG, Takeyuki TAMURA

tic;
ori_model=model;
n=size(model.rxns,1);
for i=1:n
   if model.ub(i)>9999
       model.ub(i)=1000;
   end
   if model.lb(i)<-9999
       model.lb(i)=-1000;
   end
end

gvalue=[];
gr=-1;pr=-1;it=0;success=0;
%big=10;
%big=1;
options=cplexoptimset('cplex');
options.mip.tolerances.integrality=10^(-12);
%options=cplexoptimset('TolXInteger',10^(-12));

sss=sprintf('DBgDel%d.mat');

[ori_model,targetRID,extype] = modelSetting(ori_model,targetMet)
[model,targetRID,extype] = modelSetting(model,targetMet)

m=size(model.mets,1);
n=size(model.rxns,1);
g=size(model.genes,1);
vg=zeros(g,1);
gid=find(model.c);
pid=targetRID;
model2=model;
model2.c(gid)=0;
model2.c(targetRID)=1;

[optPre.x, optPre.f, optPre.stat, optPre.output] = ...
        cplexlp(-model2.c, [],[], model2.S, zeros(m,1),model2.lb, model2.ub);
if optPre.stat<=0 
    display('no solution 1')
    save(sss);
    return;
elseif -optPre.f < PRLB
    display('TMPR < PRLB')
    save(sss);
    return;
end
    
model2=model;
model.lb(pid)=PRLB;
model.lb(gid)=GRLB;

[opt0.x, opt0.f, opt0.stat, opt0.output] = ...
    cplexlp(-model.c, [],[], model.S, zeros(m,1),model.lb, model.ub);
TMGR=-opt0.f;
big=TMGR; 
[term,ng,nt,nr,nko,reactionKO,reactionKO2term] = readGeneRules(model);
 [f,intcon,A,b,Aeq,beq,lb,ub,xname] = geneReactionMILP_DB(model,term,ng,nt,nr,nko,reactionKO,initial_remaining_gene_pool);
 
 lp.Aeq=[model.S zeros(m,ng+nt+nko);
               zeros(size(Aeq,1),n) Aeq];
lp.beq=[zeros(m,1); zeros(size(Aeq,1),1)];
j=1;
for i=1:size(model.grRules,1)
    if isempty(model.grRules{i,:})==0
        ind(1,j)=i;
        j=j+1;
    end
end
z1=-diag(model.ub);
z2=diag(model.lb);
z3=eye(n);

lp.A=[zeros(size(A,1),n) A;
           z3(ind,:) zeros(size(ind,2),ng+nt) z1(ind,ind);
               -z3(ind,:) zeros(size(ind,2),ng+nt) z2(ind,ind)];
lp.b=[b; zeros(2*nko,1)];
lp.lb=[model.lb; lb];
lp.ub=[model.ub; ub];
lp.f=[-model.c; zeros(ng+nt,1); big*ones(nko,1)];
for i=1:n
    s1=repelem('C',n);
    s2=repelem('B',ng+nt+nko);
    lp.ctype=sprintf('%s%s',s1,s2);
end
A2=lp.A;
b2=lp.b;

it =1;
while it<=maxLoop
    it
    gr=-1;pr=-1;
    [opt.x, opt.f, opt.stat, opt.output] = ...
        cplexmilp(lp.f, lp.A, lp.b, lp.Aeq, lp.beq,[],[],[],lp.lb, lp.ub,lp.ctype,[],options);

    if opt.stat>0
        for i=1:n
            vx(i,it)=opt.x(i);
            result{i,it+1}=opt.x(i);
            result{i,1}=model.rxns{i};
        end
        for i=1:ng
            vg(i,it)=opt.x(n+i);
            result{n+i,1}=xname{i};
            result{n+i,it+1}=vg(i,it);
            gvalue{i,1}=xname{i};
            gvalue{i,2}=vg(i,it)>0.1;
        end
        for i=1:nt
            vt(i,1)=opt.x(n+ng+i);
            result{n+ng+i,1}=xname{ng+i};
            result{n+ng+i,it+1}=opt.x(n+ng+i);
        end
        for i=1:nko
            vko(i,it)=opt.x(n+ng+nt+i);
            result{n+ng+nt+i,1}=xname{ng+nt+i};
            result{n+ng+nt+i,it+1}=opt.x(n+ng+nt+i);
        end
    %save('a.mat');return;
    else
        system('rm -f clone*.log');
        if it==1
            display('no solution 2')
            save('a.mat');
            return;
        end
        display('no more candidates')
        system('rm -f clone*.log');
        save(sss);
        return;
    end
    [grRules] = calculateGR(ori_model,gvalue);
    
    lb2=ori_model.lb;
    ub2=ori_model.ub;
    for i=1:nr
        if grRules{i,4}==0
            lb2(i)=0;
            ub2(i)=0;
        end
    end
    [opt2.x, opt2.f, opt2.stat, opt2.output] = ...
        cplexlp(-ori_model.c, [],[], model.S, zeros(m,1),lb2, ub2);
    %save('a.mat');system('rm -f clone*.log');return;
    grprList(it,:)=[opt2.x(gid) opt2.x(pid)];
    gr=opt2.x(gid); pr=opt2.x(pid);
    result2(:,it)=opt2.x;
    
    if (opt2.x(gid)>=GRLB) &&  (opt2.x(pid)>=PRLB)
        [opt2.x(gid) opt2.x(pid)]
        vg(:,it)
        success=1;
        time=toc;
        system('rm -f clone*.log');
        save(sss);
        return;
    end
    
    zeroList(:,it)=vko(:,it)<0.01;
    dA=[zeros(1,nr+ng+nt) -(zeroList(:,it))'];
    db=-1;
    lp.A=[lp.A; dA];
    lp.b=[lp.b; db];

    it=it+1;
    system('rm -f clone*.log');
    save(sss);
end
vg=vg>0.1;
save(sss);
end

