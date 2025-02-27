% % %Download model files in mat first if necessary
% % modelFile=sprintf('%s.mat',modelName);
% % downloadLink=sprintf('http://bigg.ucsd.edu/static/models/%s',modelFile);
% % websave('e_coli_core.mat',downloadLink);
% % load(modelFile);
% % loadModel=sprintf('model=%s',modelName);
% % eval(loadModel);

%Prepare CobraToolbox[1]
% [1] Heirendt L, Arreckx S, Pfau T, et al. Creation and analysis of biochemical constraint-based models using 
% the COBRA Toolbox v. 3.0[J]. Nature protocols, 2019, 14(3): 639-702.
initCobraToolbox;
%Use'ibm_cplex' or 'gurobi' Solver
changeCobraSolver('ibm_cplex', 'all', 1);


%Model name, pls change to test on iML1515,iMM904 and other models.
load('e_coli_core.mat');
model = e_coli_core;

%Setup parameters and input
id_biomass=find(model.c);
id_glucose=find(strcmp('EX_glc__D_e',model.rxns));
id_oxygen=find(strcmp('EX_o2_e',model.rxns));
met_oxygen=find(strcmp('o2_e',model.mets));
met_glucose=find(strcmp('glc__D_e',model.mets));
oxygen_bound=10;
glucose_bound=15;
model.lb(id_oxygen)=-1*oxygen_bound;
model.lb(id_glucose)=-1*glucose_bound;
model.ub(id_glucose)=0;
model.ub(id_oxygen)=0;


% Baseline test on GDLS [2]
% [2] Lun D S, Rockwell G, Guido N J, et al. Large‐scale identification of genetic design 
% strategies using local search[J]. molecular systems biology, 2009, 5(1): 296.
count=0;
nMet=size(model.mets,1);
stat=cell(nMet,5);
environment=getEnvironment();
%CoreUsesageSet
%parpool=(64);
for i=1:nMet
    if i~=met_oxygen && i~=met_glucose
        exid=find(strcmp(model.rxns,sprintf('EX_%s',model.mets{i})));
        dxid=find(strcmp(model.rxns,sprintf('DX_%s',model.mets{i})));
        if ~isempty(exid)
           id_target=exid;
           newModel=model;
        elseif ~isempty(dxid)
           id_target=dxid;
           newModel=model;
        else
           newModel=addExchangeRxn(model,model.mets{i});
           id_target=size(newModel.rxns,1);
           newModel.S(i,id_target)=-1;
        end
        model1=newModel;
        model1.c(id_biomass)=0;
        model1.c(id_target)=1;
        restoreEnvironment(environment);
        changeCobraSolver('ibm_cplex','LP',0,-1);
        opt=optimizeCbModel(model1);
        if opt.f>0
            restoreEnvironment(environment);
            changeCobraSolver('ibm_cplex','MILP',0,-1);
            try
                [gdlsSolution,bilevelMILPProblem,gdlsSolutionStructs]=GDLS(newModel,newModel.rxns{id_target},'timeLimit',600,'minGrowth',0.01); 
                stat(i,:)={'TMGR',gdlsSolution,bilevelMILPProblem,gdlsSolutionStructs,t};
                if gdlsSolution.minTargetProd>0
                    count=count+1;
                end
            catch MException
                stat(i,:)={'TMGR',[],[],[],MException};
            end
        end
    end
end
delete(gcp('nocreate'));
filename=sprintf('%s_GDLS_%s',model.description,datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
save(filename);

% Baseline test on optGene [3]
% [3] Rocha I, Maia P, Rocha M, et al. OptGene: a framework for in silico metabolic engineering[C]//10th International 
% Conference on Chemical and Biological Engineering. Portugal: University of Minho, 2008: 218-219.
count=0;
nMet=size(model.mets,1);
stat=cell(nMet,6);
environment=getEnvironment();
%parpool=(64);
for i=1:nMet
    if i~=met_oxygen && i~=met_glucose
        exid=find(strcmp(model.rxns,sprintf('EX_%s',model.mets{i})));
        dxid=find(strcmp(model.rxns,sprintf('DX_%s',model.mets{i})));
        if ~isempty(exid)
           id_target=exid;
           newModel=model;
        elseif ~isempty(dxid)
           id_target=dxid;
           newModel=model;
        else
           newModel=addExchangeRxn(model,model.mets{i});
           id_target=size(newModel.rxns,1);
           newModel.S(i,id_target)=-1;
        end
        model1=newModel;
        model1.c(id_biomass)=0;
        model1.c(id_target)=1;
        restoreEnvironment(environment);
        changeCobraSolver('ibm_cplex','LP',0,-1);
        opt=optimizeCbModel(model1);
        if opt.f>0
            restoreEnvironment(environment);
            changeCobraSolver('ibm_cplex','MILP',0,-1);
            try
                [x,population,scores,optGeneSol]=optGene(newModel,newModel.rxns{id_target},newModel.rxns{id_glucose},newModel.genes,'TimeLimit',600); 
                stat(i,:)={'TMGR',x,population,scores,optGeneSol,t};
                if optGeneSol.obj>0
                    count=count+1;
                end
            catch MException
                stat(i,:)={'TMGR',[],[],[],[],MException};
            end
        end
    end
end
delete(gcp('nocreate'));
filename=sprintf('%s_optGene_%s',model.description,datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd'));
save(filename);