function [rmcs, full_mcs, full_cnap, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, status, obj] = ...
    CNAgeneMCSEnumerator2(cnap,T,t,D,d,rkoCost,rkiCost,maxSolutions,maxCost,gkoCost,gkiCost,gpr_rules,options,verbose)
%
% ------------------------------------------------
% CellNetAnalyzer API function 'CNAgeneMCSEnumerator2'
% ------------------------------------------------
% --> Minimal cut set enumeration or (iterative) search with multiple target and 
%     desired regions, gene or reaction knockouts, additions and individual
%     intervention costs.
%
% Usage: [rmcs, full_mcs, full_cnap, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, status, obj] = ...
%          CNAgeneMCSEnumerator2(cnap,T,t,D,d,koCost,kiCost,maxSolutions,maxCost,gkoCost,gkiCost,gpr_rules,options,verbose)
% 
% Computes gene-based minimal cut sets (MCS) for metabolic strain design 
% from a constraint-based model and a number of design specifications. 
% The MCS approach can be used to find interventions that turn wild type
% microbes into production hosts. In contrast to bi-level optimization
% techniques, the design specifications don't state an optimization goal
% but denote the desired (e.g. growth) and undesired flux-states (e.g.
% inferior product yield). The MCS algorithm will find interventions that
% shape the solution space of steady state flux-states according to those
% specifications.
% The design specifications are presented in the way of one or more target 
% "regions" (defined by one or multiple matrices T and vectors t) that denote 
% the subspace of undesired steady flux states which should be rendered
% inaccessible by applying a minimal cut set and optionally by a number of 
% desired "regions" (defined by one or multiple matrices D and vectors d) that 
% describe the behavior that should be maintained. MCS can consist of gene
% or reaction additions or deleted, each intervention attributed with an
% individual cost factor. If only Target and no Desired constraints are
% assigned, the entire model might be rendered infeasible by deletions from
% a minimal cut set.
%
% This function is a wrapper to the CNAMCSEnumerator2 that extends the metabolic
% network with gene-protein-reaction rules before the actual MCS computation. 
% Those rules are compressed prior to their integration as pseudoreactions
% and pseudometabolites into the network.
% 
% Example, how desired behavior can be expressed as a region:
%                             D          * r <=  d
%                   ( 0 -1 0 0 0 0 0 0 ) * r <= -0.1
%            <=>                        -r_2 <= -0.1
%            <=>                         r_2 >=  0.1
%    (also: More than one row is allowed for the definition of each Desired/Target region)
%
% The MCS ensure that all target flux vectors r obeying
%
%       cnap.stoichMat * r = 0
%       cnap.reacMin(r) <= r <= cnap.reacMax(r)
%       T * r <= t
%
% will be removed from the solution space and thus the undesired behaviors
% are blocked. 
%       Target{1}: r_prod <= 1 will, for example, lead to an elimination of all
%                              flux states with production rates inferior to 1.
% 
% Whereas (one or more) flux vectors r fulfilling
%
%       cnap.stoichMat * r = 0
%       cnap.reacMin(r) <= r <= cnap.reacMax(r)
%       D * r <= d
%
% must be kept in the solution space, maintaining desired
% functions/behaviors.
%       Desired{1}: r_growth >= 0.1 will, for example, maintain the capability to
%                                   grow with rates higher than 0.1.
%
% ------------------------------------------------
% Input:
% cnap        : a CellNetAnalyzer mass flow project
% T           : <cell>[numT x 1] Matrices specifying (with vector t) the target flux
%                regions as given above. T = {T1,T2, ...}. T1: <double>[rowsT1 x cnap.numr]
% t           : <cell>[numT x 1] Vectors specifying (with matrix T) the target flux
%                regions as given above. t = {t1,t2, ...}. t1: <double>[rowsT1 x 1]
% D           : <cell>[numD x 1] Matrices specifying (with vector d) the desired flux
%                regions as given above. D = {D1,D2, ...}. D1 <double>[rowsD1 x cnap.numr]
% d           : <cell>[numD x 1] Vectors specifying (with matrix D) the desired flux
%                regions as given above. d = {d1,d2, ...}. d1 <double>[rowsD1 x 1]
% koCost      : <double>[cnap.numr x 1] vector that indicates REACTIONS that 
%                can be knocked out. Notknockable reactions carry NaN. Knockable
%                reactions carry a value that represents their knockout cost. By
%                default, no reaction is knockable, however all genes are. If
%                gkoCost is not provided, the function will derive the knockable
%                genes and reactions from the koCost vector, also keeping non-gene
%                associated reactions knockable.
% kiCost      : <double>[cnap.numr x 1] vector that indicates REACTIONS that 
%               can be added or "knocked in". Notknockable reactions carry NaN. 
%               Addable reactions carry a value that represents their addition cost.
%               If a reaction is set knockable and addable at the same time. It
%               will be assumed that the reaction is addable. By default, no
%               reactions are addable. If gkiCost vector is not provided, addable
%               reactions will be anticipated from kiCost in the same way as for
%               koCost. KIs disable KOs.
% maxSolutions:<double> max. number of solutions. Aborts enumeration after 
%               this number of (compressed) solutions have been found. Inf 
%               is allowed (default).
% maxCost     : <double> Upper cost limit for the sum of KO/KI costs in one
%                MCS. If all interventions have a cost of 1 this value 
%                represents the total number of interventions. Inf is
%                allowed (default) but not advised in networks with 100+ 
%                reactions. For full enumerations in genome scale networks 
%                this is typically between 5 and 8. In an iterative search
%                this number is not limited but would still be set to a
%                value that seems feasible in practice (e.g. is the 
%                implementation of 15 gene interventions manageable?).
% gkoCost     : <double>[num_genes x 1] vector that indicates GENES that 
%                can be knocked out. Notknockable genes carry NaN. Knockable
%                genes carry a value that represents their knockout cost. By
%                default all genes are knockable.
% gkiCost     : <double>[num_genes x 1] vector that indicates GENES that 
%                can be added or "knocked in". Notknockable genes carry NaN. 
%                Addable genes carry a value that represents their addition cost.
%                If a genes is set knockable and addable at the same time. It
%                will be assumed that the genes is addable.
% gpr_rules   : Structure that contains GPR rules. See CNAgenerateGPRrules
%                function for details. Necessary when GPR rules are not provided 
%                in the 'generic reaction data' field 'geneProductAssociation'.
% options     : <struct> controls different pre- and post-processing steps and. 
%                MILP parameters. Set the according fields if you want to change 
%                the default settings.
%       options.mcs_search_mode              : <double> 
%                                                1: will iteratively search for single solutions. 
%                                                   This is preferred when a full enumeration of intervention
%                                                   strategies is not necessary or too time consuming. Consider
%                                                   this setting also if you chose a high value for maxSumItvCost.
%                                                2: (default)  will do a full bottom-up stepwise enumeration of MCS.
%       options.pre_GPR_network_compression  : <logical> Network compression (prior to GPR rule integration).
%                                                        (default FALSE)
%       options.compression_GPR              : <logical> GPR rule compression (prior to GPR rule integration).
%                                                        (default TRUE)
%       options.preproc_compression          : <logical> Network compression (after GPR rule extension and before running the MILP).
%                                                        (default TRUE)
%       options.preproc_check_feas           : <logical> Check if all target and desired regions are feasible in the original model.
%                                                        (default TRUE) Deactivate if you are sure D and T are feasible and you want
%                                                        to reduce the computation time overhead.
%       options.preproc_D_violations         : <double> Preprocessing step, where minimal knockout combinations (with a size up to the given number)
%                                               are computed  that disrupt one or more desired regions (if provided). These knockout combinations
%                                               violating the desired regions are excluded as feasible combinations in the MCS-MILP. 
%                                               Ideally this speeds up the search for MCS. Depending on the case, it might yet 
%                                               slow down MCS computation.  (default 0)
%       options.postproc_verify_mcs          : <logical> Verify MCS after the computation to ensure that no false positives are returned.
%                                                        (default TRUE)
%       options.milp_solver                  : <char> 'matlab_cplex' (default), 'java_cplex', 'intlinprog'
%       options.milp_time_limit              : <double> in seconds. Inf is allowed (default).
%       options.milp_bigM                    : <logical or double> 
%                                                 DOUBLE: The provided value is used as big M.
%                                                 TRUE: M will be set to 110% of the highest reaction boundary. 
%                                                 FALSE or zero: indicator constraints are used (default).
%       options.milp_nullspace               : <logical> Parameter for MILP construction (default FALSE) supported only for matlab_cplex and java_cplex
%     == These MILP settings can only be used with the MATLAB-CPLEX API (see options.milp_solver) ==
%       options.milp_split_level             : <logical> Parameter for MILP construction (default TRUE) supported only for matlab_cplex
%       options.milp_reduce_constraints      : <logical> Parameter for MILP construction (default TRUE) supported only for matlab_cplex
%       options.milp_combined_z              : <logical> Parameter for MILP construction (default TRUE) supported only for matlab_cplex
%       options.milp_irrev_geq               : <logical> Parameter for MILP construction (default TRUE) supported only for matlab_cplex
% verbose      : <double> 1 yes (default), 0 no. Select 0 for a silent run.
%
%
% ------------------------------------------------
% Output:
% rmcs        : <double>[numMCS x cnap.numr] Reaction-representation of gene-MCS. 
%                Notation: 0: reaction was left unchanged; 1: reactions was added; 
%                -1: reaction was removed; NaN: reaction was addition candidate 
%                but not added. 
%                These MCS refer to the original reaction-based CNA project/network 
%                that was passed as a function parameter. rMCS are derived from
%                the gene- (and reaction) deletions and additions and are ment to
%                support an easier evaluation of the found gene-mcs with the original 
%                model. They describe the 'phenotypical' metabolic changes of 
%                each mutant towards the wild type. For the minimal sets of genetic 
%                interventions (gene-MCS) relevant for the experimental 
%                implementation, use the output parameters 'full_mcs', and 
%                'full_cnap'.
% full_mcs    : <double>[numMCS x full_cnap.numr] Gene-based MCS. Notation:
%                0: reaction was left unchanged; 1: reactions was added; 
%                -1: reaction was removed; NaN: reaction was addition candidate 
%                but not added.  These MCS refer to the "full" CNA project which 
%                is the uncompressed network with all genes and GPR rules added 
%                as pseudoreactions and metabolites. gene-KOs are therefore 
%                indicated at the position of the respective gene-pseudoreaction.
% full_cnap   : CNA mass-flow project that contains the original uncompressed
%                network with the integrated GPR rules. The MCS setup
%                expanded to the gpr-integrated model is provided in the fields 
%                cmp_cnap.mcs.*. The names of the knockouts and additions can be 
%                obtained by mapping a single MCS onto the vector of reaction names.
%                The remapped MCS setup is provided in the fields full_cnap.mcs.*  .
% cmp_mcs     : <double>[cmp_numMCS x cmp_cnap.numr] MCS that were found on the 
%                compressed setup.
% cmp_cnap    : compressed CNA mass-flow project used for MCS computation.
%                The MCS setup is provided in the fields cmp_cnap.mcs.*  .
% mcs_idx_cmp_full: <double>[1 x numMCS] A vector that maps the decompressed 
%                   mcs to the original compressed ones.
% status      : <double> 0 successful, 1 timeout without any solution, 
%                2 infeasible setup, 3 timeout but (some) solutions were found
% obj         : Object that contains the MILP setup.
%

%
% This file is part of CellNetAnalyzer. Please visit
% http://www.mpi-magdeburg.mpg.de/projects/cna/cna.html
% for more information and the latest version of CellNetAnalyzer.
%
% Copyright (C) 2000-2020 by Steffen Klamt and Axel von Kamp,
% Max Planck Institute for Dynamics of Complex Technical Systems, Magdeburg, Germany.
%
% Contributors are listed in CONTRIBUTORS.txt.
%
% This software can be used under the terms of our CellNetAnalyzer License.
% A copy of the license agreement is provided in the file named "LICENSE.txt"
% included with this software distribution. The license is also available online at
% http://www2.mpi-magdeburg.mpg.de/projects/cna/license.html
%
% For questions please contact: cellnetanalyzer@mpi-magdeburg.mpg.de

if isa(T,'double') && isa(t,'double')  % compatibility to old MCSEnumerator
    T = {T};
    t = {t};
end
if nargin < 5 || isempty(D)
    D = {};
    d = {};
elseif nargin >= 5 && ~isempty(D) && isa(D,'double') && isa(d,'double') % compatibility to old MCSEnumerator
    D = {D};
    d = {d};
end
if nargin < 13 || isempty(options)
    options = struct;
end
if isa(options,'struct')
    if ~isfield(options,'pre_GPR_network_compression')
        options.pre_GPR_network_compression = false;
    end
    if ~isfield(options,'compression_GPR')
        options.compression_GPR = true;
    end
    if ~isfield(options,'preproc_check_feas')
        options.preproc_check_feas = true;
    end
    if ~isfield(options,'milp_solver')
        options.milp_solver = 'matlab_cplex';
    elseif strcmp(options.milp_solver,'matlab_cplex') || strcmp(options.milp_solver,'java_cplex') || ...
            strcmp(options.milp_solver,'intlinprog') || strcmp(options.milp_solver,'glpk')
    else
        error(['Solver ' options.milp_solver ' is not supported']);
    end
else
    error('Function parameter ''options'' must be either a struct or empty.');
end
if nargin < 14 || isempty(verbose)
    verbose = 1;
end
if nargin < 7 ||isempty(rkiCost)
    rkiCost = nan(cnap.numr,1);
end
if nargin < 8 ||isempty(maxSolutions)
    maxSolutions = inf;
end
if nargin < 9 ||isempty(maxCost)
    maxCost = inf;
end

%% 0 Prepare Setup
% 0.1) Check if D and T are feasible in original model
if options.preproc_check_feas
    displ('Verifying D and T region feasibility.',verbose);
    feas_D = testRegionFeas(cnap,D,d,2);
    feas_T = testRegionFeas(cnap,T,t,2);
    if any(~feas_T)
        error(['At least one target region (T' num2str(find(~feas_T)) ') is infeasible in the original model.']);
    end
    if any(~feas_D)
        error(['At least one desired region (D' num2str(find(~feas_D)) ') is infeasible in the original model.']);
    end
end

displ('== gene MCS Computation ==',verbose);
if nargin < 12 || isempty(gpr_rules)
    if nargin >= 10 && all(isnan(gkoCost)) && all(isnan(gkiCost)) && ~isempty(gkoCost) % if entirely reaction based
        genes = {};
        gkoCost = [];
        gkiCost = [];
        gpr_rules = struct('reaction',{},'strReac',{},'genes',{},'strGene',{},'name',{});
    elseif ~isfield(cnap,'gpr_rules') || ~isfield(cnap,'genes')
        [~,~,genes,gpr_rules] = CNAgenerateGPRrules(cnap);
    else
        gpr_rules = cnap.gpr_rules;
        genes = cnap.genes;
    end
else
    [a,b] = unique([gpr_rules(:).genes]);
    gname = [gpr_rules(:).strGene]';
    genes(a) = gname(b);
end
if nargin < 6 || isempty(rkoCost) % if no reaction KO costs are defined, assume no reactions but all genes are knockable
    rkoCost = nan(1,cnap.numr);
    rkoCost([gpr_rules(:).reaction]) = 1;
end
if nargin < 10 || isempty(gkoCost) % derive from reaction knock-outs, if not defined. Set reactions with gene-association to notknockable.
    gkoCost = ones(1,length(genes));
    for i = 1:length(genes) % if all associated reactions are notknockable, so is the gene.
        rules_i = cellfun(@(x) ismember(i,x),{gpr_rules(:).genes});  % bool vector: Rules that contain the gene
        if all(isnan(rkoCost([gpr_rules(rules_i).reaction])))
            gkoCost(i) = nan;
        end
    end
else
    gkoCost = gkoCost(:)';
end
if nargin < 11 || isempty(gkiCost) % derive from reaction knock-ins, if not defined
    gkiCost = nan(1,length(genes));
    for i = 1:length(genes)
        rules_i = cellfun(@(x) ismember(i,x),{gpr_rules(:).genes});  % bool vector: Rules that contain the gene
        if ~all(isnan(rkiCost([gpr_rules(rules_i).reaction])))
            gkiCost(i) = min(rkiCost([gpr_rules(rules_i).reaction]));
        end
    end
else
    gkiCost = gkiCost(:)';
end
% 0.2) Prepare koCost, gkoCost, kiCost and gkiCost
gkoCost(~isnan(gkiCost)) = nan; % gene knock-ins 'override' gene knock-outs
c_macro = cnap.macroDefault;
% Transform ki & ko vector to row vector with 'nan' as notknockable and a double value as 'weight'
rkiCost = rkiCost(:)';
rkiCost([gpr_rules(:).reaction]) = nan;  % gene-associations 'override' reaction knock-outs
rkoCost = rkoCost(:)';
rkoCost([gpr_rules(:).reaction]) = nan;
rkoCost(~isnan(rkiCost)) = nan; % knock-ins 'override' knock-outs

% Making backups of some variables. Needed later for mcs expansion
full_koCost = rkoCost;
full_kiCost = rkiCost;
full_gkoCost = gkoCost;
full_gkiCost = gkiCost;
full_T = T;
full_D = D;
cnap_orig = cnap;
gpr_rules_orig = gpr_rules;

%% 1. Flux Variability Analysis
% 1.1) Determine reversibilities of metabolic model
displ('FVA to find blocked reactions and irreversibilities.',verbose);
[reacMin, reacMax] = CNAfluxVariability(cnap,[],c_macro,2,1:cnap.numr,[],[],0);
% All reaction bounds could be set here, but to avoid numerical issues, only
% exact zeros are set.
cnap.reacMin(reacMin == 0) = 0;
cnap.reacMax(reacMax == 0) = 0;
reac_off = find(reacMin == 0 & reacMax == 0);
displ(['#total_reacs: ' num2str(cnap.numr) ...
      ', #blocked: ' num2str(length(reac_off))],verbose);

% 1.2) Determine boundaries in desired scenarios and identify further essentials/notknockables
displ('FVA to determine essential reactions under desired conditions.',verbose);
essential = false(cnap.numr,1);
for i = 1:length(D)
    [lb_D, ub_D] = CNAfluxVariability(cnap,[],c_macro,2,1:cnap.numr,D{i},d{i},0); %#ok<*AGROW> Allow change of size in each loop iteration
    % An approximation is done to avoid marking essential reactions due
    % to numerical effects (e.g. lb: 1e-14 ub:5 is maybe not truely essential)
    essential = essential | ...
                (sign((abs(lb_D)>cnap.epsilon).*lb_D) .* ... % when lower bound sign equals
                 sign((abs(ub_D)>cnap.epsilon).*ub_D) == 1);  % upper bound sign
    % (the number of essential reacs changes at epsilon > 1e-9 / 1e-8)
    % *could also be used to select mandatory knock-ins
end
displ(['#essential: ' num2str(sum(essential))],verbose);
rkoCost(essential) = nan; % make essential reactions "notknockable"

%% 2. Reduce and lump set of GPR rules
if options.compression_GPR && ~isempty(gpr_rules)
    displ('Reducing set of gene protein reaction (GPR) rules.',verbose);
    % 2.1 remove rules for reactions that are off and continue with the reduced set
    gpr_rules = gpr_rules(~ismember([gpr_rules.reaction],reac_off));
    rule_removelist = zeros(1,length(gpr_rules));
    [reac_abund,reac] = hist([gpr_rules(:).reaction],unique([gpr_rules(:).reaction]));
    cgenes = {gpr_rules(:).genes}; % cell array of genes
    remaining_genes = intersect(cell2mat(cgenes),1:length(genes)); % all genes that weren't only catalyzing blocked reactions
    for i = remaining_genes
        % 2.2 if a gene only catalyzes essential reactions -> set gene notknockable
        rules_i = cellfun(@(x) ismember(i,x),cgenes); % bool vector: Rules that contain the gene
        if all(essential([gpr_rules(rules_i).reaction]))
            gkoCost(i) = nan;
        end
        % 2.3 if the gene is essential to one essential reaction -> set gene notknockable
        for j = unique([gpr_rules(rules_i).reaction])
                             % if all reaction-gene-rules for a particular reaction contain the gene
            if essential(j) && sum([gpr_rules(rules_i).reaction] == j) == reac_abund(j == reac)
                gkoCost(i) = nan;
            end
        end
        % 2.5 if gene is notknockable remove it from all rules
        if isnan(gkoCost(i)) && isnan(gkiCost(i)) && any(rules_i)
            for j = find(rules_i)
                gpr_rules(j).strGene = gpr_rules(j).strGene(gpr_rules(j).genes ~= i);
                gpr_rules(j).genes = setdiff(gpr_rules(j).genes,i);
            end
        end
    end
    % 2.4 if all genes of one rule are notknockable, delete all rules for the same reaction
    % because the reaction can never be knocked out.
    % Here, through preprocessing, some rules don't contain any genes anymore (genes = []) 
    % and are thus also notknockable
    for i = 1:length(gpr_rules)
        if all(isnan(gkoCost(gpr_rules(i).genes))) && all(isnan(gkiCost(gpr_rules(i).genes))) % finds also 'empty' (genes = []) rules 
            rule_removelist([gpr_rules(:).reaction] == gpr_rules(i).reaction) = 1;
        end
    end
    gpr_rules = gpr_rules(~rule_removelist);

    % genes that don't occur in rules anymore are notknockable
    gkoCost(setdiff(1:length(gkoCost),[gpr_rules(:).genes])) = nan;
    gkiCost(setdiff(1:length(gkiCost),[gpr_rules(:).genes])) = nan;

    if ~isempty(gpr_rules)
        % 2.6 Remove non-minimal GPR rules
        gpr_rules = remove_nonminimal(gpr_rules); 
        % 2.7 Join genes, adapt reaction's ko- and kiCosts
        [gene_subst_mat_AND,gpr_rules,gkoCost,gkiCost] = unite_genes_AND(gpr_rules,gkoCost,gkiCost);
        [gene_subst_mat_AND_2,gene_subst_mat_OR, gpr_rules,gkoCost,gkiCost] = unite_genes_OR( gpr_rules,gkoCost,gkiCost,cnap);
    else
        gene_subst_mat_AND = eye(length(genes));
        gene_subst_mat_AND_2 = eye(length(genes));
    end
else
    gene_subst_mat_AND = eye(length(genes));
    gene_subst_mat_AND_2 = eye(length(genes));
    gene_subst_mat_OR = eye(length(genes));
end

if options.pre_GPR_network_compression
%% compress network if indicated
    displ('Compressing GEM model.',verbose);
    [        cmp1_cnap, cmp1_T, cmp1_D, cmp1_rkoCost, cmp1_rkiCost, cmp1_transf_mat, cmp1_gkoCost, cmp1_gkiCost, gpr_rules ] = ... Compression
    compress(cnap, T, D, rkoCost, rkiCost,           reac_off, gkoCost, gkiCost, gpr_rules);
    text = 'compressed ';
else
    text = '';
    cmp1_transf_mat = eye(cnap.numr);
    cmp1_rkiCost = rkiCost;
    cmp1_rkoCost = rkoCost;
    cmp1_gkoCost = gkoCost;
    cmp1_gkiCost = gkiCost;
    cmp1_T = T;
    cmp1_D = D;
    cmp1_cnap = cnap;
end

%% 3. incorporate genes and enzymes
displ('Generating GPR-extended metabolic model.',verbose);
if ~isempty(gpr_rules)
    cmp1_cnap.has_gui = 0 ;
    [ cmp2_cnap, M_GEM_map, ~,koCost, kiCost, cmp2_T, cmp2_D, rType ] = ...
    CNAintegrateGPRrules( cmp1_cnap, gpr_rules, cmp1_rkoCost, cmp1_rkiCost, cmp1_T, cmp1_D, cmp1_gkoCost, cmp1_gkiCost, 0);
else
    M_GEM_map = eye(cmp1_cnap.numr);
    rType = repmat('r',1,cmp1_cnap.numr);
    cnap.rType = rType;
    kiCost = cmp1_rkiCost;
    koCost = cmp1_rkoCost;
end

% test again feasibility of target and desired vectors
displ(['Verifying D and T region feasibility in ' text 'GPR-extended model.'],verbose);
feas_D = testRegionFeas(cmp2_cnap,cmp2_D,d,2);
feas_T = testRegionFeas(cmp2_cnap,cmp2_T,t,2);
if any(~feas_T)
    error(['At least one target region (T' num2str(find(~feas_T)) ') is infeasible in the original model.']);
end
if any(~feas_D)
    error(['At least one desired region (D' num2str(find(~feas_D)) ') is infeasible in the original model.']);
end

%% 4,5,6. compute MCS
[r_mcs, status, cmp_mcs, cmp_cnap, mcs_idx_cmp, obj] = ...
    CNAMCSEnumerator2(cmp2_cnap,cmp2_T,t,cmp2_D,d,koCost,kiCost,maxSolutions,maxCost,options,verbose);

%% 7. expand to gene rules
if status ~= 2 && ~all(all(isnan(r_mcs)))
    mcs_idx_cmp_full = 1:size(r_mcs,2);
    r_mcs_exp = sparse(r_mcs);
    r_mcs_exp(isnan(r_mcs)) = 0; % for expansion non-knock-ins are first set to 0
    displ('Expand MCS.',verbose);
    % 1st undo gene extension and compression (no expansion, because knockable reactions without genes were not compressed)
    mcs_reac  = sparse(cmp1_transf_mat*M_GEM_map*r_mcs_exp);
    mcs_genes = sparse(r_mcs_exp(rType == 'g',:));
    % prepare expansion:
    % replacements
    % OR -> KI are expanded, AND -> KO are substituted
    %                             number of possible ways to employ functional knock-out/in
    % substitution array {    1    x      {  1     x     n  } }
    %                                             genes needed for one particular way
    gene_subst_mat_AND_combined = logical(double(gene_subst_mat_AND)*double(gene_subst_mat_AND_2));
    subst = cell(1,sum(rType == 'g'));
    for i = 1:sum(rType == 'g')
        if ~isnan(kiCost(sum(rType ~= 'g')+i)) % KI
            % find rule: (X OR Y OR Z), KI of only one needed for guaranteeing function
            subst{i} = num2cell(  find(gene_subst_mat_OR(:,i))  )'; % KI:OR
            for j = 1:size(subst{i},2)
                subst{i}{j} = find(gene_subst_mat_AND_combined(:,subst{i}{j}))'; % KI:AND
            end
        elseif ~isnan(koCost(sum(rType ~= 'g')+i)) % KO
            % find rule: (X OR Y OR Z), KO of all needed for blocking function
            subst(i) = { find(gene_subst_mat_OR(:,i))' }; % KO:OR
            % find all genes (Xa AND Xb AND Xc) for each sub-rule (X OR Y OR Z) and explore all possible
            % combinations to falsify the rule
            summands = cell(1,size(subst{i},2));
            for j = 1:size(subst{i},2)
                summands{j} = find(gene_subst_mat_AND_combined(:,subst{i}(j)))'; % KO:AND
            end
            combinations = combvec(summands{:})';
            combinations = arrayfun(@(x) unique(combinations(x,:)),1:size(combinations,1),'UniformOutput',false); % eliminate reduncancies
            isminimal = true(1,length(combinations));
            for j = 1:length(combinations) % delete non-minimal gene-ko combinations
                min_combin = find(isminimal);
                other_combinations_to_compare = min_combin(min_combin>j);  % compare only combinations that haven't been compared yet
                for k = other_combinations_to_compare
                    if ismember(combinations{k},combinations{j}) % find out if one set contains the other
                        isminimal(j) = false;
                    elseif ismember(combinations{j},combinations{k})
                        isminimal(k) = false;
                    end
                end
            end
            subst{i} = combinations(isminimal);
        end
    end
    % anticipate final size of MCS matrix
    num_mcs_expanded = nan(size(mcs_genes,2),1);
    for i = 1:size(mcs_genes,2)
        num_mcs_expanded(i) = prod(cellfun(@length,subst(logical(mcs_genes(:,i)))));
    end
    
    % expand mcs if expasion would lead to up to 10 000 MCS, otherwise return the compressed solution.
    if ~isempty(gpr_rules)
        if sum(num_mcs_expanded) <= 10000
            % prepare model output
            cnap_orig.has_gui = 0 ;
            % Preparing to return solutions
            displ('Generating full/uncompressed GPR-integrated model.',verbose);
            [ full_cnap, rmap, gmap, ~,~, full_T, full_D, full_rType ] = CNAintegrateGPRrules( cnap_orig, gpr_rules_orig, [], [], full_T, full_D );
            % translate kiCost and koCost
            [a,b] = find(rmap);
            % map ki- and ko-cost to correct splitting
            full_kiCost(b) = full_kiCost(a);
            full_koCost(b) = full_koCost(a);
            ivCost = [full_kiCost nan(1,sum(full_rType=='p')) full_gkiCost];
            ivCost(~isnan(full_koCost))  = full_koCost(~isnan(full_koCost));
            ivCost([false(1,sum(full_rType~='g')) ~isnan(full_gkoCost)]) = full_gkoCost(~isnan(full_gkoCost));
            % expand mcs
            mcs_genes2 = sparse([],[],[],size(gene_subst_mat_AND_combined,1),size(mcs_genes,2));
            for i = find(any(mcs_genes,2))'
                num_split = length(subst{i})*(mcs_genes(i,:)~=0) + double(mcs_genes(i,:)==0); % how many times each mcs is copied
                a = repelem_loc((1:size(mcs_genes,2)) .* double(mcs_genes(i,:)~=0),num_split); % map of reactions before and after splitting
                occ = occurcount( full(a) ) .* double(a'~=0); % counter for each mcs (if split 3 times: (R1)1, (R1)2, (R1)3)
                mcs_genes = repelem_loc(mcs_genes,1,num_split); % make copies of mcs that need to be split
                mcs_genes2 = repelem_loc(mcs_genes2,1,num_split); % make copies of mcs that need to be split
                for j = setdiff(unique(occ)',0)
                    mcs_genes2(subst{i}{j}, occ == j) = repmat(mcs_genes(i, occ == j),size(subst{i}{j},2),1);
                end
                mcs_idx_cmp_full = repelem_loc(mcs_idx_cmp_full,num_split);
            end
            mcs_reac2 = mcs_reac(:,mcs_idx_cmp_full); % copy also mcs of reaction part to fit with the gene part

            % Preparing output and eliminating MCS that are too expensive
            mcs_expanded = [rmap(:,full_rType=='r')'*mcs_reac2; zeros(sum(full_rType=='p'),size(mcs_reac2,2)) ; mcs_genes2];
            if maxCost == inf
                maxCost = 1e9; % workaround
            end
            cost_filter = ivCost;
            cost_filter(isnan(cost_filter))    = maxCost+1; % higher number, so that mcs containing this intervention are sortet out
            % check if knock-out and knock-in costs are still met
            mcs_cost = cost_filter*(mcs_expanded~=0 & ~isnan(mcs_expanded));
            mcs_affordable =  mcs_cost <= maxCost;
            for i = find(~isnan(full_kiCost)) % when KI-reactions were not knocked-in, set them to NaN in output
                mcs_expanded(i,mcs_expanded(i,:) == 0) = nan;
            end
            full_mcs = mcs_expanded(:,mcs_affordable);
            mcs_idx_cmp_full = mcs_idx_cmp(mcs_idx_cmp_full(mcs_affordable));
            gpr_rules = gpr_rules_orig;
            displ(['MCS found: ' num2str(size(full_mcs,2)) '. (' num2str(size(r_mcs_exp,2)) ' before GPR decompression / MCS expansion).'],verbose);
            % prepare other output
            full_cnap.mcs.kiCost = full_kiCost;
            full_cnap.mcs.koCost = full_koCost;
        else
            [ full_cnap, rmap, gmap, full_koCost, full_kiCost, full_T, full_D, full_rType ] = ...
                CNAintegrateGPRrules( cmp1_cnap, gpr_rules, cmp1_rkoCost, cmp1_rkiCost, cmp1_T, cmp1_D, cmp1_gkoCost, cmp1_gkiCost);
            mcs_reac  = sparse(cmp1_transf_mat*M_GEM_map*r_mcs_exp);
            mcs_genes = sparse(r_mcs_exp(rType == 'g',:));
            full_mcs = [rmap(:,full_rType=='r')'*mcs_reac; zeros(sum(full_rType=='p'),size(mcs_reac,2)) ; mcs_genes];
            ivCost = full_kiCost;
            ivCost(isnan(ivCost))  = full_koCost(isnan(ivCost));
            for i = find(~isnan(full_kiCost)) % when KI-reactions were not knocked-in, set them to NaN in output
                full_mcs(i,full_mcs(i,:) == 0) = nan;
            end
            displ(['MCS found: ' num2str(size(full_mcs,2)) '. (MCS were not decompressed because expansion would have returned ' ...
                num2str(sum(num_mcs_expanded)) ' MCS).'],verbose);
            % prepare other output
            full_cnap.mcs.kiCost = full_kiCost;
            full_cnap.mcs.koCost = full_koCost;
            full_cnap.mcs.T = cmp1_T;
            full_cnap.mcs.D = cmp1_D;
        end
        full_cnap.mcs.ivCost = ivCost;
        full_cnap.mcs.rmap = rmap;
        full_cnap.mcs.gmap = gmap;
        full_cnap.mcs.T = full_T;
        full_cnap.mcs.t = t;
        full_cnap.mcs.D = full_D;
        full_cnap.mcs.d = d;
        full_cnap.mcs.gpr_rules = gpr_rules;
        displ('Generating reaction representations for all gene MCS that correspond to the original model.',verbose)
        rmcs = gmcs2rmcs(full_mcs,gpr_rules,rmap,full_cnap.rType);
    else
        full_cnap = cnap_orig;
        full_mcs = r_mcs;
        ivCost = rkiCost;
        ivCost(isnan(ivCost))  = rkoCost(isnan(ivCost));
        full_cnap.mcs.ivCost = ivCost;
        full_cnap.mcs.rmap = eye(full_cnap.numr);
        full_cnap.mcs.gmap = [];
        full_cnap.mcs.kiCost = rkiCost;
        full_cnap.mcs.koCost = rkoCost;
        full_cnap.mcs.T = T;
        full_cnap.mcs.t = t;
        full_cnap.mcs.D = D;
        full_cnap.mcs.d = d;
        full_rType = repmat('r',1,full_cnap.numr);
        full_cnap.rType = full_rType;
    end

else
    full_cnap = CNAintegrateGPRrules( cmp1_cnap, gpr_rules, cmp1_rkoCost, cmp1_rkiCost, cmp1_T, cmp1_D, cmp1_gkoCost, cmp1_gkiCost);
    switch status
        case 1
            displ('No gene MCS found: Timeout.',verbose);
        case 2
            displ('No gene MCS found: Infeasible.',verbose);
    end
    full_mcs = double.empty(full_cnap.numr,0);
    rmcs = double.empty(cnap_orig.numr,0);
    mcs_idx_cmp_full = nan;
end

end

%% 1. compress
function [cmp_cnap, cmp_T, cmp_D, cmp_koCost, cmp_kiCost, cmp_mapReac, gkoCost, gkiCost, gpr_rules] = compress(cnap,T,D,koCost,kiCost, r_off, gkoCost, gkiCost, gpr_rules)
    %% Prepare / Protect some reactions from compression
    % identify essential reactions and adapt notknockable vector
    non_compress_reacs = any([cell2mat(D') ; cell2mat(T')],1);
    non_compress_reacs(~isnan(kiCost)) = true; % don't compress knock-in-able
    % don't compress knockable genes without gene affiliation (e.g. O2 uptake)
    non_compress_reacs(setdiff(find(~isnan(koCost)),[gpr_rules(:).reaction])) = true;
    non_compress_reacs = find(non_compress_reacs);

    javastderr= java.lang.System.err;
    java.lang.System.setErr(java.io.PrintStream('cplex_stderr.log'));
    %% Compress
    [~,~,cmp_mapReac,~,cmp_cnap] = CNAcompressMFNetwork(cnap,non_compress_reacs,[],1,0,1,r_off,0);
    java.lang.System.setErr(javastderr);
    % remap MCS-region vectors
    cmp_T = cellfun(@(x) x*cmp_mapReac,T,'UniformOutput',0);
    cmp_D = cellfun(@(x) x*cmp_mapReac,D,'UniformOutput',0);
    lumpedReacs = double(cmp_mapReac ~= 0);
    lumpedReacs(lumpedReacs == 0) = nan;
    cmp_kiCost = arrayfun(@(x) sum(kiCost(~isnan(lumpedReacs(:,x)))),1:size(lumpedReacs,2)); % for KI, all lumped reactions need to be knocked in
    cmp_koCost = min(lumpedReacs.*koCost'); % otherwise min value would be 0
    %% Compress rules - get gene-reaction rules for old model
    % rule: cell array
    %      {  1      x     n{  1    x     m  }  }
    %                  n: number of terms connected by OR
    %                            m: indices of genes connected by AND 
    numgenes = length(gkoCost);
    rule = repmat({logical.empty(numgenes,0)},1,size(cmp_mapReac,1)); % empty rules. For each reaction
    for i = unique([gpr_rules(:).reaction])
        for j = 1:length(gpr_rules)
            if i == gpr_rules(j).reaction
                rule{i} = [rule{i} full(sparse(gpr_rules(j).genes,1,true,numgenes,1))]; % add rules
            end
        end
    end
    %% Compress enzymes - map to compressed model (rules of lumped reactions are joined with AND)
    emptyRules = cellfun(@isempty,rule)';
    ruleMat = repmat({[]},1,cmp_cnap.numr);
    for i = 1:cmp_cnap.numr
        reacs = find(logical(cmp_mapReac(:,i)) & ~emptyRules); % find reaction rules to join
        if ~isempty(reacs)
            ruleMat{i} = logical(rule{reacs(1)}); % start with reaction rule of first of the lumped reactions
            for j = 2:length(reacs) % multiply rule for every AND that occurrs (equivalent to expansion of (a+b)*c = a*c + b*c)
                ruleMat{i} = repelem_loc(ruleMat{i},1,size(rule{reacs(j)},2)) | repmat(rule{reacs(j)},1,size(ruleMat{i},2));
                if size(rule{reacs(j)},2) > 1 || j == length(reacs) % reduce number of rules if possible, everytime their number increases through expansion
                    ruleMat{i} = unique(ruleMat{i}','rows')';
                    ismin = false(1,size(ruleMat{i},2));
                    for k = 1:size(ruleMat{i},2) % mark and eliminate redundand and non-minimal gene rules:
                        ismin(k) = sum(all(ruleMat{i} >= ruleMat{i}(:,k)))==1; % (A and B) or B -> B
                    end                                                        % thus rule (A and B) can be eliminated
                    ruleMat{i} = ruleMat{i}(:,ismin); % only the 'minimal' rules are kept (B)
                end
            end
        end
    end
    [a,b] = unique([gpr_rules(:).genes]);
    gname = [gpr_rules(:).strGene]';
    gNames = repmat({''},1,length(gkoCost));
    gNames(a) = gname(b);
    gpr_rules = gpr_rules([]);
    c = 1;
    for i = find(~cellfun(@isempty,ruleMat))
        for j = 1:size(ruleMat{i},2)
            gpr_rules(c).reaction = i;
            gpr_rules(c).strReac = cellstr(cmp_cnap.reacID(i,:));
            gpr_rules(c).genes = find(ruleMat{i}(:,j))';
            gpr_rules(c).strGene = gNames(gpr_rules(c).genes);
            gpr_rules(c).name = {['Rule-' strjoin(gpr_rules(c).strGene,'-') '-r' strtrim(cmp_cnap.reacID(i,:))]};
            c = c+1;
        end
    end
    % sort rules for genes (actually not necessary, but more concurrent with original order)
    [~,b] = sort(cellfun(@min,{gpr_rules(:).genes}));
    gpr_rules = gpr_rules(b);
end
%% 2. generate split bounds
function [splt_lb, splt_ub] = remapBounds(lb,ub,map)
    [splt_lb, splt_ub] = deal(nan(size(map,2),1));
    for k=1:size(map,2)
        forw=find(map(:,k)>0);
        revs=find(map(:,k)<0);
        if(~isempty(forw))
            splt_lb(k) = max(lb(forw)./map(forw,k));
            splt_ub(k) = min(ub(forw)./map(forw,k));
        end
        if(~isempty(revs))
            splt_lb(k) = max([splt_lb(k);ub(revs)./map(revs,k)]);
            splt_ub(k) = min([splt_ub(k);lb(revs)./map(revs,k)]);
        end
    end
end
%% 3.1 remove non-minimal rules
% A or (A and B) = A
function gpr_rules = remove_nonminimal(gpr_rules)
    maxreac = max([gpr_rules(:).reaction]);
    reac_rule = sparse([gpr_rules(:).reaction],1:length(gpr_rules),1,maxreac,length(gpr_rules));
    [~,~,rule_reac_abund] = unique(reac_rule','rows','stable'); 
    [rule_abund,isorule_map] = hist(rule_reac_abund,unique(rule_reac_abund));
    isorule_grule_ismin = false(1,length(gpr_rules));
    % 1. for non-isoenzymes, the only existing gene rule is the shortest
    isorule_grule_ismin( ismember(rule_reac_abund', isorule_map(rule_abund == 1)') ) = 1;
    % 2. for isoenzymes, gene rules are compared to eliminate redundancies
    for i = isorule_map(rule_abund > 1)'
        genes_i = {gpr_rules(rule_reac_abund == i).genes};
        isoenz = find(rule_reac_abund == i);
        isoenz_gene_rule = sparse(  cell2mat(genes_i),...
                                    repelem_loc(1:length(isoenz), cellfun(@length,genes_i)),...
                                    1);
        % remove redundancies
        for j = 1:length(isoenz) % mark and eliminate redundant and non-minimal gene rules:
            if sum(all(repmat(isoenz_gene_rule(:,j),1,size(isoenz_gene_rule,2)) >= isoenz_gene_rule))==1  % A or (A and B) = A
                isorule_grule_ismin(isoenz(j)) = 1; % rule is minimal e.g. A
            else
                isorule_grule_ismin(isoenz(j)) = 0; % rule is not minimal e.g. (A and B)
            end
        end
    end
    gpr_rules = gpr_rules(isorule_grule_ismin); % remove redundant enzymes -> this also reduces set of genes
end
%% 3.2 unite and substitute genes if they occur always and only together with each other
% genes, connected with AND -> koCost are kept the same (lowest), kiCost are summed up
function [gene_subst_mat,gpr_rules,gkoCost,gkiCost] = unite_genes_AND(gpr_rules,gkoCost,gkiCost)
    genes_old = 1:length(gkoCost);
    [a,b] = unique([gpr_rules(:).genes]); % ,'stable'
    gname = [gpr_rules(:).strGene]';
    gNames_old = repmat({''},1,length(genes_old));
    gNames_old(a) = gname(b);
    % remove redundant AND associations:
    gene_rule = cell2mat(cellfun(@(x) sparse(x,1,1,length(genes_old),1),   {gpr_rules(:).genes}   ,'UniformOutput',0));
    gene_rule = logical(gene_rule);
    [a,~,gene_subst_mat] = unique(gene_rule,'rows','stable'); % unite equivalent genes
    gene_subst_mat = full(sparse(gene_subst_mat, 1:length(gene_subst_mat), true)); % new gene mapping    
    gene_subst_mat = gene_subst_mat(any(a,2),:); % keep genes with rule association
    genes_new = 1:size(gene_subst_mat,1); % genes_1: Genes-Vector after removing AND-redundancies
    gkoCost = arrayfun(@(x) min(gkoCost(gene_subst_mat(x,:))), genes_new); % take lowest cost for knock-out
    gkiCost = arrayfun(@(x) sum(gkiCost(gene_subst_mat(x,:))), genes_new); % take sum of costs for knock-in
    gNames_new =  arrayfun(@(x) ['(',strjoin(gNames_old(gene_subst_mat(x,:)),'_and_'),')'],genes_new,'UniformOutput',0);
    % reassign gene indices
    for i = 1:length(gpr_rules)
        gpr_rules(i).genes = unique(arrayfun(@(x) find(gene_subst_mat(:,x)), gpr_rules(i).genes),'stable');
        gpr_rules(i).strGene = gNames_new(gpr_rules(i).genes);
    end
    gene_subst_mat = gene_subst_mat';
end
%% 3.3 unite and substitute genes if they occur only in one reaction
% genes, connected with OR -> koCost are summed up, kiCost are kept the same (lowest)
function [genes_to_join_AND,genes_to_join_OR,gpr_rules,gkoCost,gkiCost] = unite_genes_OR(gpr_rules,gkoCost,gkiCost,cnap)
    numGenes_old = length(gkoCost);
    numGenes     = numGenes_old;
    genes_old = 1:numGenes_old;
    gKoCost_old_no_nan = gkoCost;
    gKoCost_old_no_nan(isnan(gkoCost)) = 0;
    [a,b] = unique([gpr_rules(:).genes]);
    gname = [gpr_rules(:).strGene]';
    gNames = repmat({''},1,length(genes_old));
    gNames(a) = gname(b);
    maxReac = max([gpr_rules(:).reaction]);
    ruleMat = repmat({[]},1,maxReac);
    for i = unique([gpr_rules(:).reaction])
        rules = [gpr_rules(:).reaction] == i;
        rule_y = {gpr_rules(rules).genes};
        rule_x = repelem_loc(1:sum(rules),cellfun(@length,rule_y));
        ruleMat{i} = sparse([rule_y{:}],rule_x,true,numGenes_old,sum(rules));
    end
    ruleMat_all = [ruleMat{:}];
    ruleMat_new = ruleMat;
    gene_abund_all = sum(ruleMat_all,2);
    % extend genes and rules vector first, then remove obsolete entries
    genes_to_join_AND = logical(eye(numGenes_old));
    genes_to_join_OR  = logical(eye(numGenes_old));
    for i = unique([gpr_rules(:).reaction]) % for all reactions with rules
        if size(ruleMat{i},2) > 1  % if reaction has multiple rules ('or')
            gene_abund_1 = sum(ruleMat{i},2);
            % check if all gene occurrencies are with the rules of this reaction. If yes, genes can
            % be joined.
            if gene_abund_1(logical(gene_abund_1)) == gene_abund_all(logical(gene_abund_1))
                for k = find(sum(ruleMat{i},1)>1) % first join genes that are connected with AND
                    % create new gene
                    numGenes = numGenes+1;
                    genes_to_join_AND( ruleMat{i}(:,k) ,numGenes) = 1;
                    gNames{numGenes} = ['[',strjoin(gNames(ruleMat{i}(:,k)),'_and_'),']'];
                    gkoCost(numGenes) = min(gkoCost(ruleMat{i}(:,k)));
                    gkiCost(numGenes) = sum(gkoCost(ruleMat{i}(:,k)));
                    % substitute genes in gene rule
                    ruleMat_new{i}(  :  ,k) = 0;
                    ruleMat_new{i}(numGenes,k) = 1;
                end
                % join genes that are connected with 'or'
                numGenes = numGenes+1;
                genes_to_join_OR( any(ruleMat_new{i},2), numGenes ) = 1;
                gNames{numGenes} = ['[',strjoin(gNames(any(ruleMat_new{i},2)),'_or_'),']'];
                ruleMat_new{i} = false(numGenes,1);
                ruleMat_new{i}(numGenes) = true;
                % calculate KO- or KI-costs
                originalRules  = genes_to_join_AND(:,genes_to_join_OR(:,ruleMat_new{i}));
                if all(~isnan(gkoCost(logical(gene_abund_1)))) % KO -> minimal hitting set problem with exhaustive search
                    occurringGenes = any(originalRules,2);
                    gkiCost(numGenes) = nan;
                    if ~isempty(which('intlinprog'))
                        % solve ILP with linintprog to find minimal costs to knock out all rules connected with 'or'.
                        [~,gkoCost(numGenes)] =  intlinprog(gKoCost_old_no_nan,... f
                                                            1:size(occurringGenes,1), ... int
                                                            -originalRules', ... A
                                                            -ones(1,size(originalRules,2)),[],[],... b
                                                            zeros(size(occurringGenes,1),1),... lb
                                                            ones(size(occurringGenes,1),1),[],... ub
                                                            optimoptions('intlinprog','Display','off'));
                    else % if intlinprog is not available find out lowest cost by generating all knockout combinations
                        geneCombinations = false(2^sum(occurringGenes),numGenes);
                        geneCombinations(1:2^sum(occurringGenes),occurringGenes) = dec2bin(2^sum(occurringGenes)-1:-1:0)-'0';
                        cost = nan(size(geneCombinations,1),1);
                        for j = 1:size(geneCombinations,1)
                            if all(any(originalRules(geneCombinations(j,:),:),1)) % if all rules are hit
                                cost(j) = sum(gkoCost(geneCombinations(j,:)));
                            end
                        end
                        gkoCost(numGenes) = min(cost);
                    end
                elseif all(~isnan(gkiCost(logical(gene_abund_1)))) % KI -> minimal hitting set problem
                    gkoCost(numGenes) = nan;
                    cost = nan(size(originalRules,2),1);
                    for j = 1:size(originalRules,2)
                        cost(j) = sum(gkiCost(originalRules(:,j)));
                    end
                    gkiCost(numGenes) = min(cost);
                else
                    warning(['Ambigous knockability indication. Reaction ' num2str(i) ...
                             ' has genes that can be knocked in and others that can be knocked out.']);
                end
            end
        end
    end
    gpr_rules = gpr_rules([]);
    c = 1;
    for i = find(~cellfun(@isempty,ruleMat_new))
        for j = 1:size(ruleMat_new{i},2)
            gpr_rules(c).reaction = i;
            gpr_rules(c).strReac = cellstr(cnap.reacID(i,:));
            gpr_rules(c).genes = find(ruleMat_new{i}(:,j))';
            gpr_rules(c).strGene = gNames(gpr_rules(c).genes);
            gpr_rules(c).name = {['Rule-' strjoin(gpr_rules(c).strGene,'-') '-r' strtrim(cnap.reacID(i,:))]};
            c = c+1;
        end
    end
    % sort rules for genes
    [~,b] = sort(cellfun(@min,{gpr_rules(:).genes}));
    gpr_rules = gpr_rules(b);
    gkoCost(setdiff(1:length(gkoCost),[gpr_rules(:).genes])) = nan; % genes that don't occur in enzymes anymore
    gkiCost(setdiff(1:length(gkiCost),[gpr_rules(:).genes])) = nan; % are set to notknockable (and are not added to the model later).
end

%% 4. Convert gene-MCS to reaction-MCS that match the original model
function rmcs = gmcs2rmcs(gmcs,gpr_rules,rmap,rType)
    ext2old = arrayfun(@(x) find(rmap(x,:),1),1:size(rmap,1));
    rmcs2 = gmcs(ext2old,:);

    %% translate geneMCS
    numr = size(rmap,1);
    gmcs = gmcs(rType == 'g',:);

    % identified knocked-out or knocked-in rules
    rmcs1 = zeros(numr,size(gmcs,2));
    for i = 1:size(gmcs,2)
        rule_ki_ko = zeros(1,length(gpr_rules));
        for j = 1:length(gpr_rules)
            if any(gmcs(gpr_rules(j).genes,i) == -1) && ~any(gmcs(gpr_rules(j).genes,i) == 1) && all(~isnan(gmcs(gpr_rules(j).genes,i)))
                rule_ki_ko(j) = -1;
            elseif all(gmcs(gpr_rules(j).genes,i) == 1) && ~any(isnan(gmcs(gpr_rules(j).genes,i)))
                rule_ki_ko(j) = 1;
            elseif any(isnan(gmcs(gpr_rules(j).genes,i)))
                rule_ki_ko(j) = nan;
            end
        end
        % trace back knocked-out or knocked-in reactions
        rlist = [gpr_rules(:).reaction];
        for j = unique(rlist)
            ko_ki_state = rule_ki_ko(j == rlist);
            if all(ko_ki_state == -1)
                rmcs1(j,i) = -1;
            elseif any(ko_ki_state == 1)
                rmcs1(j,i) = 1;
            elseif any(isnan(ko_ki_state)) && ~any(ko_ki_state == 1)
                rmcs1(j,i) = nan;
            elseif all(ko_ki_state == 0 | ko_ki_state == -1)
                rmcs1(j,i) = 0;
            else
                warning('this statement should not be reached');
            end
        end
    end
    rmcs = rmcs1+rmcs2;
end

%% === utilities ===

%% occurrency counter [1 2 2 3 3] -> [1 1 2 1 2]
function occ_count = occurcount(y)
    occ_count = sum(cell2mat(arrayfun(@(x) cumsum(abs(y)==x).*(abs(y)==x),unique(y),'UniformOutput',0)'),1)';
end

%% local implementation of the MATLAB repelem function.
function U = repelem_loc(V,varargin)
% This function is used to guarantee compatibility with MATLAB2014b and below
    if exist('repelem','builtin')
        U = repelem(V,varargin{:});
        return;
    end
    % if V = 1xn or nx1 and N = 1xn
    if size(V,1) == 1 && length(varargin) == 1
        varargin{2} = 1;
        varargin = flip(varargin);
    elseif size(V,2) == 1 && length(varargin) == 1
        varargin{2} = 1;
    end
    % if V = 1xn and N = 1x1
    for i = find(cellfun(@length,varargin) == 1)
        varargin{i} = varargin{i}*ones(1,getelements(size(V)',i));
    end
	[reps{1:numel(varargin)}] = ndgrid(varargin{:});
	U = cell(cellfun(@length,varargin));
	U(:) = arrayfun(@(x) V(x)*ones(cellfun(@(y) y(x),reps)) ,1:numel(V),'UniformOutput',false);
	U = cell2mat(U);
end

%% 7. local implementation of the MATLAB combvec function.
function y = combvec(varargin)
    % Mark Beale, 12-15-93
    % Copyright 1992-2010 The MathWorks, Inc.
    if isempty(varargin)
        y = [];
    else
        y = varargin{1};
        for i=2:length(varargin)
            z = varargin{i};
            y = [copy_blocked(y,size(z,2)); copy_interleaved(z,size(y,2))];
        end
    end
end

function b = copy_blocked(m,n)
    [mr,mc] = size(m);
    b = zeros(mr,mc*n);
    ind = 1:mc;
    for i = (0:(n-1))*mc
        b(:,ind+i) = m;
    end
end
function b = copy_interleaved(m,n)
    [mr,mc] = size(m);
    b = zeros(mr*n,mc);
    ind = 1:mr;
    for i = (0:(n-1))*mr
        b(ind+i,:) = m;
    end
    b = reshape(b,mr,n*mc);
end