function [cnap, enzymes, genes, gpr_rules] = CNAgenerateGPRrules( cnap, grRules, bool_enz_only )
%
% ------------------------------------------------
% CellNetAnalyzer API function 'CNAgenerateGPRrules'
% ------------------------------------------------
% --> Creates and enzyme or GRP-rules struct-array.
%
% Usage: [cnap, enzymes, genes, gpr_rules] = CNAgenerateGPRrules( cnap, grRules, bool_enz_only )
% 
% generates a struct array that contains the information on
% gene-protein-reaction-associations. Gene-information can either be
% provided in the cnap-structure as the value of the field 'geneProductAssociation' 
% (see CNAsetGenericReactionData), or provided seperately as a cell array.
% Syntax: gene names don't contain space characters, 
%         'and', 'or' and curved brackets are reserved delimiters.
%    e.g. (g01 and g02) or g03
%
% This fuction analyzes the logical structure of the
% gene-reaction-associations to predict "enzymes" or "pseudo-enzymes".
% To detect enzymes correctly it is necessary that all gene-rules are
% fully expanded to DNF. That means no 'or' operator within brackets:
%  if ORIGINAL:      ((g01 and g02) or (g01 and g03)) and g04
% YES (expanded):    (g01 and g02 and g04) or (g01 and g03 and g04)
% NO  (shortest):     g01 and g04 and (g02 or g03)
%
% See 'Output' for details in the format.
%
% ------------------------------------------------
% Input:
%   cnap           : CellNetAnalyzer mass flow project
%   grRules        : <cell>[numReac x 1] If GPR-rules are not contained in the generic reaction data 
%                     (see CNAgetGenericReactionData) in the field 'geneProductAssociation', 
%                     they can be provided in a cell array in text form for each reaction
%                     (see above)
%   bool_enz_only  : If this flag is set, only enzyme pools are returned instead
%                	  of genes (with enzyme pools only as intermediate
%                     step). Not used for Gene-MCS, but useful if only
%                     KO-feasibility should be assessed.
%   
% ------------------------------------------------
% Output:
%   cnap           : CellNetAnalyzer mass flow project with 'genes', 'enzymes' and
%                     'gpr_rules' as additional struct fields.
%   enzymes: struct array. Every entry represents a different enzyme.
%            The following fields define the enzyme:
%               - name: Automatically generated Name for enzyme (from genes and reactions)
%               - reactions: double array that contains the reaction
%                            indices of the catalyzed reactions (referring
%                            to the CNA project)
%               - genes: double array that contains the indices of the
%                        participating genes (referring to the genes cell-array)
%               - strReac: cell array that contains reaction names (for verification)
%               - strGene: cell array that contains gene names (for verification)
%            example 1: r1 (g1 & g2), r2 (g1 & g2) translates to one enzyme:
%                       Enz1 - reactions: [1,2], genes: [1,2]
%            example 2: r1 (g1 | g2), r2 (g2) translates to 2 enzymes:
%                       Enz1 - reacs: [1], genes: [1]
%                       Enz2 - reacs: [1,2], genes: [2]
%
%   genes             : Cell array with gene names. Order matches with the indices in 
%                        "enzymes" and "gpr_rule" structs.
%   gpr_rules         : Almost identical to enzymes struct array. But every entry
%                        represents a sub-GPR-rule of ONLY ONE reaction. Fields:
%                         - name
%                         - reaction: only contains one reaction index
%                         - genes
%                         - strReac: only contains one reaction name
%                         - strGene
%            example 1: r1 (g1 & g2), r2 (g1 & g2) translates to 2 GPR rules:
%                       GPR1 - reaction: [1], genes: [1,2]
%                       GPR2 - reaction: [2], genes: [1,2]
%            example 2: r1 (g1 | g2), r2 (g2) translates to 3 GPR rules:
%                       GPR1 - reac: [1], genes: [1]
%                       GPR2 - reac: [1], genes: [2]
%                       GPR3 - reac: [2], genes: [2]
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
%
%% if not provided separately as cellstr array, read out gene rules from reacNotes
if nargin == 1
    grRules = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
    bool_enz_only = 0;
end

if nargin >= 2
    if isempty(grRules)
        grRules = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
    end
end

% enumerate all occurring genes, generate a vector where all of them occur
% only once.
genesStr = {};
for i = 1:length(grRules)
    genesStr = [genesStr; split_loc(regexprep(grRules(i),'\(|\)',''))]; % regexprep removed brackets
end
genesStr = unique(genesStr);
genesStr = genesStr(~strcmp(genesStr,'and')&~strcmp(genesStr,'or')&~strcmp(genesStr,''));

% find all enzymes. Gene rules are evaluated. Genes that are connected
% through an 'and' are concatenated (with '-' as delimiter) and regarded as one enzyme.
enzymeReacAssoc = {};
enzymeGeneAssoc = {};
enzRules = regexprep(regexprep(regexprep(grRules,'  *and  *','-'),'  *or  *',' '),'\(|\)','');
for i = 1:length(grRules)
     eRule = {split_loc(enzRules(i))}; % get separate enzymes that catalyze this reaction
     for j=1:length(eRule{:})
        enzymeReacAssoc(i,j) = eRule{1}(j); % enter enzymes into a Matrix 
                                            % (reacs -> list of catalyzing enzymes)
        enzymeGeneAssoc(i,j) = {findStrPos(genesStr,split_loc(eRule{1}{j},'-'))}; % enter
                                              % corresponding gene-enzyme-rules at same
                                              % position in the matrix
     end
end
% generate a vector containing all occurring enzymes only once.
enzymesStr = reshape(enzymeReacAssoc,1,size(enzymeReacAssoc,1)*size(enzymeReacAssoc,2));
enzymesStr = unique(enzymesStr(~cellfun(@isempty,enzymesStr)));
enzymesStr = enzymesStr(~strcmp(enzymesStr,''));
emptyVec = cell(length(enzymesStr),1);
enzymes = struct('name',strcat('Enz-', enzymesStr(:)),...
                 'reactions',emptyVec,'genes',emptyVec,'strReac',emptyVec,'strGene',emptyVec);
             
if isempty(gcp('nocreate'))
    parforArg = 0;
else
    parforArg = inf;
end
rIDs = cnap.reacID;
parfor (i=1:length(enzymesStr) , parforArg)
    % find enzyme-reaction association
    [reac,enzInReac]   = find(strcmp(enzymeReacAssoc,enzymesStr(i)));
    enzymes(i).reactions = reac';
    % find corresponding gene-enzyme association (on 
    genes = enzymeGeneAssoc{reac(1),enzInReac(1)}
    enzymes(i).genes     = genes;
    % additionally also save gene-enzyme-reaction associations as text
    enzymes(i).strReac   = cellstr(rIDs(reac',:))';
    enzymes(i).strGene   = cellstr(genesStr(genes))';
end

if bool_enz_only
    genesStr = {};
    for i = 1:length(enzymes)
        name = char(enzymes(i).name);
        genesStr(i,:) = [name {enzymes(i).strGene}];
        enzymes(i).strGene = cellstr(name(5:end));
        enzymes(i).genes = i;
    end
end

genes = genesStr;
if ~isempty(enzymes)
    rulecount = repelem_loc(1:length(enzymes),cellfun(@length,{enzymes(:).reactions}));
    countOcc = occurcount(rulecount);
    for i = 1:length(rulecount)
        gpr_rules(i).reaction  = enzymes(rulecount(i)).reactions(countOcc(i));
        gpr_rules(i).strReac   = enzymes(rulecount(i)).strReac(  countOcc(i));
        gpr_rules(i).genes     = enzymes(rulecount(i)).genes;
        gpr_rules(i).strGene   = enzymes(rulecount(i)).strGene;
        gpr_rules(i).name    = {['Rule-' enzymesStr{rulecount(i)} '-' num2str(gpr_rules(i).strReac{:})]};
    end
else
    gpr_rules = struct('name',{},'reaction',{},'genes',{},'strReac',{},'strGene',{});
end

cnap.enzymes  = enzymes;
cnap.genes    = genesStr;
cnap.gpr_rules = gpr_rules;
end

function indices = findStrPos( str , pattern )
% str       the space that is searched
% pattern   the search keyword or pattern
% convert str type to cell
switch class(str)
    case 'string'
        str = cellstr(strtrim(str));
    case 'char'
        str = cellstr(str);
    case 'cell'
        
    otherwise
        error('input 1 of findStrPos doesn''t correct type');
end
% convert pattern type to cell
switch class(pattern)
    case 'string'
        pattern = char(pattern);
    case 'char'
        pattern = cellstr(pattern);
    case 'cell'

    otherwise
        error('input 2 of findStrPos doesn''t correct type');
end
[rows,cols] = size(pattern);
for k = 1:(rows*cols)
    ind                      = find(strcmp(str, pattern(k)));
    indices(1:length(ind),k) = ind;
end
if (any(size(indices)==0))
    indices = [];
end
end
%% concurrency counter [1 2 2 3 3] -> [1 1 2 1 2]
function occ_count = occurcount(y)
    occ_count = sum(cell2mat(arrayfun(@(x) cumsum(abs(y)==x).*(abs(y)==x),unique(y),'UniformOutput',0)'),1)';
end

% local implementation of the MATLAB repelem function.
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
        sizeV = size(V);
        varargin{i} = varargin{i}*ones(1,sizeV(i));
    end
	[reps{1:numel(varargin)}] = ndgrid(varargin{:});
	U = cell(cellfun(@length,varargin));
	U(:) = arrayfun(@(x) V(x)*ones(cellfun(@(y) y(x),reps)) ,1:numel(V),'UniformOutput',false);
	U = cell2mat(U);
end

function [ parts ] = split_loc( str,delimiter )
% splits char or cell into multiple parts
% returns them as cell array
% str: original string
% delimiter: character (array) where to split (default is whitespace ' ')
% prepare parameters
if nargin == 1
    delimiter = ' ';
end
if ischar(delimiter)
elseif iscellstr(delimiter)
    delimiter = char(delimiter);
else
    error('split: invalid type of delimiter. Enter char or cellstr');
end
if ischar(str)
elseif iscellstr(str)
    str = char(str);
else
    error('split: invalid type of input string. Enter char or cellstr');
end
if isempty(str)
    parts = {''};
    return;
end
% replace multiple delimiters placed in line.
while ~isempty(strfind(str, [delimiter delimiter]))
    str = strrep(str,[delimiter delimiter],delimiter);
end
% start splitting
pos = strfind(str,delimiter)';
pos      = [1;pos+length(delimiter)];
pos(:,2) = [pos(2:end)-length(delimiter)-1;length(str)];
for i=1:size(pos,1)
    parts(i) = {str(pos(i,1):pos(i,2))};
end
parts = parts';
end