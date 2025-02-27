% D           : <cell>[numD x 1] Matrices specifying (with vector d) the desired flux
%                regions as given above. D = {D1,D2, ...}. D1 <double>[rowsD1 x cnap.numr]
%(metabolism must be able to grow)
% d           : <cell>[numD x 1] Vectors specifying (with matrix D) the desired flux
%                regions as given above. d = {d1,d2, ...}. d1 <double>[rowsD1 x 1]
% Desired Region (D*r<=d): growth >= GR_threshold. 
% D*r<=d (1xn matrix and 1x1 vector): -growth <= -GR_threshold

function [D, d] = initializeDesiredRegion(cnap_numr, growth_idx, GR_threshold)
    % Initialize D and d as cells
    D = cell(1,1);
    d = cell(1,1);

    % Create a 1 × numr matrix for D (one constraint)
    D{1} = zeros(1, cnap_numr);
    D{1}(1, growth_idx) = -1; % -growth

    % Create a 1 × 1 vector for d
    d{1} = -GR_threshold;
end