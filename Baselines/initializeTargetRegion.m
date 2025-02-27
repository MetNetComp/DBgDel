% %(metabolism must no longer be able to grow while not producing any product)
% % T           : <cell>[numT x 1] Matrices specifying (with vector t) the target flux
% %                regions as given above. T = {T1,T2, ...}. T1: <double>[rowsT1 x cnap.numr]
% % t           : <cell>[numT x 1] Vectors specifying (with matrix T) the target flux
% %                regions as given above. t = {t1,t2, ...}. t1: <double>[rowsT1 x 1]
% % Target Region (T*r<=t): growth >= GR_threshold, production <= PR_threshold; 
% %T*r<=t format (2xn matrix and 2x1 vector): -growth <= -GR_threshold, production <= PR_threshold

function [T, t] = initializeTargetRegion(cnap_numr, growth_idx, production_idx, GR_threshold, PR_threshold)
    % Initialize T and t as cells
    T = cell(1,1);
    t = cell(1,1);

    % Create a 2 × numr matrix for T (two constraints)
    T{1} = zeros(2, cnap_numr);
    T{1}(1, growth_idx) = -1; % -growth
    T{1}(2, production_idx) = 1; % production

    % Create a 2 × 1 vector for t
    t{1} = [-GR_threshold; PR_threshold];
end