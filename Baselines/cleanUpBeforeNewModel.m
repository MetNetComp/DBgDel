% cleanUpBeforeNewModel() - Cleans up memory and releases the parallel pool
%
% This function clears variables, closes figures, releases Java and MATLAB memory,
% and shuts down any active parallel pool before running a new model.
%
% Example:
%   cleanUpBeforeNewModel();  % Cleans up memory and shuts down parallel pool
% 
% Note:
%   It is recommended to run this cleanup before starting a new model to ensure that
%   the memory is released and any parallel pool is properly closed. This will help 
%   avoid memory issues and ensure smooth execution when transitioning between models.
% 
%   Additionally, releasing the parallel pool ensures that resources are freed, 
%   which may improve performance and prevent unnecessary resource consumption in subsequent runs.
% 
%   Recommended usage:
%   - Call this function before starting a new computational model to reset the environment.

function cleanUpBeforeNewModel()
    % Clear workspace, close figures, and release memory
    clear;      % Clear variables
    close all;  % Close all figure windows
    clc;        % Clear command window

    % Release Java memory
    java.lang.System.gc();   % Request Java garbage collection

    % Release MATLAB memory
    pack;   % Forces MATLAB to free memory from large variables

    % Optionally clear global variables
    clear global;

    % Release parallel pool
    pool = gcp('nocreate');  % Get current parallel pool (if any)
    if ~isempty(pool)
        delete(pool);  % Delete the pool if it exists
    end

    disp('Memory released, parallel pool closed.');
end