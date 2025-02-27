function setJavaHeapSize(sizeGB)
    % setJavaHeapSize(sizeGB) - Set Java heap memory for EFMTool before computation.
    %
    % Inputs:
    %   sizeGB - Desired Java heap size in GB (e.g., 16 for 16GB)
    %
    % Example:
    %   setJavaHeapSize(32)  % Sets Java heap memory to 32GB
    % 
    % Note:
    %   Please use the Recommended Heap Size Settings for EFMTool (An error such as 'GC overhead limit exceeded' may cause the operation 
    %   to pauses or terminate and return a false status value.)

    % Total_RAM	  Recommended_HeapSize
    %   16GB RAM	-Xmx4G to -Xmx6G
    %   32GB RAM	-Xmx8G to -Xmx12G
    %   64GB RAM	-Xmx16G to -Xmx24G
    %   128GB RAM	-Xmx32G to -Xmx48G

    % Convert size to string format for Java options
    heap_size = sprintf('-Xmx%dG', sizeGB);  
    
    % Get MATLAB root directory
    matlab_root = matlabroot;
    
    % Define the correct path to java.opts on Linux
    java_opts_file = fullfile(matlab_root, 'bin', 'glnxa64', 'java.opts');

    % Check if running on Linux (only execute on Linux systems)
    if isunix
        % Construct the shell command to set Java heap size
        command = sprintf('echo "%s" | sudo tee %s', heap_size, java_opts_file);

        % Execute the command
        [status, cmdout] = system(command);

        if status == 0
            disp(['Java heap size set to ', heap_size]);
            disp('Restart MATLAB for changes to take effect.');
        else
            error('Failed to set Java heap size. Error: %s', cmdout);
        end
    else
        warning('This function is intended for Linux systems. Use MATLAB Preferences on Windows/Mac.');
    end

    % Display verification instructions
    disp('Run the following command after restarting MATLAB to check the new heap size:');
    disp('java.lang.Runtime.getRuntime.maxMemory / 1024^3');
end