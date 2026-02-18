function rel_error = cellular_model_run_simulations(x)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set working directory
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    base_path = fileparts(mfilename('fullpath'));
    template_folder = fullfile(base_path, 'Original_Folder');
    iterIdx = char(java.util.UUID.randomUUID());
    work_path = fullfile(base_path, sprintf('Attempt_%s',iterIdx));
    % Copy original folder
    copyfile(template_folder, work_path);

    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Write the candidate gamma values to a new CSV file (overwrite).
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        csvPath = fullfile(work_path, 'Modules', 'tracheids_gamma_attempts.csv');
        writeGammaCSVFromInput(csvPath, x(1:4));
        fprintf('CSV file written at: %s\n', csvPath);
        pause(5);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate the Fortran include file from input values
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        incFilePath = fullfile(work_path, 'Modules', 'tracheids_gamma_values.inc');
        generateIncludeFileFromValues(x(1:4), incFilePath);
        fprintf('Include file generated at: %s\n', incFilePath);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run simulations by calling the Python script.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        input_python = 'run_for_MatSuMoTo.py';    
        system_call = sprintf('cd "%s" && python3 %s %s', work_path, input_python, work_path);
        [status, cmdout] = system(system_call, '-echo');    
        if status ~= 0
            error('Error running Python script: %s', cmdout);
        else
            disp('Python script ran successfully.');
        end
        pause(5);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read simulation results 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gr_folder = fullfile(work_path, 'GR');
        sim_res_path = fullfile(gr_folder, 'GR_gamma.csv');
        if ~exist(sim_res_path,'file')
            error('Missing results file: %s', sim_res_path);
        end
        sim_tau = [0.1, 1, 10, 100];
        t1 = 0 : 0.1 : 1;

        % second segment: 1 to 150 in steps of 1
        t2 = 1 : 1 : 150;

        % concatenate, and (optionally) remove the duplicate “1” at the join
        sim_time = [t1, t2(2:end)];
        sim_results = get_creep_curves(sim_res_path, sim_time, sim_tau);
        % ----------------------------------------------------------------
        % write simulation results out
        sim_out_path = fullfile(gr_folder, 'sim_creep_curves.csv');
        %save(sim_out_path, 'sim_results');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read target results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        target_res_path = fullfile(base_path, 'macro_master_gamma.csv');
        macro_tau = [1, 10, 100, 1000];
        target_results = get_target_curves(target_res_path, sim_time, macro_tau);
        % ----------------------------------------------------------------
        % write target results out
        target_out_path  = fullfile(base_path, 'macro_creep_curves.csv');
        %save(target_out_path, 'target_results');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate error
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        weights = [1, 0, 0, 0, 0, 0, 0, 0, 0];
        [~, rel_error, ~, ~] = compareCreepCurves(target_results, sim_results, weights);
        
        fprintf('\nRelative Euclidean error norm: %f\n\n', rel_error);
    
    catch ME
       fprintf('\nSomething went wrong...');
       rethrow(ME);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeGammaCSVFromInput(csvPath, x)
    % writeGammaCSVFromInput Writes a CSV file containing a single row of gamma values.
    %   This function overwrites any existing CSV file at csvPath.
    %
    %   Input:
    %       csvPath - Full path to the CSV file.
    %       x       - 1x4 vector of gamma values.
    %
    %   The file will have a header row and one row of data.
    if numel(x) ~= 4
        error('Input x must be a vector of 4 elements.');
    end

    fid = fopen(csvPath, 'a');  % use 'w' mode to overwrite file each time.
    if fid == -1
        error('Could not open CSV file "%s" for writing.', csvPath);
    end

    % Write header row.
    %fprintf(fid, 'gamma1\tgamma2\tgamma3\tgamma4\n');
    % Write the gamma values row.
    fprintf(fid, '%.5f\t%.5f\t%.5f\t%.5f\n', x(1), x(2), x(3), x(4));
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function generateIncludeFileFromValues(values, includePath)
    % generateIncludeFileFromValues Generates a Fortran include file using gamma values.
    %   The output file is written in a format that umat expects.
    %
    %   Input:
    %       values      - Structure with fields: gamma1, gamma2, gamma3, gamma4.
    %       includePath - Full path (including filename) for the include file.
    fid = fopen(includePath, 'w');
    if fid == -1
        error('Cannot open file "%s" for writing.', includePath);
    end

    fprintf(fid, '      ! Gamma values generated from CSV file\n');
    fprintf(fid, '      KW_Gamma_i(1) = %.5fD0\n', values(1));
    fprintf(fid, '      KW_Gamma_i(2) = %.5fD0\n', values(2));
    fprintf(fid, '      KW_Gamma_i(3) = %.5fD0\n', values(3));
    fprintf(fid, '      KW_Gamma_i(4) = %.5fD0\n', values(4));
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pronyData = get_target_curves(csv_path, time, tau)
    % Reads prony coefficients from a CSV file and computes the prony fit of 
    % the creep compliance.
    T = readtable(csv_path, 'ReadRowNames', true);    
    pronyData = struct('name', {}, 'gammas', {}, 'time', {}, 'pronyFit', {});    
    DNames = T.Properties.RowNames;
    
    for i = 1:length(DNames)
        gammas = table2array(T(i, :));       
        pronyFit = zeros(size(time));
        for j = 1:length(gammas)
            pronyFit = pronyFit + gammas(j) * (1- exp(-time ./ tau(j)));
        end
        pronyData(i).name = DNames{i};
        pronyData(i).gammas = gammas;
        pronyData(i).time = time;
        pronyData(i).pronyFit = pronyFit;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pronyData = get_creep_curves(csv_path, time, tau)
    % Reads prony coefficients from a CSV file and computes the prony fit of 
    % the creep compliance.
    %T = readtable(csv_path, 'ReadRowNames', true);
    T = readtable(csv_path,'ReadVariableNames',true,'ReadRowNames',false);
    pronyData = struct('name', {}, 'gammas', {}, 'time', {}, 'pronyFit', {});    
    DNames = size(T,1); %T.Properties.RowNames;
    
    for i = 1:size(T,1)%length(DNames)
        gammas = table2array(T(i, :));       
        pronyFit = zeros(size(time));
        for j = 1:length(gammas)
            pronyFit = pronyFit + 1/gammas(j) * (1- exp(-time ./ tau(j)));
        end
        %pronyData(i).name = DNames{i};
        pronyData(i).gammas = gammas;
        pronyData(i).time = time;
        pronyData(i).pronyFit = pronyFit;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [overallError, overallRelativeError, errors, ...
    relativeErrors] = compareCreepCurves(set1, set2, weights)
    nComponents = length(set1);
    if length(weights) ~= nComponents
        error('Weights vector length must equal the number of components (%d).', nComponents);
    end
    errors = zeros(nComponents, 1);
    relativeErrors = zeros(nComponents, 1);
    
    for i = 1:nComponents
        curve1 = set1(i).pronyFit; %target
        curve2 = set2(i).pronyFit; % simulation

        errors(i) = norm(curve1 - curve2);
        if norm(curve1) > 0
            relativeErrors(i) = errors(i) / norm(curve1);
        else
            relativeErrors(i) = NaN;
        end
    end
    disp(relativeErrors);
    overallError = sum(weights(:) .* errors) / sum(weights);
    overallRelativeError = sum(weights(:) .* relativeErrors) / sum(weights);
end
