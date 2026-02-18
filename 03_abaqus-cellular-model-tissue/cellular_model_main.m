clear, clc;

% Set working directory
cwd = fileparts(mfilename('fullpath'));
%{
% Define pathsto delete
tasks = {
    'simulation_lock.lock'
    'Results.mat'
    'Job-*'
    'abaqus*'
    'ABQ*'
    fullfile('TW','TW_gamma.csv')
    fullfile('TW','TW_creep_compliance_coeffs.csv')
    fullfile('EW','EW_gamma.csv')
    fullfile('EW','EW_creep_compliance_coeffs.csv')
    fullfile('LW','LW_gamma.csv')
    fullfile('LW','LW_creep_compliance_coeffs.csv')
    fullfile('GR','GR_gamma*.csv')
    fullfile('GR','GR_creep_compliance_coeffs.csv')
    fullfile('Modules','*.inc')
    fullfile('Modules','layers_gamma_attempts.csv')};

% Delete files
for i = 1:numel(tasks)
    pattern = tasks{i};
    files = dir(fullfile(cwd, pattern));
    for f = files'
        fullpath = fullfile(f.folder, f.name);
        delete(fullpath);
        fprintf('File "%s" found and deleted.\n', fullpath);
    end
end
%}

% Delete existing attempts
folders = dir(fullfile(cwd,'Attempt_*'));
for k = 1:numel(folders)
    if folders(k).isdir
        fullP = fullfile(cwd, folders(k).name);
        rmdir(fullP,'s');                    % remove folder + contents
        fprintf('Deleted old attempt folder: %s\n', folders(k).name);
    end
end

% Add MATSuMoTo Toolbox to path
addpath(genpath(fullfile(cwd,'MATSuMoTo-master')));
% Define MATSuMoTo input parameters
data_file = 'cellular_model_datainput'; % file containing the problem setup
maxeval = 50; % maximum number of function evaluations
surogate_model = 'POLYcubr'; % surrogate model
sampling_technique = 'CANDglob'; % sampling strategy for selecting new candidate points
initial_design = 'LHS'; % type of initial experimental design
number_startpoints = []; % no extra start points are provided
NumberNewSamples = 7; % new sample points to add at each iteration of optimization
% Read initial values from csv file
csv_filename = fullfile(cwd,'Original_Folder','Modules','tracheids_gamma_initial_values.csv');
T = readtable(csv_filename);
lastRow = T(end, :);
starting_point = [lastRow.gamma1, lastRow.gamma2, lastRow.gamma3, lastRow.gamma4];

% Limit parallel pool size
desiredWorkers = NumberNewSamples;
p = gcp('nocreate');
if isempty(p)
    parpool('local', desiredWorkers);
elseif p.NumWorkers ~= desiredWorkers
    delete(p);
    parpool('local', desiredWorkers);
end

% Run MATSuMoTo
MATSuMoTo(data_file,maxeval,surogate_model,sampling_technique,...
    initial_design,number_startpoints,starting_point,NumberNewSamples);
