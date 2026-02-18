function Data = cellular_model_datainput   
    Data.xlow = [10, 0, 10, 0.]; % lower variable bounds
    Data.xup = [30, 10, 30, 10]; % upper variable bounds
    Data.dim = 4; % problem dimension = 4 gamma_i
    Data.integer = []; %indices of integer variables
    Data.continuous = (1:4); % indices of continuous variables   
    Data.objfunction=@cellular_model_run_simulations; % handle to objective function
end