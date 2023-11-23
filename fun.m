function fitness = fun(pop, p_train, t_train)

%%  Extraction of optimal parameters
n_trees = round(pop(1));
n_layer = round(pop(2));
nodesize = round(pop(3));
%%  Parameters of the data
num_size = length(t_train);

%%  Cross-validation procedures
indices = crossvalind('Kfold', num_size, 7);

for i = 1 : 7
    
    % Get the logical value of the index for the ith data copy
    valid_data = (indices == i);
    
    % Invert to get the logical value of the index of the ith training data
    train_data = ~valid_data;
    
    % 1 test, 7 trainings
    pv_train = p_train(train_data, :);
    tv_train = t_train(train_data, :);
    
    pv_valid = p_train(valid_data, :);
    tv_valid = t_train(valid_data, :);

   % Build the model
    model = classRF_train (pv_train, tv_train, n_trees, n_layer, nodesize);

    % Simulation test
    t_sim = classRF_predict (pv_valid, model);

    % fitness value
    error(i) = 1 - sum(t_sim == tv_valid) ./ length(tv_valid);

end

%% Acquisition of fitness
fitness = mean(error);