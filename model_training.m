%%  Empty environment variables
warning off             % turn off warning messages
close all               % Close open windows
clear                   % Clear variables
clc                     % clear command line


%%  Add Path
 addpath('RF_Toolbox\');

%%  retrieve data
res = xlsread('feature_20.csv');%Requires you to complete data feature extraction to save to .cvs
res=res([1:8113],:);%Pick the data you think is appropriate
%%  Analyzing data

num_class = length(unique(res(:, end)));  % number of categories (Excel last column to put categories)
num_res = size(res, 1);                   % number of samples (each row, is a sample)
num_size = 0.7;                           % proportion of training set to data set
res = res(randperm(num_res), :);          % disrupt the dataset (comment the row when not disrupting the data)
f_ = size(res, 2) - 1;                    % number of features
flag_conusion = 1;                        % lag bit is 1 to open confusion matrix (requires version 2018 and above)

%%  Setting variables to store data
P_train = []; P_test = [];
T_train = []; T_test = [];

%%  Divide the data set
for i = 1 : num_class
    mid_res = res((res(:, end) == i), :);           % Loop through the samples of different categories
    mid_size = size(mid_res, 1);                    % get the number of samples of different categories
    mid_tiran = round(num_size * mid_size);         % get the number of training samples in the category

    P_train = [P_train; mid_res(1: mid_tiran, 1: end - 1)];       % Input for training set
    T_train = [T_train; mid_res(1: mid_tiran, end)];              % Training set output

    P_test  = [P_test; mid_res(mid_tiran + 1: end, 1: end - 1)];  % Test set input
    T_test  = [T_test; mid_res(mid_tiran + 1: end, end)];         % test set output
end

%%   data transposition
P_train = P_train'; P_test = P_test';
T_train = T_train'; T_test = T_test';

%%  Get the training set and number of test samples
M = size(P_train, 2);
N = size(P_test , 2);

%% data normalization
[p_train, ps_input] = mapminmax(P_train, 0, 1);
p_test = mapminmax('apply', P_test, ps_input );
t_train = T_train;
t_test  = T_test ;

%%  Transpose to fit the model
p_train = p_train'; p_test = p_test';
t_train = t_train'; t_test = t_test';

%%  Parameter initialization
c1      =  1.49;            % learning factor 
c2      =  1.49;            % learning factor
maxgen  =     50;             % Number of population renewals 
sizepop =     5;             % population size
Vmax    = [ 100.0,  2.0 , 1];    % Maximum speed
Vmin    = [-100.0, -2.0 , -1];   % Minimum velocity
popmax  = [  600, f_ , 10] ;      % Maximum boundary
popmin  = [  100,  1  ,1];       %  Minimum boundary

%%   population initialization
for i = 1 : sizepop
    pop(i, :) = rand(1, length(popmax)) .* (popmax - popmin) + popmin;% Initialize the particle's position
    V(i, :) = rand(1, length(Vmax)) .* Vmax + 1;                      
    fitness(i) = fun(pop(i, :), p_train, t_train);
end
%%  Individual best and global best
[fitnesszbest, bestindex] = min(fitness);
zbest = pop(bestindex, :);     %  global best
gbest = pop;                    % individual best
fitnessgbest = fitness;         % individual best fit value
BestFit = fitnesszbest;          % global best fitness value
w = [0.6 1.0];                  % inertia weights 
alpha = 1-1./maxgen;             % global best fitness value
tic
%%  Iterative optimization
for i = 1 : maxgen
    for j = 1 : sizepop
       W = w(2)-(1 - 1./(1+i).^alpha)*(w(2)-w(1));
        V(j, :) =W*V(j, :) + c1  * rand * (gbest(j, :) - pop(j, :)) + c2 * rand * (zbest - pop(j, :));%rand 0-1
        if (sum(V(j, :) > Vmax) > 0 || sum(V(j, :) < Vmin) > 0)% Over-speed and slow processing
            [~, vmax_index] = find(V(j, :) > Vmax);
            [~, vmin_index] = find(V(j, :) < Vmin);
            V(j, vmax_index) = Vmax(vmax_index);
            V(j, vmin_index) = Vmin(vmin_index);
        end
        
        % Population renewal and location renewal
        pop(j, :) = pop(j, :) + 0.2 * V(j, :);

        if (sum(pop(j, :) > popmax) > 0 || sum(pop(j, :) < popmin) > 0)% border processing
            [~, pmax_index] = find(pop(j, :) > popmax);
            [~, pmin_index] = find(pop(j, :) < popmin);
            pop(j, pmax_index) = popmax(pmax_index);
            pop(j, pmin_index) = popmin(pmin_index);
        end
        
        % adaptive variation 
        if rand > 0.95
            pop(j, :) = rand(1, length(popmax)) .* (popmax - popmin) + 1;
        end
        
        % fitness value
        fitness(j) = fun(pop(j, :), p_train, t_train); %Calculate adaptation once after updating
    end
    
    for j = 1 : sizepop
        
        % Individual best Update
        if fitness(j) < fitnessgbest(j)
            gbest(j, :) = pop(j, :);
            fitnessgbest(j) = fitness(j);
        end

        % global best Update 
        if fitness(j) < fitnesszbest
            zbest = pop(j, :);
            fitnesszbest = fitness(j);
        end

    end
   if mod(i,2)  == 0
        i                  % Output convergence progress
    end
    pause(0.1) 
    train_fit(i)=max(fitnesszbest);
    valid_fit(i)=mean(fitnesszbest);
    BestFit = [BestFit, fitnesszbest];    
end
toc
%%   Extraction of optimal parameters
n_trees = round(zbest(1));
n_layer = round(zbest(2));
nodesize = round(zbest(3));
%% Create the model
model = classRF_train(p_train, t_train, n_trees, n_layer, nodesize); % importance = model.importance; % importance of features
importance = model.importance; %% importance of features

%% Simulation test
[T_sim1, Vote1] = classRF_predict(p_train, model);
tic
[T_sim2, Vote2] = classRF_predict(p_test , model);
toc
%%  performance evaluation
error1 = sum((T_sim1' == T_train)) / M * 100 ;
error2 = sum((T_sim2' == T_test ))  / N * 100 ;

%% plot
figure
plot(1: M, T_train, 'r-*', 1: M, T_sim1, 'b-o', 'LineWidth', 1)
legend('True label', 'Predicted label')
xlabel('Predicted sample')
ylabel('Predicted Result')
string = {'Comparison of training set prediction results'; ['Accuracy=' num2str(error1) '%']};
title(string)
grid

figure
plot(1: N, T_test, 'r-*', 1: N, T_sim2, 'b-o', 'LineWidth', 1)
legend('True label', 'Predicted label')
xlabel('Predicted sample')
ylabel('Predicted Result')
string = {'Comparison of test set prediction results'; ['Accuracy=' num2str(error2) '%']};
title(string)
grid

%% Iterative plot of error curves
figure;
plot(1 : length(BestFit), BestFit, 'LineWidth', 1.5);
xlabel('Number of particle swarm iterations');
ylabel('Fitness value');
xlim([1, length(BestFit)])
string = {'Model iteration error variation'};
title(string)
grid on

%% Mapping the significance of features
figure
bar(importance)
legend('importance')
xlabel('feature')
ylabel('importance')

%%  confusion matrix
if flag_conusion == 1

    figure
    cm = confusionchart(T_train, T_sim1);
    cm.Title = 'Confusion Matrix for Train Data';
    cm.ColumnSummary = 'column-normalized';
    cm.RowSummary = 'row-normalized';
    
    figure
    cm = confusionchart(T_test, T_sim2);
    cm.Title = 'Confusion Matrix for Test Data';
    cm.ColumnSummary = 'column-normalized';
    cm.RowSummary = 'row-normalized';
end
%%
save('.../your_own_model.mat','model'); 
