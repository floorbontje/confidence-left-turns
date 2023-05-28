% Training_Decision.m 
% By Floor Bontje (2023)

% Creating a decision model for go/wait decisions based on DDM modelling
% with a generalized gap for left-turn gap acceptance. The model has a
% dynamic drift-rate as well as dynamic boundaries. 

% Based on the findings by Zgonnikov, A. et al. (2020, July 7)
% "Should I stay or should I go? Evidence accumulation drives decision
%making in human drivers"  
     
% This code searches for the optimal values of the model parameters 
% for different potential decision models. The difference between the
% potential models is defined by the amount of free parameters, describing 
% the evidence accumulator(s) ("go"/"wait"). 

clear all
clc 

%% Initialize
% Initial parameters: 
% parameters found for the initial drift-diffusion decision model 
parameters_decision_RT_WLSVincent = ...
    readtable('.\Modelling\Parameters Initial\subj_all_RT_WLSVincent.csv');
pos = find(parameters_decision_RT_WLSVincent.loss == min(parameters_decision_RT_WLSVincent.loss));
parameters = parameters_decision_RT_WLSVincent(pos,:);


% Initial parameters
x0(1) = parameters.alpha;           x0(2) = parameters.beta; 
x0(3) = parameters.b_0;             x0(4) = parameters.k; 
x0(5) = parameters.ndt_location;    x0(6) = parameters.ndt_scale; 
x0(7) = parameters.theta;           x0(8:14)= x0; 

%Boundary conditions model 
b0(1,:) = [0.1, 2];         % alpha
b0(2,:) = [0.001, 1];       % beta 
b0(3,:) = [0.5, 3];         % b_0
b0(4,:) = [0.1, 2];         % k 
b0(5,:) = [0, 2.5];         % ndt_location
b0(6,:) = [0.01, 1];        % ndt_scale 
b0(7,:) = [4, 25];          % theta_c

lb0 = [b0(:,1);b0(:,1)]; 
ub0 = [b0(:,2);b0(:,2)]; 

% Storage location fitted parameters 
file_loc= '.\Modelling\decision\Parameters\'; 
    
%% Training 

% Model conditions:
% 0: one decision accumulator for both decisions (DDM)
% 1: alpha differs between go/wait accumulators
% 2: alpha & theta differ between go/wait accumulators
% 3: alpha & b0 differ between go/wait accumulators
% 4: all parameters differ between go/wait accumulators 

trials = 1000; 
x_opt = x0; 
txt = '_int'; 

for tr = 1:2
     if tr == 1 
        %based on initial training 
        load ([file_loc,'2DV_decision_x0.mat'])
        load ([file_loc,'2DV_decision_para_alpha_int.mat']) 
        load ([file_loc,'2DV_decision_para_alpha_b0_int.mat'])
        load ([file_loc,'2DV_decision_para_alpha_theta_int.mat'])
        load ([file_loc,'2DV_decision_para_all_int.mat'])
         
    else 
        % based on x0
        load ([file_loc,'2DV_decision_para_alpha_x0.mat']) 
        load ([file_loc,'2DV_decision_para_alpha_b0_x0.mat'])
        load ([file_loc,'2DV_decision_para_alpha_theta_x0.mat'])
        load ([file_loc,'2DV_decision_para_all_x0.mat'])
     end 
     
for kind_of_decision = 1:4 
      
    if kind_of_decision == 1 
        fun = @(x)function_Decision_Model(trials, [x(1,1:8),x(1,2:7)], 0, 'Alpha');
        
        for i = 1 :10 
            x_alpha (i,:) = fmincon(fun,x_opt(1:8), [],[],[],[],lb0(1:8),ub0(1:8)); 
            loss_alpha(i) = function_Decision_Model(trials, [x_alpha(i,1:8),x_alpha(i,2:7)],0, 'Alpha'); 
        end
        pos_min_alpha = find(loss_alpha == min(loss_alpha));
        save([file_loc,'2DV_decision_para_alpha', txt, '.mat'],'x_alpha','loss_alpha', 'pos_min_alpha')
        
    elseif kind_of_decision == 2 
         
         fun = @(x)function_Decision_Model(trials, [x(1,1:8),x(1,2:6),x(1,9)], 0, 'Alpha & theta_c');
        
        for i = 1 :10 
            x_alpha_theta (i,:) = fmincon(fun,[x_opt(1:8),x_opt(14)], [],[],[],[],[lb0(1:8);lb0(14)],[ub0(1:8);ub0(14)]); 
            loss_alpha_theta(i) = function_Decision_Model(trials, [x_alpha_theta(i,1:8),x_alpha_theta(i,2:6),x_alpha_theta(i,9)],0, 'Alpha & theta_c'); 
        end
        pos_min_alpha_theta = find(loss_alpha_theta == min(loss_alpha_theta));
        save([file_loc,'2DV_decision_para_alpha_theta', txt, '.mat'],'x_alpha_theta','loss_alpha_theta', 'pos_min_alpha_theta')
         
    elseif kind_of_decision == 3 
        
          fun = @(x)function_Decision_Model(trials, [x(1,1:8),x(1,2),x(1,9),x(1,4:7)], 0, 'Alpha & b0');
        
        for i = 1 :10 
            x_alpha_b0 (i,:) = fmincon(fun,[x_opt(1:8),x_opt(10)], [],[],[],[],[lb0(1:8);lb0(10)],[ub0(1:8);lb0(10)]); 
            loss_alpha_b0(i) = function_Decision_Model(trials, [x_alpha_b0(i,1:8),x_alpha_b0(i,2),x_alpha_b0(i,9),x_alpha_b0(i,4:7)],0, 'Alpha&b0'); 
        end
        pos_min_alpha_b0 = find(loss_alpha_b0 == min(loss_alpha_b0));
        save([file_loc,'2DV_decision_para_alpha_b0', txt, '.mat'],'x_alpha_b0','loss_alpha_b0', 'pos_min_alpha_b0')
    
    elseif kind_of_decision == 4 
          
        fun = @(x)function_Decision_Model(trials, x, 0, 'All');
        
        for i = 1 :10 
            x_all (i,:) = fmincon(fun,x_opt, [],[],[],[],lb0,ub0); 
            loss_all(i) = function_Decision_Model(trials, x_all(i,:),0, 'All'); 
        end
        pos_min_all= find(loss_all == min(loss_all));
        save([file_loc,'2DV_decision_para_all', txt, '.mat'],'x_all','loss_all', 'pos_min_all')

    end 
end 

    if tr == 1
    fun = @(x)function_Decision_Model(trials, [x(1,1:7),x(1,1:7)], 0, 'x0');
    for i = 1 :10 
        x_0 (i,:) = fmincon(fun,x0(1:7), [],[],[],[],lb0(1:7),ub0(1:7)); 
        loss_0(i) = function_Decision_Model(trials, [x_0(i,1:7),x_0(i,1:7)],0, 'x0'); 
    end
    pos_min_0 = find(loss_0 == min(loss_0));
    save([file_loc,'2DV_decision_x0.mat'],'x_0','loss_0', 'pos_min_0')
    x_opt = [x_0(pos_min_0,:),x_0(pos_min_0,:)]; 
    txt = '_x0'; 
    end 
end 

%% Plotting results 
close all 
trials = 1000; 
plot_fig = 0; 

for tr = 1:2
      if tr == 1 
        %based on initial training 
        load ([file_loc,'2DV_decision_x0.mat'])
        load ([file_loc,'2DV_decision_para_alpha_int.mat']) 
        load ([file_loc,'2DV_decision_para_alpha_b0_int.mat'])
        load ([file_loc,'2DV_decision_para_alpha_theta_int.mat'])
        load ([file_loc,'2DV_decision_para_all_int.mat'])
        txt = '(int)';
    else 
        % based on x0
        load ([file_loc,'2DV_decision_para_alpha_x0.mat']) 
        load ([file_loc,'2DV_decision_para_alpha_b0_x0.mat'])
        load ([file_loc,'2DV_decision_para_alpha_theta_x0.mat'])
        load ([file_loc,'2DV_decision_para_all_x0.mat'])
        txt = '(x0)';
      end 
     

    for i = 1:10 
        if tr == 1 
            WLS_int(i)  = function_Decision_Model(trials, x0,plot_fig, 'intital input'); 
            WLS_x0(i)   = function_Decision_Model(trials, [x_0(pos_min_0,1:7),x_0(pos_min_0,1:7)],plot_fig, 'x0'); 
        end 
            WLS_x_alpha(i) = function_Decision_Model(trials, [x_alpha(pos_min_alpha,1:8),x_alpha(pos_min_alpha,2:7)],plot_fig, ['Alpha', txt]); 
            WLS_x_alpha_b0(i) = function_Decision_Model(trials, ...
               [x_alpha_b0(pos_min_alpha_b0,1:8),x_alpha_b0(pos_min_alpha_b0,2),x_alpha_b0(pos_min_alpha_b0,9),x_alpha_b0(pos_min_alpha_b0,4:7)],plot_fig, ['Alpha&b0', txt]); 
            WLS_x_alpha_theta(i) = function_Decision_Model(trials, [x_alpha_theta(pos_min_alpha_theta,1:8),...
                x_alpha_theta(pos_min_alpha_theta,2:6),x_alpha_theta(pos_min_alpha_theta,9)],plot_fig, ['Alpha&theta', txt]);
            WLS_x_all(i) = function_Decision_Model(trials, x_all(pos_min_all,:),plot_fig, ['All', txt]);
    end 

    WLS = table(mean(WLS_int), 1.96*std(WLS_int)/10, mean(WLS_x0), 1.96* std(WLS_x0)/10, mean(WLS_x_alpha), 1.95*std(WLS_x_alpha)/10,...
    mean(WLS_x_alpha_b0), 1.96*std(WLS_x_alpha_b0)/10, mean(WLS_x_alpha_theta), 1.96*std(WLS_x_alpha_theta)/10,...
    mean(WLS_x_all), 1.96*std(WLS_x_all)/10, 'VariableNames',...
    {'WLS int', 'WLS int CI', 'WLS mean x0', 'WLS CI x0', ...
    'WLS mean alpha',  'WLS CI alpha',   'WLS mean alpha b0',  'WLS CI alpha b0',... 
    'WLS mean alpha theta',  'WLS CI alpha theta', 'WLS mean all',    'WLS CI all'});

    if tr == 1 
        WLS_intitial = WLS 
        writetable(WLS, [file_loc, '\WLS_best_fit_int.csv'])
    else 
        WLS_x0 = WLS
        writetable(WLS, [file_loc, '\WLS_best_fit_x0.csv'])
    end 

end 
