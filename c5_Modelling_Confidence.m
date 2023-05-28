% Modelling_Confidence.m
% By Floor Bontje 

% Four potential confidenc models are evaluated, all models are build on a
% decision model. Two different decision models are taken into accoun the
% race model(RT) and the drift-diffusion model (DDM). Moreover, the
% influence of post decision evidence accumulation on confidence is
% investigated. 

close all
clear all 
clc 

file_loc = '.\Modelling\decision\Parameters\';
file_loc_conf = '.\Modelling\Confidence\'; 

%% Decision models
close all; clc 
trials = 1000;  plot_fig = 1; 


% Initial parameters ("go" decisions)
parameters_decision_RT_WLSVincent = ...
    readtable('.\Modelling\Parameters Initial\subj_all_RT_WLSVincent.csv');
pos = find(parameters_decision_RT_WLSVincent.loss == min(parameters_decision_RT_WLSVincent.loss));
parameters = parameters_decision_RT_WLSVincent(pos,:);

x_int(1) = parameters.alpha;           x_int(2) = parameters.beta; 
x_int(3) = parameters.b_0;             x_int(4) = parameters.k; 
x_int(5) = parameters.ndt_location;    x_int(6) = parameters.ndt_scale; 
x_int(7) = parameters.theta;           x_int(8:14)= x_int(1:7); 


% Trained parameters  ("Training_Decision.m") 
load([file_loc,'2DV_decision_x0.mat'])
p_DDM = x_0(pos_min_0,:); 

load([file_loc,'2DV_decision_para_alpha_theta_x0.mat'])
p_RM = [x_alpha_theta(pos_min_alpha_theta,1:8),...
            x_alpha_theta(pos_min_alpha_theta,2:6),x_alpha_theta(pos_min_alpha_theta,9)];


% DDM, initial parameters  
WLS_int_DDM = function_Decision_Model(trials, x_int,plot_fig, 'int');  

% DDM, trained paramters 
WLS_DDM = function_Decision_Model(trials, [p_DDM,p_DDM],plot_fig, 'DDM'); 

% Race model (alpha & theta)
WLS_RM =function_Decision_Model(trials, p_RM, plot_fig, 'Alpha and theta');



%% Confidence models with-out additional evidence accumulation  
trials = 1000; 
plot_fig = 0; 

% Initial confidence parameters 
[coeff_go_m1_0, coeff_wait_m1_0] = function_Conf_DDM_x0(trials, p_DDM, 0); 
[coeff_go_m3_0, coeff_wait_m3_0] = function_Conf_Race_Model_x0(trials, p_RM, 0); 

 
% Training coefficients 
fun_m1 = @(x)function_Conf_DDM(trials, p_DDM, 0, plot_fig, x(1:2), x(3:4)); 
fun_m3 = @(x)function_Conf_Race_Model(trials, p_RM, 0, plot_fig, x(1:2), x(3:4));
 

for i = 1:10
        coeff_m1_fit(i,:) = fmincon(fun_m1,[coeff_go_m1_0;coeff_wait_m1_0], [],[],[],[],[-10, 0, -10, 0],10*ones(4,1))'; 
        RMSE_m1(i) = function_Conf_DDM(trials, p_DDM, 0, plot_fig, coeff_m1_fit(i,1:2), coeff_m1_fit(i,3:4)); 
           
        coeff_m3_fit(i,:) = fmincon(fun_m3,[coeff_go_m3_0;coeff_wait_m3_0], [],[],[],[],[-10,0,-10,0],10*ones(4,1))'; 
        RMSE_m3(i) = function_Conf_Race_Model(trials, p_RM, 0, plot_fig, coeff_m3_fit(i,1:2), coeff_m3_fit(i,3:4));
end 

pos_m1 = find(RMSE_m1 == min(RMSE_m1));
pos_m3 = find(RMSE_m3 == min(RMSE_m3)); 

coeff_go_m1   = coeff_m1_fit(pos_m1,1:2); 
coeff_wait_m1 = coeff_m1_fit(pos_m1,3:4);

coeff_go_m3   = coeff_m3_fit(pos_m3,1:2); 
coeff_wait_m3 = coeff_m3_fit(pos_m3,3:4);


save ([file_loc_conf,'RMSE_models.mat'],  'RMSE_m1', 'RMSE_m3')
save ([file_loc_conf,'Coeff_m_fit.mat'],  'coeff_m1_fit','coeff_m3_fit') 
save ([file_loc_conf,'Coeff_conf_model.mat'],'coeff_go_m1','coeff_wait_m1','coeff_go_m3','coeff_wait_m3')

%% Determination of extra evidence accumulation time (tau)
% Plot of effect CT values for different values of tau 
% tau is the time evidence accumulaion continues after the decision is made

t = 0.0:0.01:2.5; 
plot_fig=0; 
for ii = 1:length(t)
    tau = t(ii); 
  
    [coeff_go_m2_ct, coeff_wait_m2_ct] = function_Conf_DDM_x0(trials, p_DDM, tau);
    [coeff_go_m4_ct, coeff_wait_m4_ct] = function_Conf_Race_Model_x0(trials, p_RM);
    
    
    if coeff_wait_m2_ct(2) < 0 
        coeff_wait_m2_ct(2) = 0; 
    end 
    if coeff_wait_m4_ct(2) < 0
        coeff_wait_m4_ct(2) = 0; 
    end 
    coeff_m2_go (ii,:)=coeff_go_m2_ct'; coeff_m2_wait (ii,:)=coeff_wait_m2_ct'; 
    coeff_m4_go (ii,:)=coeff_go_m4_ct'; coeff_m4_wait (ii,:)=coeff_wait_m4_ct'; 
    for i = 1:5
        Loss_m2 (ii,i)= function_Conf_DMM(trials, p_DDM, tau, 0, coeff_go_m2_ct,coeff_wait_m2_ct); 
        Loss_m4 (ii,i)= function_Conf_Race_Model(trials, p_RM, tau, 0, coeff_go_m4_ct, coeff_wait_m4_ct);
    end 
    Loss_conf_model2_t(ii)= mean(Loss_m2(ii,:)); 
    Loss_conf_model4_t(ii)= mean(Loss_m4(ii,:));
end 
save('CT_mean_RMSE.mat','t','Loss_conf_model4_t', 'Loss_conf_model2_t','coeff_m2_go', 'coeff_m4_go', 'coeff_m2_wait', 'coeff_m4_wait')


%% Figure: Effect of additional time (tau) RMSE 
load([file_loc_conf,'CT_mean_RMSE.mat'])

CT_add_time_fig =figure; 
plot(t, Loss_conf_model2_t,'LineWidth', 1.3)
hold on; grid on; box off
plot(t, Loss_conf_model4_t,'LineWidth', 1.3)
plot(t(t==1.2),Loss_conf_model4_t(t==1.2),'.k','MarkerSize', 20)
plot(t(t==1.2),Loss_conf_model2_t(t==1.2),'.k','MarkerSize', 20)
set(gca,'FontSize',13)
legend('DDM (model 2)', 'Race model (model 4)') 
xlabel('Inter-judgement time (\tau), s')
ylabel('RMSE')
% title ('Effect Different Additional Times (CT = RT + \tau))')

CT_time_point_fig = figure; 
plot(t(81:181),Loss_conf_model2_t(81:181),'LineWidth', 1.3)
hold on; grid on; box off
plot(t(81:181),Loss_conf_model4_t(81:181),'LineWidth', 1.3)
plot(t(t==1.2),Loss_conf_model4_t(t==1.2),'.k','MarkerSize', 20)
plot(t(t==1.2),Loss_conf_model2_t(t==1.2),'.k','MarkerSize', 20)
xlim([t(81),t(181)])
set(gca,'FontSize',13)
legend('DDM (model 2)','Race model (model 4)') 
xlabel('Inter-judgement time (\tau), s') 
ylabel('RMSE') 
% title ('Effect Different Additional Times (CT = RT + \tau)')


saveas(CT_add_time_fig, fullfile(file_loc_conf,'Effect additional time CT_version_lm.jpg'))
saveas(CT_time_point_fig, fullfile(file_loc_conf, 'Effect additional time around tau.jpg')); 



%% Confidence models
load([file_loc_conf,'Coeff_conf_model.mat'])
load([file_loc_conf,'CT_mean_RMSE.mat'])
tau_opt = 1.2; %Defined additional time  
t = 0.0:0.01:2.5; 
i = find(t==tau_opt);

for count = 10
    plot_fig = 0;
    if count == 10 
        plot_fig = 1; 
    end 
    RMSE_m1_fit(count)=function_Conf_DDM(trials, p_DDM, 0, plot_fig, coeff_go_m1,coeff_wait_m1); 
    RMSE_m2_fit(count)=function_Conf_DDM(trials, p_DDM, t(i), plot_fig, coeff_m2_go(i,:),coeff_m2_wait(i,:)); 
    RMSE_m3_fit(count)=function_Conf_Race_Model(trials, p_RM, 0, plot_fig, coeff_go_m3,coeff_wait_m3);
    RMSE_m4_fit(count)=function_Conf_Race_Model(trials, p_RM, t(i), plot_fig, coeff_m4_go(i,:),coeff_m4_wait(i,:));
end 

mean_RMSE_m1 = mean(RMSE_m1_fit); 
mean_RMSE_m2 = mean(RMSE_m2_fit); 
mean_RMSE_m3 = mean(RMSE_m3_fit); 
mean_RMSE_m4 = mean(RMSE_m4_fit); 


Loss_mean_models = table(mean_RMSE_m1,mean_RMSE_m2,mean_RMSE_m3,mean_RMSE_m4,'VariableNames', ... 
    {'Model 1', 'Model2','Model3', 'Model4'})

writetable(Loss_mean_models, fullfile(file_loc_conf, 'RMSD_models_mean.csv'))





