function  WLS = function_Decision_Model(trials, x, plot_fig, txt)
% Decision model with dynamic drift-rate and dynamic boundaries

% Weighted least squares method (Ratcliff & Tuerlinckx (2002):
% The sum of the squared difference between observed an predicted acccuracy values plus the 
% sum of the squared differences between observed predicted quantile reaction times for 
% correct and error responses is minimized. 

% If all the parameters (x) are equal the model behaves as a
% drift-diffusion model. Otherwise, does the model describe a race model. 

% Model parameters 
alpha1      = x(1);               beta1     = x(2); 
b01         = x(3);               k1        = x(4); 
mu_ND1      = x(5);               sigma_ND1 = x(6); 
theta_c1    = x(7);

alpha2      = x(8);               beta2     = x(9); 
b02         = x(10);              k2        = x(11); 
mu_ND2      = x(12);              sigma_ND2 = x(13); 
theta_c2    = x(14);


%% Initialise
% Conditions; dynamic drift depending on distance and TTA and dynamic decision boundary with collapses with TTA
tta  = [5.5,6.5]; 
dist = [70,90];

% Time 
max_t = 7;                              % Maximum time
delta_t   = 0.001;                      % time steps 
max_inter = max_t/delta_t;              % total interations
t = (0:delta_t:max_t-delta_t)';         % time vector


% Noise
accu_noise_mean = 0; accu_noise_std = 1; 
sen_noise_mean  = 0;  sen_noise_std = 0; 
 
% Evidence 
evidence_go     = nan(max_inter,trials*length(tta)*length(dist));  
evidence_wait   = nan(max_inter,trials*length(tta)*length(dist));  
decision_go     = zeros (1,trials*length(tta)*length(dist)); 
decision_wait   = zeros (1,trials*length(tta)*length(dist)); 
RT_raw_go       = ones(1,4*trials)*1000;  
RT_raw_wait     = ones(1,4*trials)*1000;                


% Non-decision time  
t_noise_go      = sigma_ND1*randn(1, trials*length(tta)*length(dist)); 
t_noise_wait    = sigma_ND2*randn(1, trials*length(tta)*length(dist)); 

if sigma_ND2 == sigma_ND1 
    t_noise_wait = t_noise_go; 
end 

t_non_decision_go   = t_noise_go+mu_ND1;
t_non_decision_go(t_non_decision_go < 0 ) = 0; 
t_non_decision_wait = t_noise_wait+mu_ND2;
t_non_decision_wait(t_non_decision_wait < 0 ) = 0; 


%% Model
% Boundaries 
tta_con     = [tta,tta];  
dist_con    = [dist(1)*ones(1,2),dist(2)*ones(1,2)]; 

TTA = tta_con - t;                    %time to arrival 
d   = dist_con - dist_con ./tta_con.*t; %distance 
generalized_gap_go      = (TTA + beta1 * d);   %generalized gap 
generalized_gap_wait    = (TTA + beta2 * d);

upperboundary = b01 ./ (1 + exp(-k1 * (generalized_gap_go - theta_c1))); 
lowerboundary = b02 ./ (1 + exp(-k2 * (generalized_gap_wait - theta_c2))); 
    
% noise 
 sen_noise   = 1 - (sen_noise_std*randn(max_inter, trials*length(tta)*length(dist)) + sen_noise_mean); %sensory noise 


if x(1:7) == x(8:14)
    % disp("DDM")
    accu_noise  = accu_noise_std*randn(max_inter, (trials*length(tta)*length(dist))) + accu_noise_mean;  %accumulation noise  
    accu_noise_go = accu_noise; 
    accu_noise_wait = accu_noise; 
else 
    % disp("Race Model")
    accu_noise_go  = accu_noise_std*randn(max_inter, (trials*length(tta)*length(dist))) + accu_noise_mean;  %accumulation noise  
    accu_noise_wait  = accu_noise_std*randn(max_inter, (trials*length(tta)*length(dist))) + accu_noise_mean;  %accumulation noise  

end 

% Drift-rate 
driftrate_c_go      = alpha1.*(generalized_gap_go - theta_c1); 
driftrate_c_wait    = alpha2.*(generalized_gap_wait - theta_c2); 



%% Calculations
driftrate_go    = driftrate_c_go;
driftrate_wait  = driftrate_c_wait; 
tta_con_trials  = tta_con; 
dist_con_trials = dist_con; 
low_boundary_trials     = lowerboundary; 
upper_boundary_trials   = upperboundary;

for i = 1 : trials - 1 
    driftrate_go(:,end+1:end+4)     = driftrate_c_go;
    driftrate_wait(:,end+1:end+4)   = driftrate_c_wait;
    
    tta_con_trials (:,end+1:end+4)       = tta_con; 
    dist_con_trials (:,end+1:end+4)      = dist_con; 
    low_boundary_trials(:,end+1:end+4)   = lowerboundary; 
    upper_boundary_trials(:,end+1:end+4) = upperboundary;
end 


for i = 1:max_inter
        dx_go = sen_noise(i,:).*driftrate_go(i,:)*delta_t + accu_noise_go(i,:)*sqrt(delta_t); 
        dx_wait = -(sen_noise(i,:).*driftrate_wait(i,:)*delta_t + accu_noise_wait(i,:)*sqrt(delta_t)); 
        if i == 1
            evidence_go(i,:)    = dx_go;
            evidence_wait(i,:)  = dx_wait;
        else 
            evidence_go(i,:)   = evidence_go(i-1,:) + dx_go; 
            evidence_wait(i,:) = evidence_wait(i-1,:) + dx_wait; 
        end
        
        
        RT_raw_go(evidence_go(i,:)> upper_boundary_trials(i,:) & decision_go == 0) = t(i);
        decision_go(evidence_go(i,:)> upper_boundary_trials(i,:) & decision_go == 0) = 1; 

        RT_raw_wait(evidence_wait(i,:)> low_boundary_trials(i,:) & decision_wait == 0) = t(i);
        decision_wait(evidence_wait(i,:)> low_boundary_trials(i,:) & decision_wait == 0) = -1; 
end 

% Correction of reaction time, taking into account non decision time 
RT_go = RT_raw_go + t_non_decision_go; 
RT_wait = RT_raw_wait + t_non_decision_wait; 

decision = nan(size(RT_raw_wait)); 
decision(RT_raw_go < RT_raw_wait) = 1;
decision(RT_raw_go > RT_raw_wait) = -1; 
decision(RT_raw_go == RT_raw_wait) = nan; 
prs_no_decision = length(find(isnan(decision)))/length(decision)*100; 

%% Statistics
con = ['c1';'c2';'c3';'c4']; 
for i = 1 :length(tta)*length(dist)
    pos_go.(con(i,:))   = find(decision == 1 & tta_con_trials == tta_con(i) & dist_con_trials == dist_con(i)); 
    pos_wait.(con(i,:)) = find(decision == -1 & tta_con_trials == tta_con(i) & dist_con_trials == dist_con(i));    
    
    
    pos_go_m = zeros(1,trials); 
    pos_go_m(pos_go.(con(i,:)))=1; 
    
    pos_wait_m = zeros(1,trials); 
    pos_wait_m(pos_wait.(con(i,:)))= 1; 
    
    CI_go_model(i) = 1.96 * std(pos_go_m)/sqrt(trials); 
    CI_wait_model(i) = 1.96 * std(pos_wait_m)/sqrt(trials); 
    
    
    go_pr(i)    = length(pos_go.(con(i,:)))/trials;
    wait_pr(i)  = length(pos_wait.(con(i,:)))/trials;
    
    RT_mean_go(i)   = mean(RT_go(pos_go.(con(i,:))));
    RT_mean_wait(i) = mean(RT_wait(pos_wait.(con(i,:)))); 
    
    CI_RT_go(i)   = 1.96 * std(RT_go(pos_go.(con(i,:))))/sqrt(length(pos_go.(con(i,:))));
    CI_RT_wait(i) = 1.96 * std(RT_wait(pos_wait.(con(i,:))))/sqrt(length(pos_wait.(con(i,:))));
end


%% Measures 
F_inter = 0:0.01:1;
data = readtable('.\Modelling\data model feeting\data_output_RT.csv');
for i = 1: 4
    pos_go_exp.(con(i,:)) = find(data.num_decision ==1 & data.tta_condition == tta_con(i) & data.d_condition == dist_con(i)); 
    pos_wait_exp.(con(i,:)) = find(data.num_decision ==2 & data.tta_condition == tta_con(i) & data.d_condition == dist_con(i)); 
    length_con = length(find(data.tta_condition == tta_con(i) & data.d_condition == dist_con(i))); 
    
    go_pr_exp(i) = length(pos_go_exp.(con(i,:)))/length_con; 
    wait_pr_exp(i)= length(pos_wait_exp.(con(i,:)))/length_con; 
    
    pos_go_con = zeros(1,length_con); 
    pos_go_con(pos_go_exp.(con(i,:)))=1; 
    
    pos_wait_con = zeros(1,length_con); 
    pos_wait_con(pos_wait_exp.(con(i,:)))= 1; 
    
    CI_go_exp(i) = 1.96 * std(pos_go_con)/sqrt(length_con); 
    CI_wait_exp(i) = 1.96 * std(pos_wait_con)/sqrt(length_con); 
    
    mean_rt_exp_go(i) = mean(data.RT(pos_go_exp.(con(i,:)))); 
    mean_rt_exp_wait(i) = mean(data.RT(pos_wait_exp.(con(i,:)))); 
    
    CI_RT_go_exp(i) = 1.96 * std(data.RT(pos_go_exp.(con(i,:))))/sqrt(length(pos_go_exp.(con(i,:))));
    CI_RT_wait_exp(i) = 1.96 * std(data.RT(pos_go_exp.(con(i,:))))/sqrt(length(pos_go_exp.(con(i,:)))); 
end 

RT_exp = data.RT; 

% CDf
for i = 1:4
    [Fgo_exp.(con(i,:)),RTgo_exp.(con(i,:))]=ecdf(RT_exp(pos_go_exp.(con(i,:))));
    [Fwait_exp.(con(i,:)),RTwait_exp.(con(i,:))]=ecdf(RT_exp(pos_wait_exp.(con(i,:))));
end 

% Interpolation CDF 
for i = 1:4 
    RTgo_exp_inter.(con(i,:))   = interp1(Fgo_exp.(con(i,:)),RTgo_exp.(con(i,:)),F_inter);
    RTwait_exp_inter.(con(i,:)) = interp1(Fwait_exp.(con(i,:)),RTwait_exp.(con(i,:)),F_inter);
end 

%% Weighted Least Squares Fitting 
% Quantiles and weights 
q = [0.1, 0.3, 0.5, 0.7, 0.9]; 
wt = [2, 2, 1, 1, 0.5]; 


for i = 1: length(q)
    pos_q(i) = find(F_inter == q(i)); 
end 

%Accuracy
pr_th_go = go_pr;
pr_ex_go = go_pr_exp; 
pr_th_wait = wait_pr;
pr_ex_wait = wait_pr_exp; 

%Quantile reaction time (experiment) [sec]
for i = 1:4
    Q_exp_go.(con(i,:)) = RTgo_exp_inter.(con(i,:))(pos_q);
    Q_exp_wait.(con(i,:)) = RTwait_exp_inter.(con(i,:))(pos_q);
end 

for i = 1:4
    WLS(i) = 4*(pr_th_go(i)-pr_ex_go(i))^2 + ...
        4 * (pr_th_wait(i) -pr_ex_wait(i))^2; 
    
    
    if pr_th_go(i) > 0 && pr_ex_go(i) > 0.001
        [Fgo_model.(con(i,:)),RTgo_model.(con(i,:))]=ecdf(RT_go(pos_go.(con(i,:))));
        RTgo_model_inter.(con(i,:)) = interp1(Fgo_model.(con(i,:)),RTgo_model.(con(i,:)),F_inter);
        Q_th_go.(con(i,:))= RTgo_model_inter.(con(i,:))(pos_q); 

        WLS(i) = WLS(i) +  dot((Q_th_go.(con(i,:))  -  Q_exp_go.(con(i,:))).^2, wt) * pr_ex_go(i);
    end 
    
    if pr_th_wait(i) > 0 && pr_ex_wait(i) > 0.001
        [Fwait_model.(con(i,:)),RTwait_model.(con(i,:))]=ecdf(RT_wait(pos_wait.(con(i,:))));
        RTwait_model_inter.(con(i,:)) = interp1(Fwait_model.(con(i,:)),RTwait_model.(con(i,:)),F_inter);
        Q_th_wait.(con(i,:)) = RTwait_model_inter.(con(i,:))(pos_q);

        WLS(i) = WLS(i) +  dot((Q_th_wait.(con(i,:))  -  Q_exp_wait.(con(i,:))).^2, wt) * pr_ex_wait(i);
    end   
end 

WLS = sum(WLS); 
%% Figures
fname = '.\Modelling\decision\figures';

if plot_fig == 1 
  
    
    % Decision behaviour 
    figure0= figure; 
    subplot(1,2,1)
    errorbar(tta, go_pr(1:2)*100, CI_go_model(1:2)*100, '.-b','LineWidth',1.3,'MarkerSize', 20)
    hold on; box off 
    errorbar(tta, go_pr(3:4)*100, CI_go_model(3:4)*100,'.-r','LineWidth',1.3,'MarkerSize', 20)
    errorbar(tta, go_pr_exp(1:2)*100, CI_go_exp(1:2)*100,'.--','color', [.73 .73 1], 'LineWidth',1.3,'MarkerSize', 20)
    errorbar(tta, go_pr_exp(3:4)*100, CI_go_exp(3:4)*100,'.--','color', [1 .73 .73], 'LineWidth',1.3,'MarkerSize', 20)
    ax = gca;             ax.FontSize = 12;
    title('Go', 'FontSize', 14)
    legend('model: 70m', 'model: 90m', 'exp: 70m', 'exp: 90m', 'Location', 'northwest','FontSize', 10)
    legend boxoff
    xlim([5.3, 6.7]);  xticks([5.5, 6.5]);    ylim([0, 100]) 
    xlabel('Time-to-arrival (TTA), s','FontSize', 13)
    ylabel('Percentage of going','FontSize', 13)
    
    
    subplot(1,2,2)
    errorbar(tta, wait_pr(1:2)*100,CI_wait_model(1:2)*100, '.-b','LineWidth',1.3,'MarkerSize', 20)
    hold on; box off
    errorbar(tta, wait_pr(3:4)*100, CI_wait_model(3:4)*100, '.-r','LineWidth',1.3,'MarkerSize', 20)
    errorbar(tta, wait_pr_exp(1:2)*100, CI_wait_exp(1:2)*100,'.--','color', [.73 .73 1], 'LineWidth',1.3,'MarkerSize', 20)
    errorbar(tta, wait_pr_exp(3:4)*100, CI_wait_exp(3:4)*100,'.--','color', [1 .73 .73], 'LineWidth',1.3,'MarkerSize', 20)
    ax = gca;             ax.FontSize = 12;
    title('Wait', 'FontSize', 14)
    xlim([5.3, 6.7]);  xticks([5.5, 6.5]);    ylim([0, 100])
    xlabel('Time-to-arrival (TTA), s','FontSize', 13)
    ylabel('Percentage of waiting','FontSize', 13)
    saveas(figure0, fullfile(fname, [txt, ' decision behaviour modelled.jpg']))
    saveas(figure0, fullfile(fname, [txt, ' decision behaviour modelled.pdf']))
    
    %Predicted RTs by model 
    figure1= figure; 
    subplot(1,2,1)

    errorbar(tta, RT_mean_go(1:2), CI_RT_go(1:2),'.-b','LineWidth',1.3,'MarkerSize', 20)
    hold on; box off
    errorbar(tta, RT_mean_go(3:4), CI_RT_go(3:4), '.-r','LineWidth',1.3,'MarkerSize', 20)
    errorbar(tta, mean_rt_exp_go(1:2),CI_RT_go_exp(1:2),'.--','color', [.73 .73 1], 'LineWidth',1.3,'MarkerSize', 20)
    errorbar(tta, mean_rt_exp_go(3:4),CI_RT_go_exp(3:4),'.--','color', [1 .73 .73], 'LineWidth',1.3,'MarkerSize', 20)
    ax = gca;             ax.FontSize = 12;
    title('Go','FontSize', 14)
    legend('model: 70m', 'model: 90m', 'exp: 70m', 'exp: 90m', 'Location', 'northwest','FontSize', 10)
    legend boxoff
    xlim([5.3, 6.7]);  xticks([5.5, 6.5]);    ylim([1,3.5])
    xlabel('Time-to-arrival (TTA), s','FontSize', 13)
    ylabel('Response time (RT), s','FontSize', 13)
    
    
    subplot(1,2,2)
    errorbar(tta, RT_mean_wait(1:2), CI_RT_wait(1:2), '.-b','LineWidth',1.3,'MarkerSize', 20)
    hold on; box off 
    errorbar(tta, RT_mean_wait(3:4), CI_RT_wait(3:4),'.-r','LineWidth',1.3,'MarkerSize', 20)
    errorbar(tta, mean_rt_exp_wait(1:2), CI_RT_wait_exp(1:2),'.--','color', [.73 .73 1], 'LineWidth',1.3,'MarkerSize', 20)
    errorbar(tta, mean_rt_exp_wait(3:4), CI_RT_wait_exp(3:4),'.--','color', [1 .73 .73], 'LineWidth',1.3,'MarkerSize', 20)
    ax = gca;             ax.FontSize = 12;
    title('Wait','FontSize', 14)
    xlim([5.3, 6.7]);  xticks([5.5, 6.5]);    ylim([1,3.5])
    xlabel('Time-to-arrival (TTA), s','FontSize', 13)
    ylabel('Response time (RT), s','FontSize', 13)
    
    saveas(figure1, fullfile(fname, [txt, ' RT modelled.jpg']))
    saveas(figure1, fullfile(fname, [txt, ' RT modelled.pdf']))
 
    
    %Illustration of model 
    max_bound = 0.5+max(upper_boundary_trials(1,1:4)); 
    txt2 = 'tau 0';
    for conf = 1:2 
    if decision(3) == 1 
        stop_pos = find(t==RT_raw_go(3)); 
        conf_pos = find(t>RT_raw_go(3)+1,1);
        
    else 
        stop_pos = find(t==RT_raw_wait(3)); 
        conf_pos = find(t>RT_raw_wait(3)+1,1);
    end 
    
    if abs(alpha1-alpha2) > 0  
        figure2=figure;  
        hold on
        plot(t(1:stop_pos), evidence_go(1:stop_pos,3),'-','color',[0.3 0.2 0.9],'LineWidth',1.3)
        plot(t(1:stop_pos), -evidence_wait(1:stop_pos,3),'-','color',[0.32 0.7 0.9],'LineWidth',1.3)
        plot(t,upper_boundary_trials(:,3),'--k','LineWidth',1.3)
        plot(t,-low_boundary_trials(:,3),':k','LineWidth',1.3)
        if decision(3) == 1 
            plot(RT_raw_go(3),upper_boundary_trials(stop_pos,3),'.r','MarkerSize', 20,'LineWidth',1.3)
        else 
            plot(RT_raw_wait(3),-low_boundary_trials(stop_pos,3),'.r','MarkerSize', 20,'LineWidth',1.3)
        end
        plot(t, zeros(size(t)),'-k')
        
        if conf == 2
            plot(t(1:conf_pos), evidence_go(1:conf_pos,3),':','color',[0.3 0.2 0.9],'LineWidth',1.3)
            plot(t(1:conf_pos), -evidence_wait(1:conf_pos,3),':','color',[0.32 0.7 0.9],'LineWidth',1.3)
            title('Race Model, CT = RT + \tau')
            txt2 = 'tau CT';
        else 
            title('Race Model, CT = RT')
        end 
        text(t(stop_pos), evidence_go(stop_pos,3)+0.3, 'go','color',[0.3 0.2 0.9], 'FontSize', 14)
        text(t(stop_pos), -evidence_wait(stop_pos,3)-0.3, 'wait','color',[0.32 0.7 0.9], 'FontSize', 14)
        text(0.2,0.85, 'go', 'FontSize', 14)
        text(0.2,-0.85, 'wait', 'FontSize', 14)
        ax = gca;             ax.FontSize = 14;             
        box off 
        ylabel('Evidence (EV)','FontSize', 16);       
        xlabel('time, s','FontSize', 16)
        ylim([-1,max_bound]);           xlim([0,t(conf_pos)+0.2])  
        saveas(figure2, fullfile(fname, [txt, txt2, ' illustration model race.jpg']))
        saveas(figure2, fullfile(fname, [txt, txt2, ' illustration model race.pdf']))
    else 
       
        figure2=figure;  
        if decision(3) == 1 
            plot(t(1:stop_pos), evidence_go(1:stop_pos,3),'-','color',[0.6350 0.0780 0.1840],'LineWidth',1.3)
        else 
            plot(t(1:stop_pos), -evidence_wait(1:stop_pos,3),'-','color',[0.6350 0.0780 0.1840],'LineWidth',1.3)
        end 
        hold on 
        plot(t, upper_boundary_trials(:,3),'--k','LineWidth',1.3)
        plot(t, -low_boundary_trials(:,3),':k','LineWidth',1.3)
        if decision(3) == 1 
            plot(RT_raw_go(3),upper_boundary_trials(stop_pos,3),'.r','MarkerSize', 20,'LineWidth',1.3)
        else 
            plot(RT_raw_wait(3),-low_boundary_trials(stop_pos,3),'.r','MarkerSize', 20,'LineWidth',1.3)
        end
        plot(t, zeros(size(t)), '-k')
        
        if conf == 2 
            if decision(3) == 1 
            plot(t(1:conf_pos), evidence_go(1:conf_pos,3),':','color',[0.6350 0.0780 0.1840],'LineWidth',1.3)
            else 
            plot(t(1:conf_pos), -evidence_wait(1:conf_pos,3),':','color',[0.6350 0.0780 0.1840],'LineWidth',1.3)
            end 
            title('DDM, CT = RT + \tau')
            txt2 = 'tau CT';
        else
            title('DDM, CT = RT')
        end 
        ax = gca;             ax.FontSize = 14;
        text(0.2,0.95, 'go', 'FontSize', 14)
        text(0.2, -0.95, 'wait', 'FontSize', 14)                 
        box off
        ylim([-max_bound,max_bound]);   xlim([0,t(conf_pos)+0.2])
        ylabel('Evidence (EV)','FontSize', 16);        
        xlabel('time, s','FontSize', 16)
        saveas(figure2, fullfile(fname, [txt, txt2, ' illustration model DDM.jpg']))
        saveas(figure2, fullfile(fname, [txt, txt2, ' illustration model DDM.pdf']))
    end   
    end 
end 











