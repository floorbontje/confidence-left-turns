function  [coeff_go, coeff_wait] = function_Conf_Race_Model_x0(trials, p, conf_t)
% Race model 
% Two independent different drift-rates for go/wait decisions 
% For both go/wait decisions the decision boundary is positive
% Dynamic boundary and drift based on the generalized gap (TTA + beta * d)
close all 

if conf_t == 0 
    % CT = RT 
    model = 'model3'; 
    txt = 'Model 3'; 
else  
    % CT > RT 
    model = 'model4'; 
    txt = ['Model 4 (t_e_x_t_r_a =', num2str(conf_t),')']; 
end 

%% Measures data 
data_loc = '.\Data analysis\general data analysis\'; 
load([data_loc,'Confidence_mean_RT.mat']) 

%% Initialise

% Model parameters 
alpha_go = p(1); beta_go=p(2); b0_go=p(3); k_go = p(4);  
mu_ND_go= p(5); sigma_ND_go =p(6); theta_c_go = p(7);

alpha_wait = p(8); beta_wait=p(9); b0_wait=p(10); k_wait = p(11);  
mu_ND_wait = p(12); sigma_ND_wait =p(13); theta_c_wait = p(14);

% Conditions 
tta = [5.5,6.5]; 
dist = [70,90];

% Time 
max_t = 7;                          % Maximum time
delta_t   = 0.001;                  % time steps 
max_inter = max_t/delta_t;          % total interations
t = (0:delta_t:max_t-delta_t)';     % time vector


% Noise
accu_noise_mean = 0; accu_noise_std = 1; 
sen_noise_mean = 0;  sen_noise_std = 0; 
 
% Evidence, decision, RT  
evidence_go = nan(max_inter,trials*length(tta)*length(dist));  
evidence_wait = evidence_go;  
decision_go = zeros (1,trials*length(tta)*length(dist)); 
decision_wait = decision_go; 
RT_raw_go = ones(1,4*trials)*1000;  
RT_raw_wait = RT_raw_go; 



%% Model
% Time-to-arrival & distance gap 
tta_con = [tta,tta];  
dist_con = [dist(1)*ones(1,2),dist(2)*ones(1,2)]; 

TTA = tta_con - t;                                                          
d = dist_con - dist_con ./tta_con.*t;                                        

% generalized gap 
generalized_gap_go   = (TTA + beta_go * d);                                  
generalized_gap_wait = (TTA + beta_wait * d);  

% boundaries
b_go   = b0_go ./ (1 + exp(-k_go * (generalized_gap_go - theta_c_go)));     
b_wait = b0_wait ./ (1 + exp(-k_wait * (generalized_gap_wait - theta_c_wait)));

    
% noise: accumulation noise & sensory noise
accu_noise = accu_noise_std*randn(max_inter, (trials*length(tta)*length(dist))) + accu_noise_mean;  %accumulation noise  
sen_noise  = 1 - (sen_noise_std*randn(max_inter, trials*length(tta)*length(dist)) + sen_noise_mean); %sensory noise 

% Drift-rate 
driftrate_c_go   = alpha_go.*(generalized_gap_go - theta_c_go); 
driftrate_c_wait = alpha_wait.*(generalized_gap_wait - theta_c_wait); 


%% Calculations
driftrate_go    = driftrate_c_go;
driftrate_wait  = driftrate_c_wait; 
tta_con_trials  = tta_con; 
dist_con_trials = dist_con; 
wait_boundary_trials = b_wait; 
go_boundary_trials   = b_go;

for i = 1 : trials - 1 
    driftrate_go(:,end+1:end+4) = driftrate_c_go;
    driftrate_wait(:,end+1:end+4) = driftrate_c_wait;
    
    tta_con_trials (:,end+1:end+4)= tta_con; 
    dist_con_trials (:,end+1:end+4)= dist_con; 
    wait_boundary_trials(:,end+1:end+4) = b_wait; 
    go_boundary_trials(:,end+1:end+4) = b_go; 
end 

for i = 1:max_inter
        dx_go = sen_noise(i,:).*driftrate_go(i,:)*delta_t + accu_noise(i,:)*sqrt(delta_t); 
        dx_wait = -(sen_noise(i,:).*driftrate_wait(i,:)*delta_t + accu_noise(i,:)*sqrt(delta_t)); 
        if i == 1
            evidence_go(i,:) = dx_go;
            evidence_wait(i,:) = dx_wait;
        else 
            evidence_go(i,:) = evidence_go(i-1,:) + dx_go; 
            evidence_wait(i,:) = evidence_wait(i-1,:) + dx_wait; 
        end
        
        RT_raw_go(evidence_go(i,:)> go_boundary_trials(i,:) & decision_go == 0) = t(i);
        decision_go(evidence_go(i,:)> go_boundary_trials(i,:) & decision_go == 0) = 1; 
               
        RT_raw_wait(evidence_wait(i,:) > wait_boundary_trials(i,:) & decision_wait == 0) = t(i);
        decision_wait(evidence_wait(i,:) > wait_boundary_trials(i,:) & decision_wait == 0) = -1;      
end

decision = nan(1,4*trials); 
decision(RT_raw_go < RT_raw_wait) = 1;
decision(RT_raw_go > RT_raw_wait) = -1; 

%% Confindence judgments 
% Time moment of confidence judgment (CT)
CT = nan(size(RT_raw_wait)); 
CT(decision == 1)  = RT_raw_go(decision == 1) + conf_t; 
CT(decision == -1) = RT_raw_wait(decision == -1) + conf_t; 

% Position of CT
int_CT = nan(size(RT_raw_wait));
for i = 1:length(RT_raw_wait)
     int_CT(i) = find(t > CT(i),1); 
end 

% Evidence, decision boundaries: positive decision boundary both decisions
DV_dec = nan(1,4*trials);  DV_alt = DV_dec; 
b_dec = nan(1,4*trials);   b_alt = b_dec ; 

for i = find(decision == 1)
    DV_dec(i) = evidence_go(int_CT(i),i); 
    DV_alt(i) = evidence_wait(int_CT(i),i); 
    b_dec(i)  = go_boundary_trials(int_CT(i),i); 
    b_alt(i)  = wait_boundary_trials(int_CT(i),i); 
end 

for i = find(decision == -1)
    DV_dec(i) = evidence_wait(int_CT(i),i); 
    DV_alt(i) = evidence_go(int_CT(i),i); 
    b_dec(i)  = wait_boundary_trials(int_CT(i),i); 
    b_alt(i)  = go_boundary_trials(int_CT(i),i); 
end 


for i = 1 :length(tta)*length(dist)
    pos_go = find(decision == 1 & tta_con_trials == tta_con(i) & dist_con_trials == dist_con(i)); 
    pos_wait =find(decision == -1 & tta_con_trials == tta_con(i) & dist_con_trials == dist_con(i));    
    
    DV_mean_loss_go(i)   = mean(-DV_alt(pos_go));
    DV_mean_loss_wait(i) = mean(-DV_alt(pos_wait)); 
end 

%% Loss function
% Value losing DV 
Vc= [DV_mean_loss_go';DV_mean_loss_wait']; 


% Linear model
table_output = table(Vc,[Conf_mean_c_go';Conf_mean_c_wait'],...
    categorical([ones(4,1);zeros(4,1)],[1,0],{'Go' 'Wait'}),...
    'VariableNames', {'Vc', 'Conf', 'dec'}); 
model_conf_0 = fitlm(table_output, 'Conf ~ Vc*dec'); 

coeff_go = model_conf_0.Coefficients.Estimate([1,2]); 
coeff_wait = model_conf_0.Coefficients.Estimate([1,2])+model_conf_0.Coefficients.Estimate([3,4]); 











