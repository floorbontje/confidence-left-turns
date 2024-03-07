function  [coeff_go, coeff_wait] =function_Conf_DDM_x0(trials, p, conf_t)

% Dynamic Drift Diffusion Model 
% Dynamic boundary and drift based on the generalized gap (TTA + beta * d)
close all 

if conf_t == 0 
    % CT = RT 
    txt = 'Model 1';
    model = 'model1'; 
else 
    % CT > RT 
    model = 'model2';
    txt = ['Model 2 (t_e_x_t_r_a =', num2str(conf_t),')']; 
end 

%% Measures experiment 
load('.\Data analysis\general data analysis\Confidence_mean_RT.mat')
data_loc = '.\Data anlysis\general data analysis\'; 
% load([data_loc,'Confidence_mean_RT.mat']) 

%% Initialise
% Model parameters 
alpha =p(1); beta=p(2); b0= p(3); k = p(4);  
mu_ND= p(5); sigma_ND =p(6); theta_c = p(7);


% Conditions
tta = [5.5,6.5]; 
dist = [70,90];

% Time 
max_t = 7;                            % Maximum time
delta_t   = 0.001;                    % time steps 
max_inter = max_t/delta_t;            % total interations
t = (0:delta_t:max_t-delta_t)';       % time vector

% Noise 
accu_noise_mean = 0;  accu_noise_std = 1; 
sen_noise_mean = 0;   sen_noise_std = 0;

% Evidence, decision, RT 
evidence = nan(max_inter,trials*length(tta)*length(dist)); 
decision = zeros (1,trials*length(tta)*length(dist));     
RT_raw = ones(1,4*trials)*1000;



%% Model 
% Time-to-arrival & distance gap 
tta_con = [tta,tta];  
dist_con = [dist(1)*ones(1,2),dist(2)*ones(1,2)]; 

TTA = tta_con - t;                                                          
d = dist_con - dist_con ./tta_con.*t;                                       

% generalized gap
generalized_gap = (TTA + beta * d);   

% boundary 
b = b0 ./ (1 + exp(-k * (generalized_gap - theta_c)));                      
upperboundary = b; 
lowerboundary = -b; 
 
    
% Noise: acculumatio noise & sensory noise 
accu_noise = accu_noise_std*randn(max_inter, (trials*length(tta)*length(dist))) + accu_noise_mean; 
sen_noise  = 1 - (sen_noise_std*randn(max_inter, trials*length(tta)*length(dist)) + sen_noise_mean); 


% Drift-rate 
driftrate_c = alpha.*(generalized_gap - theta_c); 


%% Calculations 
driftrate       = driftrate_c; 
tta_con_trials  = tta_con; 
dist_con_trials = dist_con; 
low_boundary_trials   = lowerboundary; 
upper_boundary_trials = upperboundary;

for i = 1 : trials - 1 
    driftrate(:,end+1:end+4) = driftrate_c;
    tta_con_trials (:,end+1:end+4)= tta_con; 
    dist_con_trials (:,end+1:end+4)= dist_con; 
    low_boundary_trials(:,end+1:end+4) = lowerboundary; 
    upper_boundary_trials(:,end+1:end+4) = upperboundary;
end 


for i = 1:max_inter
        dx = sen_noise(i,:).*driftrate(i,:)*delta_t + accu_noise(i,:)*sqrt(delta_t); 
        if i == 1
            evidence(i,:) = dx; 
        else 
            evidence(i,:) = evidence(i-1,:) + dx; 
        end 
        
        RT_raw(evidence(i,:)<low_boundary_trials(i,:) & decision == 0) = t(i);
        RT_raw(evidence(i,:)>upper_boundary_trials(i,:) & decision == 0) = t(i);
        decision(evidence(i,:)<low_boundary_trials(i,:) & decision == 0) = -1; 
        decision(evidence(i,:)>upper_boundary_trials(i,:) & decision == 0) = 1; 
end

% Evidence for confidence judgment 
for i = 1: length(RT_raw)
    if abs(decision(i))>0
        pos(i) = find(t > (RT_raw(i)+conf_t),1);
        EV(i) = evidence(pos(i),i);
    else 
        EV(i) = 0; 
    end 
end 

% Evidence
for i = 1 :length(tta)*length(dist)
    pos_go = find(decision == 1 & tta_con_trials == tta_con(i) & dist_con_trials == dist_con(i)); 
    pos_wait =find(decision == -1 & tta_con_trials == tta_con(i) & dist_con_trials == dist_con(i));    
    
    EV_mean_go(i) = mean(EV(pos_go));
    EV_mean_wait(i) = mean(-EV(pos_wait)); 
    Vc(pos_go) = EV(pos_go); 
    Vc(pos_wait) = -EV(pos_wait); 
    
    CI_EV_go(i) = 1.96 * std(EV(pos_go))/sqrt(length(pos_go));
    CI_EV_wait(i) = 1.96 * std(-EV(pos_wait))/sqrt(length(pos_wait));
end

%% Loss function
Vc_mean= [EV_mean_go';EV_mean_wait']; 

% Linear model
table_output = table(Vc_mean,[Conf_mean_c_go';Conf_mean_c_wait'],...
    categorical([ones(4,1);zeros(4,1)],[1,0],{'Go' 'Wait'}),...
    'VariableNames', {'Vc', 'Conf', 'dec'}); 
model_conf_0 = fitlm(table_output, 'Conf ~ Vc*dec'); 

coeff_go = model_conf_0.Coefficients.Estimate([1,2]); 
coeff_wait = model_conf_0.Coefficients.Estimate([1,2])+model_conf_0.Coefficients.Estimate([3,4]); 

Rsq  = model_conf_0.Rsquared.Ordinary; 
RMSE = model_conf_0.RMSE; 

Loss_conf = RMSE; 

Conf_go = coeff_go(1)+EV_mean_go*coeff_go(2); 
Conf_wait = coeff_wait(1)+EV_mean_wait*coeff_wait(2); 




