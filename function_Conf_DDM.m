function  RMSE =function_Conf_DDM(trials, p, conf_t, plot_fig, coeff_go, coeff_wait)
% Dynamic Drift Diffusion Model 
% Dynamic boundary and drift based on the generalized gap (TTA + beta * d)
close all 

if conf_t == 0 
    % CT = RT 
    model = 'Model1';
    txt = 'DDM, CT = RT'; 
else 
    % CT > RT 
    model = 'Model2';
    txt = 'DDM, CT = RT + \tau'; 
end 

%% Measures experiment 
data_loc = '.\Data analysis\general data analysis\'; 
load([data_loc,'Confidence_mean_RT.mat']) 

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

% non decision time 
t_non_decision = sigma_ND*randn(1, trials*length(tta)*length(dist)) +mu_ND;
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


%% Model predictions 

for i = 1: length(decision)
    if decision(i) == 1 
        conf_raw(i) = coeff_go(1)+(Vc(i))*coeff_go(2);
    else 
        conf_raw(i)= coeff_wait(1)+(Vc(i))*coeff_wait(2);
    end 
end 
 

conf = conf_raw; 

conf(conf<1)=1; 
conf(conf>5)=5; 

for i = 1 :length(tta)*length(dist)
    pos_go_c = find(decision == 1 & tta_con_trials == tta_con(i) & dist_con_trials == dist_con(i)); 
    pos_wait_c =find(decision == -1 & tta_con_trials == tta_con(i) & dist_con_trials == dist_con(i));    
    Conf_go(i) = mean(conf(pos_go_c)); 
    Conf_wait(i) = mean(conf(pos_wait_c)); 
    CI_conf_go(i) = 1.96*std(conf(pos_go_c))/length(pos_go_c); 
    CI_conf_wait(i) = 1.96*std(conf(pos_wait_c))/length(pos_wait_c);
end

conf_model = [Conf_go'; Conf_wait']; 
conf_measured = [Conf_mean_c_go';Conf_mean_c_wait']; 

RMSE = sqrt(sum((conf_model - conf_measured).^2)/length(conf_model));
 



 %% Figures 
    if  plot_fig == 1 
        linear_regression_fig=figure; 
        subplot(1,2,1) 
        errorbar(tta, Conf_go(1:2),CI_conf_go(1:2), '.-b', 'LineWidth',1.3,'MarkerSize', 20)
        hold on; box off
        errorbar(tta, Conf_go(3:4),CI_conf_go(3:4), '.-r', 'LineWidth',1.3,'MarkerSize', 20)
        errorbar(tta, Conf_mean_c_go(1:2), CI_conf_mean_c_go(1:2), '.--', 'color', [.73 .73 1], 'LineWidth',1.3,'MarkerSize', 20)
        errorbar(tta, Conf_mean_c_go(3:4), CI_conf_mean_c_go(3:4),'.--', 'color',[1 .73 .73], 'LineWidth',1.3,'MarkerSize', 20) 
        text(5.5,4.7,['c_{0}^{go}: ', num2str(coeff_go(1))], 'FontSize', 12)
        text(5.5,4.4,['c_{1}^{go}: ', num2str(coeff_go(2))], 'FontSize', 12)
        ax = gca; 
        ax.FontSize = 12;
        xlabel('tta [sec]','FontSize', 14);         
        ylabel('Confidence','FontSize', 14) 
        title('Go','FontSize', 16)
        ylim([1,5]); xlim([5.3;6.7]); xticks([5.5, 6.5])
        legend('model 70m', 'model 90m', 'data 70m', 'data 90m', 'Location', 'southeast', 'Fontsize', 12)
        

        subplot(1,2,2) 
        errorbar(tta, Conf_wait(1:2),CI_conf_wait(1:2), '.-b', 'LineWidth',1.3,'MarkerSize', 20)
        hold on; box off
        errorbar(tta, Conf_wait(3:4),CI_conf_wait(3:4), '.-r', 'LineWidth',1.3,'MarkerSize', 20)
        errorbar(tta, Conf_mean_c_wait(1:2),CI_conf_mean_c_wait(1:2), '.--', 'color',[.73 .73 1], 'LineWidth',1.3,'MarkerSize', 20)
        errorbar(tta, Conf_mean_c_wait(3:4),CI_conf_mean_c_wait(3:4), '.--', 'color',[1 .73 .73], 'LineWidth',1.3,'MarkerSize', 20) 
        text(5.5,4.7,['c_{0}^{wait}: ', num2str(coeff_wait(1))], 'FontSize', 12)
        text(5.5,4.4,['c_{1}^{wait}: ', num2str(coeff_wait(2))],'FontSize', 12)
        ax = gca; 
        ax.FontSize = 12;
        xlabel('tta [sec]','FontSize', 14)
        ylabel('Confidence','FontSize', 14) 
        title('Wait','FontSize', 16)
        ylim([1,5]); xlim([5.3;6.7]); xticks([5.5, 6.5])
        dim = [.6 .1 .1 .1];
        sgtitle(txt, 'FontSize', 18) 
       


        figure_Vc = figure; 
        subplot(2,2,1)
        plot(tta, EV_mean_go(1:2),'.-b', 'MarkerSize',15)
        hold on; grid on 
        plot(tta, EV_mean_go(3:4), '.-r', 'MarkerSize',15)
        title ('Go')
        xlim([5.3 6.7]); xticks([5.5, 6.5])
        xlabel('Time-to-arrival (TTA), s')
        ylabel('V_c')

        subplot(2,2,2)
        plot(tta, EV_mean_wait(1:2),'.-b', 'MarkerSize',15)
        hold on; grid on 
        plot(tta, EV_mean_wait(3:4), '.-r', 'MarkerSize',15)
        title ('Wait')
        xlabel('Time-to-arrival (TTA), s')
        xlim([5.3 6.7]); xticks([5.5, 6.5])
        ylabel('V_c')

        subplot(2,2,3)
        plot(tta,Conf_mean_c_go(1:2),'.-b', 'MarkerSize',15)
        hold on; grid on
        plot(tta,Conf_mean_c_go(3:4),'.-r', 'MarkerSize',15)
        legend('70 m', '90m', 'Location', 'Southwest')
        title ('Go')
        xlabel('Time-to-arrival (TTA), s')
        ylim([1,5]); xlim([5.3 6.7]); xticks([5.5, 6.5])
        ylabel('Confidence') 

        subplot(2,2,4)
        plot(tta,Conf_mean_c_wait(1:2),'.-b', 'MarkerSize',15)
        hold on; grid on
        plot(tta,Conf_mean_c_wait(3:4),'.-r', 'MarkerSize',15)
        title ('Wait decisions')
        xlabel('Time-to-arrival (TTA), s')
        ylabel('Confidence')
        ylim([1,5]); xlim([5.3 6.7]); xticks([5.5, 6.5])
        sgtitle(txt)




   %% Save figures 
        fname = ['.\Modelling\Confidence\', model];
        saveas(linear_regression_fig, fullfile(fname,['linear regression model ' , model, '.jpg']))
        saveas(figure_Vc, fullfile(fname,['Raw evidence and conf ', model, '.jpg']))
        saveas(linear_regression_fig, fullfile(fname,['linear regression model ' , model, '.pdf']))
        saveas(figure_Vc, fullfile(fname,['Raw evidence and conf ', model, '.pdf']))
        
        
        RT = RT_raw + t_non_decision ;
        

        is_go_decision=categorical(decision,[1 -1],{'True' 'False'});

        data_output= table(RT',is_go_decision',tta_con_trials', dist_con_trials',Vc', conf',...
            'VariableNames',{'RT','is_go_decision','tta_condition','d_condition','Vc', 'confidence'});  

        save_name = [fname, '\data_output_',model,'.csv']; 
        writetable(data_output, save_name)
    end 





