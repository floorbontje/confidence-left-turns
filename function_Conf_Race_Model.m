function  [RMSE]  = function_Conf_Race_Model(trials, p, conf_t, plot_fig, coeff_go, coeff_wait)
    % Race model
    % Two independent different drift-rates for go/wait decisions 
    % For both go/wait decisions the decision boundary is positive
    % Dynamic boundary and drift based on the generalized gap (TTA + beta * d)

    % Race model 
    % Two independent different drift-rates for go/wait decisions 
    % For both go/wait decisions the decision boundary is positive
    % Dynamic boundary and drift based on the generalized gap (TTA + beta * d)
    close all 

    if conf_t == 0 
        % CT = RT 
        model = 'Model3'; 
        txt = 'Race Model, CT = RT'; 
    else  
        % CT > RT 
        model = 'Model4'; 
        txt = 'Race Model, CT = RT + \tau'; 
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

    % Non-decision time  
    t_noise_go      = sigma_ND_go*randn(1, trials*length(tta)*length(dist)); 
    t_noise_wait    = sigma_ND_wait*randn(1, trials*length(tta)*length(dist)); 

    if sigma_ND_wait == sigma_ND_go 
        t_noise_wait = t_noise_go; 
    end 

    t_non_decision_go   = t_noise_go+mu_ND_go;
    t_non_decision_go(t_non_decision_go < 0 ) = 0; 
    t_non_decision_wait = t_noise_wait+mu_ND_wait;
    t_non_decision_wait(t_non_decision_wait < 0 ) = 0; 



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


%% Model predictions 

    Vc= -DV_alt; 

    for i = 1: length(decision)
        if decision(i) == 1 
            conf(i) = coeff_go(1)+(Vc(i))*coeff_go(2);
        else 
            conf(i) = coeff_wait(1)+(Vc(i))*coeff_wait(2);
        end 
    end 
    
    % Confidence boundaries (min/max values)
    conf(conf<1)=1; 
    conf(conf>5)=5; 

    % Mean confidence jdugments 
    for i = 1 :length(tta)*length(dist)
        pos_go_c = find(decision == 1 & tta_con_trials == tta_con(i) & dist_con_trials == dist_con(i)); 
        pos_wait_c =find(decision == -1 & tta_con_trials == tta_con(i) & dist_con_trials == dist_con(i));    
        Conf_go(i) = mean(conf(pos_go_c)); 
        Conf_wait(i) = mean(conf(pos_wait_c)); 
        CI_conf_go(i) = 1.96*std(conf(pos_go_c))/sqrt(length(pos_go_c)); 
        CI_conf_wait(i) = 1.96*std(conf(pos_wait_c))/sqrt(length(pos_wait_c)); 
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
        text(5.5,4.7,['c_{0,Go}: ', num2str(coeff_go(1))], 'FontSize', 12)
        text(5.5,4.4,['c_{go}   : ', num2str(coeff_go(2))], 'FontSize', 12)
        ax = gca;          ax.FontSize = 12;
        xlabel('tta [sec]','FontSize', 14);         
        ylabel('Confidence','FontSize', 14) 
        title('Go','FontSize', 16)
        ylim([1,5]); xlim([5.3;6.7])
        legend('model 70m', 'model 90m', 'data 70m', 'data 90m', 'Location', 'southeast', 'Fontsize', 12)
        

        subplot(1,2,2) 
        errorbar(tta, Conf_wait(1:2),CI_conf_wait(1:2), '.-b', 'LineWidth',1.3,'MarkerSize', 20)
        hold on; box off
        errorbar(tta, Conf_wait(3:4),CI_conf_wait(3:4), '.-r', 'LineWidth',1.3,'MarkerSize', 20)
        errorbar(tta, Conf_mean_c_wait(1:2),CI_conf_mean_c_wait(1:2), '.--', 'color',[.73 .73 1], 'LineWidth',1.3,'MarkerSize', 20)
        errorbar(tta, Conf_mean_c_wait(3:4),CI_conf_mean_c_wait(3:4), '.--', 'color',[1 .73 .73], 'LineWidth',1.3,'MarkerSize', 20) 
        ax = gca;         ax.FontSize = 12;
        text(5.5,4.7,['c_{0,wait}: ', num2str(coeff_wait(1))], 'FontSize', 12)
        text(5.5,4.4,['c_{wait}   : ', num2str(coeff_wait(2))],'FontSize', 12)
        xlabel('tta [sec]','FontSize', 14)
        ylabel('Confidence','FontSize', 14) 
        title('Wait','FontSize', 16)
        ylim([1,5]); xlim([5.3;6.7])
        dim = [.6 .1 .1 .1];
        sgtitle(txt, 'FontSize', 18) 
       


        figure_Vc = figure; 
        subplot(2,2,1)
        plot(tta, DV_mean_loss_go(1:2),'.-b', 'MarkerSize',15)
        hold on; grid on 
        plot(tta, DV_mean_loss_go(3:4), '.-r', 'MarkerSize',15)
        title ('Go')
        xlim([5.3 6.7])
        xlabel('Time-to-arrival (TTA), s')
        ylabel('V_c')

        subplot(2,2,2)
        plot(tta, DV_mean_loss_wait(1:2),'.-b', 'MarkerSize',15)
        hold on; grid on 
        plot(tta, DV_mean_loss_wait(3:4), '.-r', 'MarkerSize',15)
        title ('Wait')
        xlabel('Time-to-arrival (TTA), s')
        xlim([5.3 6.7])
        ylabel('V_c')

        subplot(2,2,3)
        plot(tta,Conf_mean_c_go(1:2),'.-b', 'MarkerSize',15)
        hold on; grid on
        plot(tta,Conf_mean_c_go(3:4),'.-r', 'MarkerSize',15)
        legend('70 m', '90m', 'Location', 'Southwest')
        title ('Go')
        xlabel('Time-to-arrival (TTA), s')
        ylim([1,5]); xlim([5.3 6.7])
        ylabel('Confidence') 

        subplot(2,2,4)
        plot(tta,Conf_mean_c_wait(1:2),'.-b', 'MarkerSize',15)
        hold on; grid on
        plot(tta,Conf_mean_c_wait(3:4),'.-r', 'MarkerSize',15)
        title ('Wait decisions')
        xlabel('Time-to-arrival (TTA), s')
        ylabel('Confidence')
        ylim([1,5]); xlim([5.3 6.7])
        sgtitle(txt)




        %% Save figures 
        fname = ['.\Modelling\Confidence\', model];
        saveas(linear_regression_fig, fullfile(fname,['linear regression model ' , model, '.jpg']))
        saveas(figure_Vc, fullfile(fname,['Raw evidence and conf ', model, '.jpg']))
        saveas(linear_regression_fig, fullfile(fname,['linear regression model ' , model, '.pdf']))
        saveas(figure_Vc, fullfile(fname,['Raw evidence and conf ', model, '.pdf']))
        
        RT = size(RT_raw_go); 
        RT(decision == -1) = RT_raw_wait(decision == -1) +t_non_decision_wait(decision == -1); 
        RT(decision == 1) = RT_raw_go (decision == 1) + t_non_decision_go(decision ==1); 

        is_go_decision=categorical(decision,[1 -1],{'True' 'False'});

        data_output= table(RT',is_go_decision',tta_con_trials', dist_con_trials',Vc', conf',...
            'VariableNames',{'RT','is_go_decision','tta_condition','d_condition','Vc', 'confidence'});  

        save_name = [fname, '\data_output_',model,'.csv']; 
        writetable(data_output, save_name)
    end 
 






