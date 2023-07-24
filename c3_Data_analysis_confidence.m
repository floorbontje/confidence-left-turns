% Data analysis - Confidence Durning Left Turns  
% By Floor Bontje  

close all
clear 
clc 

fname_lme  = '.\Data analysis\figures\analysis'; 
fname_dist = '.\Data analysis\figures\data distribution';
fname_output = '.\Data analysis\general data analysis';

%%  Load data 
load Conf_all
load Log_all 

 
[Go, Wait, change_of_mind,no_button,time, start_ego_inter, end_ego_inter,change_of_mind_kind] = All_inter(Log_all, Conf_all, []);
[~, ~, Go_c, Wait_c,Go_c_pr, Wait_c_pr]= Conf_global_op(Log_all, Conf_all, Go, Wait, change_of_mind, no_button); 

ID      = unique(Conf_all(:,1))'; 
tta     = unique(Conf_all(:,5))'; 
dist    = unique(Conf_all(:,6))';
nbr_int = length(Conf_all);  

%% No button presses and changes of mind 
nbr_change_of_mind = length(change_of_mind); 
nbr_no_button = length(no_button); 

nbr_go_nb = 0; nbr_wait_nb = 0;  
for i = 1:nbr_no_button
    if ~isempty(find(Go == no_button(i), 1))
    nbr_go_nb = nbr_go_nb + 1;  
    end 
    if ~isempty (find(Wait == no_button(i),1))
    nbr_wait_nb = nbr_wait_nb + 1;
    end 
end 

Variable = {'nbr_change_of_mind'; 'nbr_no_button'; 'GO: nbr_no_button'; 'Wait: nbr_no_button'; ...
    'prs_changes_of_mind_and_no_button'; 'prs_changes_of_mind'; 'prs_no_button'; ...
    'Go: prs_no_button';'Wait: prs_no_button'; 'Go_to_wait'; 'Wait_to_go';  ...
    'Go_to_wait: no indication'; 'Go_to_wait: indication'; 'Go_to_wait: False indication';...
    'Wait_to_go: no indication'; 'Wait_to_go: indication'; 'Wait_to_go: False indication'};

nbr_count  = nan(size(Variable)); 
percentage = nan(size(Variable)); 

nbr_count(1:4)      = [nbr_change_of_mind; nbr_no_button; nbr_go_nb; nbr_wait_nb];
percentage(5:7)     = [(nbr_change_of_mind + nbr_no_button); nbr_change_of_mind;nbr_no_button]/nbr_int*100;
percentage(8:9)     = [nbr_go_nb; nbr_wait_nb]/nbr_no_button*100; 
percentage(10:17)   = [length(find(change_of_mind_kind==1 | change_of_mind_kind==2));length(find(change_of_mind_kind==-1 | change_of_mind_kind==-2));...
    length(find(change_of_mind_kind==1)); length(find(change_of_mind_kind==2)); length(find(change_of_mind_kind==3)); ...
    length(find(change_of_mind_kind==-1));length(find(change_of_mind_kind==-2));length(find(change_of_mind_kind==-3))]/nbr_change_of_mind*100; 

table_excluded_inter= table(Variable, nbr_count, percentage, ...
    'VariableNames',{'Variable','nbr','percentage'});  
writetable(table_excluded_inter, [fname_output '\excluded_intersections.csv'])


%% Input for modelling 
prs_decision = ones(length(Conf_all(:,1)),1);
prs_decision(Wait)=0; 
decision = ones(length(Conf_all(:,1)),1);
decision(Wait)= 2; 
decision = categorical(decision,[1 2],{'Go' 'Wait'});

time_first_throttle_use = NaN(length(start_ego_inter),1); 
for i = 1:length(start_ego_inter)
    pos_first_throttle=start_ego_inter(i)+find(Log_all(start_ego_inter(i):end_ego_inter(i),29)>0,1);
    time_first_throttle_use(i)=Log_all(pos_first_throttle,7)-Log_all(start_ego_inter(i),7);
end 

table_all= table(Conf_all(:,1),Conf_all(:,4),Conf_all(:,5),Conf_all(:,6),...
    decision,time.delta_1, time.delta_2, time.delta_3, prs_decision, time_first_throttle_use,...
    'VariableNames',{'ID','conf','tta','dist','decision','reaction_time','turn_time','judg_time','prs_decision','first_thr'});  
table_all([change_of_mind,no_button],:) = []; 


% Mean confidence and reaction times (Go/Wait) 
Mean_conf_go=mean(table_all.conf(table_all.prs_decision == 1));
Mean_conf_wait=mean(table_all.conf(table_all.prs_decision == 0));
std_conf_go=std(table_all.conf(table_all.prs_decision == 1));
std_conf_wait=std(table_all.conf(table_all.prs_decision == 0));
Mean_RT_go=mean(table_all.reaction_time(table_all.prs_decision == 1));
Mean_RT_wait=mean(table_all.reaction_time(table_all.prs_decision == 0)); 
std_RT_go=std(table_all.reaction_time(table_all.prs_decision == 1)); 
std_RT_wait=std(table_all.reaction_time(table_all.prs_decision == 0)); 


% Positions in table, mean values & CIs 
con = ['c1'; 'c2'; 'c3'; 'c4'];
tta_con = [tta, tta]; 
dist_con = [dist(1),dist(1),dist(2),dist(2)]; 



for i = 1: 4 
    pos.(con(i,:))= find(table_all.tta== tta_con(i) &  table_all.dist==dist_con(i)); 
    pos_go.(con(i,:)) = find(table_all.tta== tta_con(i) &  table_all.dist==dist_con(i)& table_all.decision=='Go');
    pos_wait.(con(i,:)) = find(table_all.tta== tta_con(i) &  table_all.dist==dist_con(i) & table_all.decision=='Wait');  
    
    Conf_mean_c_go_exp(i) = mean(Conf_all(Go_c.(con(i,:)),4));
    Conf_mean_c_wait_exp(i) = mean(Conf_all(Wait_c.(con(i,:)),4));
    
    CI_mean_conf_go_exp(i) = 1.96*std(Conf_all(Go_c.(con(i,:)),4))/sqrt(length(Go_c.(con(i,:)))); 
    CI_mean_conf_wait_exp(i) = 1.96*std(Conf_all(Wait_c.(con(i,:)),4))/sqrt(length(Wait_c.(con(i,:))));
    
    Conf_mean_rt_go_exp(i) = mean(table_all.reaction_time(pos_go.(con(i,:))));
    Conf_mean_rt_wait_exp(i) = mean(table_all.reaction_time(pos_wait.(con(i,:))));
    
    CI_mean_RT_mean_c_go(i) = 1.96*std(table_all.reaction_time(pos_go.(con(i,:))))/sqrt(length(pos_go.(con(i,:)))); 
    CI_mean_RT_mean_c_wait(i) = 1.96*std(table_all.reaction_time(pos_wait.(con(i,:))))/sqrt(length(pos_wait.(con(i,:)))); 
    
    for indiv = 1: length(ID)
        RT_mean_go_C(indiv,i)=mean(table_all.reaction_time(pos_go.(con(i,:))(table_all.ID(pos_go.(con(i,:)))==ID(indiv))));
        RT_mean_wait_C(indiv,i)=mean(table_all.reaction_time(pos_wait.(con(i,:))(table_all.ID(pos_wait.(con(i,:)))==ID(indiv))));
    end 
    
    SEM_Go_thr(i)=1.96*std(time_first_throttle_use(pos_go.(con(i,:))))/sqrt(length(pos_go.(con(i,:)))); 
    Mean_Go_thr(i) = mean(time_first_throttle_use(pos_go.(con(i,:)))); 
end 

pos_go_all = find(table_all.decision == 'Go');
pos_wait_all = find(table_all.decision == 'Wait'); 

% Throttle use w.r.t RT 
delta_time_thr_RT = table_all.first_thr - table_all.reaction_time; 
prs_thr_before_RT = length(find(delta_time_thr_RT <0))/length(delta_time_thr_RT)*100; 
prs_thr_before_go_RT = length(find(delta_time_thr_RT(pos_go_all) <0))/length(pos_go_all)*100; 
prs_thr_before_wait_RT = length(find(delta_time_thr_RT(pos_wait_all) <0))/length(pos_wait_all)*100; 


%% Linear mixed effects models 
close all
clc
% Decision behaviour 
lme_all_go_prs=fitlme(table_all, 'prs_decision ~ tta + dist + (1|ID)','CheckHessian', true);
[random_effect_inter_decision,~] = ...
    function_lme_results ('Probability of go decision', lme_all_go_prs, pos, [], [Go_c_pr.c1,Go_c_pr.c2,Go_c_pr.c3,Go_c_pr.c4], zeros(1,4), ...
    [],[]); 

% Response time
lme_all_rt = fitlme(table_all, 'reaction_time ~ tta *decision + dist *decision + (decision|ID)','CheckHessian', true);

[random_effect_inter_RT,random_effect_wait_RT] = ...
    function_lme_results ('Response time (RT), s', lme_all_rt, pos_go, pos_wait, Conf_mean_rt_go_exp, CI_mean_RT_mean_c_go, ...
    Conf_mean_rt_wait_exp,  CI_mean_RT_mean_c_wait);

% Confidence  
lme_all_conf = fitlme(table_all,'conf ~ reaction_time*decision + tta*decision +dist*decision+(decision|ID)','CheckHessian', true);
lme_all_conf_thr = fitlme(table_all, 'conf ~ first_thr*decision + tta *decision + dist *decision + (decision|ID)','CheckHessian', true); 

[random_effect_inter_conf,random_effect_wait_conf] = ...
    function_lme_results ('Conf', lme_all_conf, pos_go, pos_wait, Conf_mean_c_go_exp, CI_mean_conf_go_exp, ...
    Conf_mean_c_wait_exp,  CI_mean_conf_wait_exp); 

[random_effect_inter_conf_thr,random_effect_wait_conf_thr] = ...
    function_lme_results ('Conf (thr)', lme_all_conf_thr, pos_go, pos_wait, Conf_mean_c_go_exp, CI_mean_conf_go_exp, ...
    Conf_mean_c_wait_exp,  CI_mean_conf_wait_exp); 


 
%% Post-Hoc analysis of linear mixed effects models 

% Response time 
H_rt_go   = [eye(3,3), zeros(3,3)]; 
H_rt_wait = [eye(3,3), eye(3,3)]; 

for i = 1 : 3 
    % decision (i=1), tta (i=2), dist (i=3)
    [pVal_RT_go(i), F_RT_go(i), DF1_1dec, DF2_1dec] = coefTest(lme_all_rt, H_rt_go(i,:));
    [pVal_RT_wait(i), F_RT_wait(i), ~, ~] = coefTest(lme_all_rt, H_rt_wait(i,:));
    [pVal_RT_go_wait(i), F_RT_go_wait(i), DF1_2dec, DF2_2dec] = coefTest(lme_all_rt, [H_rt_go(i,:);H_rt_wait(i,:)]);
end 

F_crit_1dec_RT = finv(0.95, DF1_1dec, DF2_1dec); % Critical F-value for alpha = 0.05
F_crit_2dec_RT = finv(0.95, DF1_2dec, DF2_2dec);

chategory = {'go', 'wait', 'go : wait'}'; 
int_RT = table(chategory, [F_RT_go(1); F_RT_wait(1); F_RT_go_wait(1)],[pVal_RT_go(1);pVal_RT_wait(1);pVal_RT_go_wait(1)],...
    'VariableNames', {'Decision', 'F', 'pValue'});
TTA_RT  = table(chategory, [F_RT_go(2); F_RT_wait(2); F_RT_go_wait(2)],[pVal_RT_go(2);pVal_RT_wait(2);pVal_RT_go_wait(2)],...
    'VariableNames', {'Decision', 'F', 'pValue'});
dis_RT  = table(chategory, [F_RT_go(3); F_RT_wait(3); F_RT_go_wait(3)],[pVal_RT_go(3);pVal_RT_wait(3);pVal_RT_go_wait(3)],...
    'VariableNames', {'Decision', 'F', 'pValue'});

writetable(int_RT, fullfile([fname_output,'\lme models'], 'Ftest RT decision inter.csv'));
writetable(TTA_RT, fullfile([fname_output,'\lme models'], 'Ftest RT TTA.csv'));
writetable(dis_RT, fullfile([fname_output,'\lme models'], 'Ftest RT distance.csv'));

% Confidencce 
H_conf_go = [eye(4,4), zeros(4,4)]; 
H_conf_wait = [eye(4,4), eye(4,4)]; 

for i = 1 : 4 
    % decision (i=1), tta(i=2), dist(i = 3), time (i=4) 
    [pVal_rt_go(i), F_rt_go(i), DF1_1dec, DF2_1dec] = coefTest(lme_all_conf, H_conf_go(i,:));
    [pVal_rt_wait(i), F_rt_wait(i), ~, ~] = coefTest(lme_all_conf, H_conf_wait(i,:));
    [pVal_rt_go_wait(i), F_rt_go_wait(i), DF1_2dec, DF2_2dec] = coefTest(lme_all_conf, [H_conf_go(i,:);H_conf_wait(i,:)]);

    [pVal_thr_go(i), F_thr_go(i), ~, ~] = coefTest(lme_all_conf_thr, H_conf_go(i,:));
    [pVal_thr_wait(i), F_thr_wait(i), ~, ~] = coefTest(lme_all_conf_thr, H_conf_wait(i,:));
    [pVal_thr_go_wait(i), F_thr_go_wait(i), ~, ~] = coefTest(lme_all_conf_thr, [H_conf_go(i,:);H_conf_wait(i,:)]);
end 

F_crit_1dec_conf = finv(0.95, DF1_1dec, DF2_1dec); % Critical F-value for alpha = 0.05
F_crit_2dec_conf = finv(0.95, DF1_2dec, DF2_2dec);


int_conf_rt = table(chategory, [F_rt_go(1); F_rt_wait(1); F_rt_go_wait(1)],[pVal_rt_go(1);pVal_rt_wait(1);pVal_rt_go_wait(1)],...
    'VariableNames', {'Decision', 'F', 'pValue'});
TTA_conf_rt  = table(chategory, [F_rt_go(2); F_rt_wait(2); F_rt_go_wait(2)],[pVal_rt_go(2);pVal_rt_wait(2);pVal_rt_go_wait(2)],...
    'VariableNames', {'Decision', 'F', 'pValue'});
dis_conf_rt  = table(chategory, [F_rt_go(3); F_rt_wait(3); F_rt_go_wait(3)],[pVal_rt_go(3);pVal_rt_wait(3);pVal_rt_go_wait(3)],...
    'VariableNames', {'Decision', 'F', 'pValue'});
time_conf_rt = table(chategory, [F_rt_go(4); F_rt_wait(4); F_rt_go_wait(4)],[pVal_rt_go(4);pVal_rt_wait(4);pVal_rt_go_wait(4)],...
    'VariableNames', {'Decision', 'F', 'pValue'});

int_conf_thr = table(chategory, [F_thr_go(1); F_thr_wait(1); F_thr_go_wait(1)],[pVal_thr_go(1);pVal_thr_wait(1);pVal_thr_go_wait(1)],...
    'VariableNames', {'Decision', 'F', 'pValue'});
TTA_conf_thr  = table(chategory, [F_thr_go(2); F_thr_wait(2); F_thr_go_wait(2)],[pVal_thr_go(2);pVal_thr_wait(2);pVal_thr_go_wait(2)],...
    'VariableNames', {'Decision', 'F', 'pValue'});
dis_conf_thr  = table(chategory, [F_thr_go(3); F_thr_wait(3); F_thr_go_wait(3)],[pVal_thr_go(3);pVal_thr_wait(3);pVal_thr_go_wait(3)],...
    'VariableNames', {'Decision', 'F', 'pValue'});
time_conf_thr = table(chategory, [F_thr_go(4); F_thr_wait(4); F_thr_go_wait(4)],[pVal_thr_go(4);pVal_thr_wait(4);pVal_thr_go_wait(4)],...
    'VariableNames', {'Decision', 'F', 'pValue'});

writetable(int_conf_rt, fullfile([fname_output,'\lme models'], 'Ftest Conf decision inter (RT).csv'));
writetable(TTA_conf_rt, fullfile([fname_output,'\lme models'], 'Ftest Conf TTA (RT).csv'));
writetable(dis_conf_rt, fullfile([fname_output,'\lme models'], 'Ftest Conf distance (RT).csv'));
writetable(time_conf_rt, fullfile([fname_output,'\lme models'], 'Ftest Conf time (RT).csv'));

writetable(int_conf_thr, fullfile([fname_output,'\lme models'], 'Ftest Conf decision inter (Thr).csv'));
writetable(TTA_conf_thr, fullfile([fname_output,'\lme models'], 'Ftest Conf TTA (Thr).csv'));
writetable(dis_conf_thr, fullfile([fname_output,'\lme models'], 'Ftest Conf distance (Thr).csv'));
writetable(time_conf_thr, fullfile([fname_output,'\lme models'], 'Ftest Conf time (Thr).csv'));

%% Correlations 
% Response time - Confidence
[c_rt,p_corr_rt]=corrcoef(table_all.conf,table_all.reaction_time); 
[c_go_rt,p_corr_go_rt]=corrcoef(table_all.conf(pos_go_all),table_all.reaction_time(pos_go_all));
[c_wait_rt,p_corr_wait_rt]=corrcoef(table_all.conf(pos_wait_all),table_all.reaction_time(pos_wait_all));

for i = 1: 4 
    for ii = 1 : 3 
        if ii == 1
            input = 'All';
            pos_cor = pos.(con(i,:));
        elseif ii == 2
            input = 'Go';
             pos_cor = pos_go.(con(i,:));
        elseif ii == 3
            input = 'Wait';
            pos_cor = pos_wait.(con(i,:));
        end         
    [Corr_conf_rt_c,P_corr_rt_c]=corrcoef(table_all.conf(pos_cor),table_all.reaction_time(pos_cor),'Rows','pairwise');
    Corr_conf_rt.(input)(i) = Corr_conf_rt_c(2,1);      
    P_corr_rt.(input)(i) = P_corr_rt_c(2,1); 
    end 
end


% confidence - intial throttle operation moment  
[corr_conf_thr,p_corr_conf_thr]=corrcoef(table_all.conf,table_all.first_thr);
[corr_conf_thr_go,p_corr_conf_thr_go]=corrcoef(table_all.conf(pos_go_all),table_all.first_thr(pos_go_all));
[corr_conf_thr_wait,p_corr_conf_thr_wait]=corrcoef(table_all.conf(pos_wait_all),table_all.first_thr(pos_wait_all));

for i = 1: 4 
    [corr_conf_thr_go_int,p_corr_conf_thr_go_int]=corrcoef(table_all.conf(pos_go.(con(i,:))),table_all.first_thr(pos_go.(con(i,:))));
    [corr_conf_thr_wait_int,p_corr_conf_thr_wait_int]=corrcoef(table_all.conf(pos_wait.(con(i,:))),table_all.first_thr(pos_wait.(con(i,:))));
    corr_conf_thr_go_c(i) = corr_conf_thr_go_int(2,1); 
    p_corr_conf_thr_go_c(i) = p_corr_conf_thr_go_int(2,1); 
    corr_conf_thr_wait_c(i) = corr_conf_thr_wait_int(2,1); 
    p_corr_conf_thr_wait_c(i) = p_corr_conf_thr_wait_int(2,1);    
end 

% response time - intial throttle operation moment  
[corr_RT_thr_,p_corr_RT_thr]=corrcoef(table_all.reaction_time,table_all.first_thr);
[corr_RT_thr_go,p_corr_RT_thr_go]=corrcoef(table_all.reaction_time(pos_go_all),table_all.first_thr(pos_go_all));
[corr_RT_thr_wait,p_corr_RT_thr_wait]=corrcoef(table_all.reaction_time(pos_wait_all),table_all.first_thr(pos_wait_all));

for i = 1: 4 
    [corr_rt_thr_go_int,p_corr_rt_thr_go_int]=corrcoef(table_all.reaction_time(pos_go.(con(i,:))),table_all.first_thr(pos_go.(con(i,:))));
    [corr_rt_thr_wait_int,p_corr_rt_thr_wait_int]=corrcoef(table_all.reaction_time(pos_wait.(con(i,:))),table_all.first_thr(pos_wait.(con(i,:))));
    corr_rt_thr_go_c(i) = corr_rt_thr_go_int(2,1); 
    p_corr_rt_thr_go_c(i) = p_corr_rt_thr_go_int(2,1); 
    corr_rt_thr_wait_c(i) = corr_rt_thr_wait_int(2,1); 
    p_corr_rt_thr_wait_c(i) = p_corr_rt_thr_wait_int(2,1);    
end 


%% Action dynamics 
clc; close all
con = ['c1'; 'c2';'c3'; 'c4'];
condition_4 = ['5.5 sec & 70m'; '6.5 sec & 70m'; '5.5 sec & 90m'; '6.5 sec & 90m']; 

if length(start_ego_inter)-length(table_all.reaction_time)>0
    start_ego_inter([change_of_mind,no_button]) = []; 
    end_ego_inter([change_of_mind,no_button]) = []; 
end 

% analysis with respect to throttle input
time_moment = 1; 
delta_turn_thr = table_all.turn_time+table_all.reaction_time-table_all.first_thr; 

% Time interval (end moment of time 75% of the turns finished)
[f_go,t_go]=ecdf(delta_turn_thr(pos_go_all)); 
[f_wait,t_wait]=ecdf(delta_turn_thr(pos_wait_all));

t_lim_go= [0,t_go(find(f_go>0.75,1))];
t_lim_wait=[0,t_wait(find(f_wait>0.75,1))]; 


% Velocity (all condition)
dynamic = 1; 
[table_time_decision_vel]=...
    analysis_action_dynamics(table_all,Log_all,start_ego_inter,end_ego_inter,...
   dynamic,t_lim_go,t_lim_wait, [],time_moment); 
pause 

% Velocity (each condition individually)
% for i = 1:4 
% [velocity_table_time_decision_thr_c.(con(i,:))]=...
%  analysis_action_dynamics(table_all(pos.(con(i,:)),:),Log_all,start_ego_inter(pos.(con(i,:))),end_ego_inter(pos.(con(i,:))),...
%    dynamic,t_lim_go,t_lim_wait, condition_4(i,:),time_moment); 
% pause
% end 

 
% Distance to centre of intersection (all conditions)
dynamic = 2; 
condition = []; 
[table_time_decision_distance]=...
    analysis_action_dynamics(table_all,Log_all,start_ego_inter,end_ego_inter,...
   dynamic,t_lim_go,t_lim_wait, [],time_moment); 
pause 

% Velocity (each condition individually)
% for i = 1:4 
% [velocity_table_time_decision_thr_c.(con(i,:))]=...
%  analysis_action_dynamics(table_all(pos.(con(i,:)),:),Log_all,start_ego_inter(pos.(con(i,:))),end_ego_inter(pos.(con(i,:))),...
%    dynamic,t_lim_go,t_lim_wait, condition_4(i,:),time_moment); 
% pause 
% end 







%% Data - Decision behaviour (Figures)
% Average decision behaviour over all paricipants versus individual
% average decision behaviour

% Confidence interval decision (95%) over participants  
load Conf_data_indiv.mat
CI_95_GO = 1.96 * std([GoC1Pr,GoC2Pr,GoC3Pr,GoC4Pr])/sqrt(length(ID));

% Figures 
figure1_decision = figure;  
subplot(1,2,1)
hold on 
errorbar(tta,[Go_c_pr.c1, Go_c_pr.c2],CI_95_GO(1:2),'.-k', 'LineWidth', 1, 'MarkerSize', 20)
plot((ones(length(GoC2Pr),2).*tta)',[GoC1Pr,GoC2Pr]','.--','color', [.7 .7 .7],'MarkerSize',10)
box off; xlim([5.3,6.7])
ax = gca;             ax.FontSize = 12;
title('Distance 70 meter')
legend('Average (CI 95%)','Individual', 'Location', 'northwest')
legend boxoff
xlabel('Time-to-arrival (TTA), s'); ylabel('Probability of go decision')
subplot(1,2,2)
hold on 
errorbar(tta,[Go_c_pr.c3, Go_c_pr.c4],CI_95_GO(3:4),'.-k', 'LineWidth', 1, 'MarkerSize', 20)
plot((ones(length(GoC2Pr),2).*tta)',[GoC3Pr,GoC4Pr]','.--','color', [.7 .7 .7],'MarkerSize',10)
ax = gca;             ax.FontSize = 12;
box off; xlim([5.3,6.7])
title('Distance 90 meter')
xlabel('Time-to-arrival (TTA), s'); ylabel('Probability of go decision')



figure2_decision = figure; 
subplot(2,2,1)
histogram(GoC1Pr,'Normalization','probability')
title('5.5 sec & 70 m')
ylabel ('Proportion'); xlabel ('Percentage')
subplot(2,2,2)
histogram(GoC2Pr,'Normalization','probability')
title('6.5 sec & 70 m')
ylabel ('Proportion'); xlabel ('Percentage')
subplot(2,2,3)
histogram(GoC3Pr,'Normalization','probability')
title('5.5 sec & 90 m')
ylabel ('Proportion'); xlabel ('Percentage')
subplot(2,2,4)
histogram(GoC4Pr,'Normalization','probability')
title('6.5 sec & 90 m')
ylabel ('Proportion'); xlabel ('Percentage')
sgtitle('Decisions Behaviour: Go')


saveas(figure1_decision, fullfile(fname_dist, 'Decision behaviour average and individual.jpg'))
saveas(figure1_decision, fullfile(fname_dist, 'Decision behaviour average and individual.pdf'))
saveas(figure2_decision, fullfile(fname_dist, 'Dribution decision behaviour (conditions).jpg'))  

%% Data - Response times (Figures)
for i = 1:2 
    if i == 1 
        mean_RT_condition = Conf_mean_rt_go_exp; 
        CI_RT = CI_mean_RT_mean_c_go; 
        RT_mean_C = RT_mean_go_C;
        txt = 'Go'; 
    elseif i == 2 
        mean_RT_condition = Conf_mean_rt_wait_exp; 
        CI_RT = CI_mean_RT_mean_c_wait; 
        RT_mean_C = RT_mean_wait_C;
        txt = 'Wait'; 
    end 

observed_rt_fig=figure; 
subplot(1,2,1)
errorbar(tta,mean_RT_condition(1:2),CI_RT(1:2),'o-k', 'LineWidth', 1.7)
hold on 
plot((ones(length(RT_mean_C(:,2)),2).*tta)',RT_mean_C(:,[1,2])','.--','color', [.7 .7 .7],'LineWidth', 1.2,'MarkerSize',10)
ax = gca;             ax.FontSize = 12;
xlim([5.3,6.7]); ylim([1,4]); box off
title('Distance 70 meter')
legend('Average (CI 95%)', 'Individual')
legend boxoff
xlabel('Time-to-arrival (TTA), s')
ylabel('Response time (RT), s')

subplot(1,2,2)
errorbar(tta,mean_RT_condition(3:4),CI_RT(3:4),'o-k', 'LineWidth', 1.7)
hold on 
plot((ones(length(RT_mean_C(:,4)),2).*tta)',RT_mean_C(:,[3,4])','.--','color', [.7 .7 .7],'LineWidth', 1.2,'MarkerSize',10)
ax = gca;             ax.FontSize = 12;
xlim([5.3,6.7]); ylim([1,4]); box off
title('Distance 90 meter')
xlabel('Time-to-arrival (TTA), s'); ylabel('Response time (RT), s')
sgtitle(txt)

saveas(observed_rt_fig, fullfile(fname_dist, ['Response time ',txt,'.jpg']))
saveas(observed_rt_fig, fullfile(fname_dist, ['Response time ',txt,'.pdf']))
end  

figure_conf_RT = figure; 
subplot(1,2,1) 
yyaxis left
hold on; grid on 
errorbar(tta,Conf_mean_rt_go_exp(1:2),CI_mean_RT_mean_c_go(1:2),'.-', 'LineWidth', 0.5, 'MarkerSize', 10)
errorbar(tta,Conf_mean_rt_go_exp(3:4),CI_mean_RT_mean_c_go(3:4),'-s', 'LineWidth',  0.5, 'MarkerSize', 5)
ax = gca;             ax.FontSize = 12;
ylim([1,3]); xlim([5.3,6.7])
xlabel('TTA [sec]'); ylabel('Reaction Time')
yyaxis right 
ylabel('Confidence Judgement')
errorbar(tta,Conf_mean_c_go_exp(1:2),CI_mean_conf_go_exp(1:2),'.--','LineWidth',  0.5,'MarkerSize',10)
errorbar(tta,Conf_mean_c_go_exp(3:4),CI_mean_conf_go_exp(3:4),'--s','LineWidth',  0.5,'MarkerSize',5)
ax = gca;             ax.FontSize = 12;
ylim([1,5])
title ('Go Decisions')

subplot(1,2,2)
yyaxis left 
hold on; grid on 
errorbar(tta,Conf_mean_rt_wait_exp(1:2),CI_mean_RT_mean_c_wait(1:2),'.-', 'LineWidth', 0.5, 'MarkerSize', 10)
errorbar(tta,Conf_mean_rt_wait_exp(3:4),CI_mean_RT_mean_c_wait(3:4),'-s', 'LineWidth',  0.5, 'MarkerSize', 5)
ax = gca;             ax.FontSize = 12;
ylim([1,3]); xlim([5.3,6.7])
xlabel('TTA [sec]'); ylabel('Reaction Time ')
yyaxis right
errorbar(tta,Conf_mean_c_wait_exp(1:2),CI_mean_conf_wait_exp(1:2),'.--','LineWidth',  0.5,'MarkerSize',10)
errorbar(tta,Conf_mean_c_wait_exp(3:4),CI_mean_conf_wait_exp(3:4),'--s','LineWidth',  0.5,'MarkerSize',5)
ax = gca;             ax.FontSize = 12;
ylim([1,5])
ylabel('Confidence Judgement')
legend('RT 70 m (CI 95%)', 'RT 90 m (CI 95%)', 'Conf. 70m (CI 95%)','Conf. 90m (CI 95%)', 'Location', 'southeast')
title ('Wait Decisions')
sgtitle(["Response Times and", "Confidence Judgements"])
saveas(figure_conf_RT, fullfile(fname_dist, 'RT and conf.jpg'))


figure_RT_go_wait = figure;  
subplot(2,1,1)
histogram((table_all.reaction_time(pos_go_all)), 'Normalization','probability')
ax = gca;             ax.FontSize = 12;
title('Go'); xlabel ('time [seconds]'); ylabel ('probability')
xlim([0,6])
subplot(2,1,2)
histogram((table_all.reaction_time(pos_wait_all)), 'Normalization','probability')
ax = gca;             ax.FontSize = 12;
title('Wait'); xlabel ('time [seconds]'); ylabel ('probability')
xlim([0,6])
sgtitle ('Response time ')
saveas(figure_RT_go_wait, fullfile(fname_dist, 'Distribution response time.jpg'))




%% Data - Confidence judgmeents (Figures)
% Average confidence judgements over all participans versus average
% confidence judgmens of individuals 
close all

for i = 1:2 
    if i == 1 
        mean_conf_condition = Conf_mean_c_go_exp; 
        CI_conf = CI_mean_conf_go_exp; 
        conf_mean_C = conf_mean_go_C;
        txt = 'Go'; 
    elseif i == 2 
        mean_conf_condition = Conf_mean_c_wait_exp; 
        CI_conf = CI_mean_conf_wait_exp; 
        conf_mean_C = conf_mean_wait_C;
        txt = 'Wait'; 
    end 

    observed_conf_fig=figure; 
    subplot(1,2,1)
    errorbar(tta,mean_conf_condition(1:2),CI_conf(1:2),'o-k', 'LineWidth', 1.7)
    hold on 
    plot((ones(length(conf_mean_C(:,2)),2).*tta)',conf_mean_C(:,[1,2])','.--','color', [.7 .7 .7],'LineWidth', 1.5,'MarkerSize',10)
    ax = gca;             ax.FontSize = 12;
    xlim([5.3,6.7]); ylim([1,5]); box off 
    title('Distance 70 meter')
    legend('average (CI 95%)', 'individual', 'Location', 'Southwest')
    legend boxoff
    xlabel('Time-to-arrival (TTA), s')
    ylabel('Confidence')

    subplot(1,2,2)
    errorbar(tta,mean_conf_condition(3:4),CI_conf(3:4),'o-k', 'LineWidth', 1.7)
    hold on 
    plot((ones(length(conf_mean_C(:,4)),2).*tta)',conf_mean_C(:,[3,4])','.--','color', [.7 .7 .7],'LineWidth', 1.5,'MarkerSize',10)
    ax = gca;             ax.FontSize = 12;
    xlim([5.3,6.7]); ylim([1,5]); box off
    title('Distance 90 meter')
    xlabel('Time-to-arrival (TTA), s'); ylabel('Confidence ')
    sgtitle(txt)
    saveas(observed_conf_fig, fullfile(fname_dist, ['Confidence ',txt,'.jpg']))
    saveas(observed_conf_fig, fullfile(fname_dist, ['Confidence ',txt,'.pdf']))
end  


% Confidence distribution for each condition 
for i = 1: 2 
    if i == 1 
        plot_pos = Go_c; 
        mean_conf_c = Conf_mean_c_go_exp;
        txt = 'Go decisions (n=17)'; 
    elseif i == 2 
        plot_pos = Wait_c;
        mean_conf_c = Conf_mean_c_wait_exp;
        txt = 'Wait decisions (n=17)'; 
    end 
    
    figure_conf_dist = figure; 
    subplot(2,2,1)
    histogram(Conf_all(plot_pos.c1,4),'Normalization','pdf')
    hold on
    plot(ones(1,2)*mean_conf_c(1),[0,0.6],'LineWidth',1.5)
    plot(0.5:0.1:5.5,pdf(fitdist(Conf_all(plot_pos.c1,4),'Normal'),0.5:0.1:5.5),'LineWidth',1.5)
    ax = gca;             ax.FontSize = 12;
    xlim([1,5]); ylim([0,0.6])
    title('tta 5.5 sec and dist 70 m')
    xlabel('Confidence')
    ylabel('Probability')
    legend('judgement', 'mean','normal dist','Location','northwest')

    subplot(2,2,2)
    histogram(Conf_all(plot_pos.c2,4),'Normalization','pdf')
    hold on
    plot(ones(1,2)*mean_conf_c(2),[0,0.6],'LineWidth',1.5)
    plot(0.5:0.1:5.5,pdf(fitdist(Conf_all(plot_pos.c2,4),'Normal'),0.5:0.1:5.5),'LineWidth',1.5)
    ax = gca;             ax.FontSize = 12;
    xlim([1,5]); ylim([0,0.6])
    title('tta 6.5 sec and dist 70 m')
    xlabel('Confidence')
    ylabel('Probability')

    subplot(2,2,3)
    histogram(Conf_all(plot_pos.c3,4),'Normalization','pdf')
    hold on
    plot(ones(1,2)*mean_conf_c(3),[0,0.6],'LineWidth',1.5)
    plot(0.5:0.1:5.5,pdf(fitdist(Conf_all(plot_pos.c3,4),'Normal'),0.5:0.1:5.5),'LineWidth',1.5)
    ax = gca;             ax.FontSize = 12;
    xlim([1,5]); ylim([0,0.6])
    title('tta 5.5 sec and dist 90 m')
    xlabel('Confidence')
    ylabel('Probability')

    subplot(2,2,4)
    histogram(Conf_all(plot_pos.c4,4),'Normalization','pdf')
    hold on
    plot(ones(1,2)*mean_conf_c(4),[0,0.6],'LineWidth',1.5)
    plot(0.5:0.1:5.5,pdf(fitdist(Conf_all(plot_pos.c4,4),'Normal'),0.5:0.1:5.5),'LineWidth',1.5)
    ax = gca;             ax.FontSize = 12;
    xlim([1,5]); ylim([0,0.6])
    title('tta 6.5 sec and dist 90 m')
    xlabel('Confidence')
    ylabel('Probability')

    sgtitle(['Confidence Judgements: ', txt])
    saveas(figure_conf_dist, fullfile(fname_dist, ['Distribution confidence ', txt,'.jpg']))
    close all
end 

figure_conf_dec = figure; 
subplot(1,2,1) 
yyaxis left
hold on; grid on 
plot(tta,[Go_c_pr.c1, Go_c_pr.c2],'.-', 'LineWidth', 0.5, 'MarkerSize', 10)
plot(tta,[Go_c_pr.c3, Go_c_pr.c4],'-s', 'LineWidth', 0.5, 'MarkerSize', 5)
ax = gca;             ax.FontSize = 12;
ylim([0,1]); xlim([5.3,6.7])
xlabel('TTA [sec]'); ylabel('Probability of Going')
yyaxis right 
ylabel('Confidence Judgement')
plot(tta,Conf_mean_c_go_exp(1:2),'.--','LineWidth', 0.5,'MarkerSize',10)
plot(tta,Conf_mean_c_go_exp(3:4),'--s','LineWidth', 0.5,'MarkerSize',5)
ylim([1,5])
legend('Pr. 70 m', 'Pr. 90 m', 'Conf. 70m','Conf. 90m', 'Location', 'northwest')
title ('Go Decisions')

subplot(1,2,2)
yyaxis left 
hold on; grid on 
plot(tta,[1-Go_c_pr.c1, 1-Go_c_pr.c2],'.-', 'LineWidth', 0.5, 'MarkerSize', 10)
plot(tta,[1-Go_c_pr.c3, 1-Go_c_pr.c4],'-s', 'LineWidth', 0.5, 'MarkerSize', 5)
ax = gca;             ax.FontSize = 12;
ylim([0,1]); xlim([5.3,6.7])
xlabel('TTA [sec]'); ylabel('Probability of Going')
yyaxis right
plot(tta,Conf_mean_c_wait_exp(1:2),'.--','LineWidth', 0.5,'MarkerSize',10)
plot(tta,Conf_mean_c_wait_exp(3:4),'--s','LineWidth', 0.5,'MarkerSize',5)
ylim([1,5])
ylabel('Confidence Judgement')
title ('Wait Decisions')
sgtitle(["Decision Behaviour and", "Confidence Judgements"])
saveas(figure_conf_dec, fullfile(fname_dist, 'decision and confidence.jpg'))

%% Data - First throttle & turn time (Figures)
close all
fig_response_dist = figure;
subplot(2,1,1)
histogram(table_all.first_thr(pos_go_all)-table_all.reaction_time(pos_go_all),'Normalization','probability')
title('Go'); xlabel('time, s'); ylabel('Proportion')
ax = gca;             ax.FontSize = 12;
ylim([0,0.42])
box off
subplot(2,1,2)
histogram(table_all.first_thr(pos_wait_all)-table_all.reaction_time(pos_wait_all),'Normalization','probability')
title('Wait'); xlabel('time, s'); ylabel('Proportion')
ax = gca;             ax.FontSize = 12;
box off
ylim([0,0.42])
saveas(fig_response_dist, fullfile(fname_dist, 'Response time distribution.jpg'))
saveas(fig_response_dist, fullfile(fname_dist, 'Response time distribution.pdf'))


figure_time_gap = figure;  
subplot(2,1,1)
histogram(table_all.turn_time(pos_go_all), 'Normalization','probability')
ax = gca;             ax.FontSize = 12;
hold on 
title('Go'); xlim([0,20])
xlabel ('time [seconds]');     ylabel ('probability')
subplot(2,1,2)
histogram(table_all.turn_time(pos_wait_all), 'Normalization','probability')
ax = gca;             ax.FontSize = 12;
hold on 
title('Wait'); xlim([0,20])
xlabel ('time [seconds]');     ylabel ('probability')
sgtitle ('Time between decision and confidence judgment all conditions')
saveas(figure_time_gap, fullfile(fname_dist, 'Turn time.jpg'))



%% Save: CI decision, mean confidence and response times 

save([fname_output, '\CI_decision.mat'],'CI_95_GO')

Conf_mean_c_go = Conf_mean_c_go_exp; Conf_mean_c_wait = Conf_mean_c_wait_exp; 
save([fname_output, '\Confidence_mean_all.mat'], 'Conf_mean_c_go', 'Conf_mean_c_wait')

RT_mean_c_go = Conf_mean_rt_go_exp; RT_mean_c_wait = Conf_mean_rt_wait_exp;
save([fname_output, '\Reaction_times_mean_all.mat'], 'RT_mean_c_go', 'RT_mean_c_wait', 'CI_mean_RT_mean_c_go', 'CI_mean_RT_mean_c_wait')

%% Export data: for parameter tuning confidence model 
close all
save_name = '.\Modelling\data\'; 

% All data 
is_go_decision = ones(length(Conf_all(:,1)),1);
is_go_decision(Wait)=2; 
is_go_decision=categorical(is_go_decision,[1 2],{'True' 'False'});
is_go_decision([change_of_mind,no_button])=[]; 

num_decision = ones(length(Conf_all(:,1)),1);
num_decision(Wait)=2; 
num_decision([change_of_mind,no_button])=[];

data_output= table(table_all.ID,table_all.reaction_time,is_go_decision,table_all.tta, table_all.dist,num_decision,table_all.conf,...
    'VariableNames',{'subj_id','RT','is_go_decision','tta_condition','d_condition','num_decision','confidence'});  
writetable(data_output, [save_name, 'data_output.csv'])


% Throttle as response time for "Go" decisions 
data_output_thr = data_output; 
data_output_thr.RT(data_output_thr.is_go_decision == 'True') = time_first_throttle_use(data_output_thr.is_go_decision == 'True');
writetable(data_output_thr, [save_name, 'data_output_thr.csv'])



% Exclude extreem RT's 
data_output_RT = data_output;
prs_lef_out_go_min   = length(find(data_output_RT.is_go_decision == 'True' & data_output_RT.RT < 1))/length(find(data_output_RT.is_go_decision == 'True'))*100;
prs_lef_out_wait_min = length(find(data_output_RT.is_go_decision == 'False' & data_output_RT.RT < 1))/length(find(data_output_RT.is_go_decision == 'False'))*100;
prs_lef_out_go_max   = length(find(data_output_RT.is_go_decision == 'True' & data_output_RT.RT > 2.5))/length(find(data_output_RT.is_go_decision == 'True'))*100;
prs_lef_out_wait_max = length(find(data_output_RT.is_go_decision == 'False' & data_output_RT.RT > 4))/length(find(data_output_RT.is_go_decision == 'False'))*100;

data_output_RT(data_output_RT.RT < 1,:) =[];
data_output_RT(data_output_RT.is_go_decision == 'True' & data_output_RT.RT > 2.5,:) =[];
data_output_RT(data_output_RT.is_go_decision == 'False'& data_output_RT.RT > 4,:) =[];

writetable(data_output_RT, [save_name, 'data_output_RT.csv'])

% Confidence and Response times excl. extreem  
table_excl_ID=table_all; 
table_excl_ID(table_excl_ID.decision == 'Go' & table_excl_ID.reaction_time > 2.5, :) = [];
table_excl_ID(table_excl_ID.decision == 'Wait' & table_excl_ID.reaction_time > 4, :) = []; 
table_excl_ID(table_excl_ID.reaction_time < 1, :) = []; 

for i = 1: 4 
    pos_go.(con(i,:)) = find(table_excl_ID.tta== tta_con(i) &  table_excl_ID.dist==dist_con(i) & table_excl_ID.decision=='Go');
    pos_wait.(con(i,:)) = find(table_excl_ID.tta== tta_con(i) &  table_excl_ID.dist==dist_con(i) & table_excl_ID.decision=='Wait');
    
    std_RT_c_go_all(i) = std(table_excl_ID.reaction_time(pos_go.(con(i,:)))); 
    std_RT_c_wait_all(i) = std(table_excl_ID.reaction_time(pos_wait.(con(i,:)))); 
    
    RT_mean_c_go(i) = mean(table_excl_ID.reaction_time(pos_go.(con(i,:)))); 
    RT_mean_c_wait(i) = mean(table_excl_ID.reaction_time(pos_wait.(con(i,:))));
    
    CI_mean_RT_mean_c_go(i) = function_CI_95(table_excl_ID.reaction_time, pos_go.(con(i,:))', 1); 
    CI_mean_RT_mean_c_wait(i) = function_CI_95(table_excl_ID.reaction_time, pos_wait.(con(i,:))', 1); 
    
    CI_conf_mean_c_go(i) = function_CI_95(table_excl_ID.conf, pos_go.(con(i,:))', 1);
    CI_conf_mean_c_wait(i) = function_CI_95(table_excl_ID.conf, pos_wait.(con(i,:))', 1);
    
    Conf_mean_c_go(i) = mean(table_excl_ID.conf(pos_go.(con(i,:)))); 
    Conf_mean_c_wait(i) = mean(table_excl_ID.conf(pos_wait.(con(i,:))));
end 


save([fname_output,'\Reaction_times_mean_RT.mat'], 'RT_mean_c_go', 'RT_mean_c_wait', 'CI_mean_RT_mean_c_go', 'CI_mean_RT_mean_c_wait','std_RT_c_go_all', 'std_RT_c_wait_all')
save([fname_output,'\Confidence_mean_RT.mat'], 'Conf_mean_c_go', 'Conf_mean_c_wait', 'CI_conf_mean_c_go', 'CI_conf_mean_c_wait')












