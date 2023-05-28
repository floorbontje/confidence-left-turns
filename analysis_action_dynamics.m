function [table_time_decision]=...
    analysis_action_dynamics(table_all,Log_all,start_ego_inter,end_ego_inter,...
   dynamic,t_lim_go,t_lim_wait, condition,time_moment)


% dynamic
% 1: velocity
% 2: distance to intersection 
% 3: throttle 

if dynamic ==1 
    text1 = 'velocity ';
    text = 'velocity, ';
    unit = 'm/s';
    y_lim_dynamic = [0,10.5];
elseif dynamic == 2 
    y_lim_dynamic = [2,9];
    text1 = 'distance to center of intersection ';
    text = 'distance, ';
    unit = 'm';
elseif dynamic == 3 
    text1 = 'throttle input ';
    text = 'Throttle input '; 
    unit = '[0 - 1]';
    y_lim_dynamic=[0,1.1]; 
    t_lim_go = [-1, t_lim_go(2)]; 
    t_lim_wait = [-1, t_lim_wait(2)];
end 

pos_go = find(table_all.decision == 'Go');  
pos_wait = find(table_all.decision == 'Wait');



for i = 1 : length(start_ego_inter)
    if time_moment == 1 %first throttle use 
        start_in(i) = find(Log_all(start_ego_inter(i): end_ego_inter(i),29)>0,1); %first value throttle > 0 
        text_time = ' (w.r.t gass throttle)'; 
    elseif time_moment == 0 %button press
        t_raw = Log_all(start_ego_inter(i):end_ego_inter(i),7) - Log_all(start_ego_inter(i),7);
        start_in(i)=find(t_raw>table_all.reaction_time(i),1);
        text_time = ' (w.r.t button press)'; 
    end 
end 


pos_in = (max(start_in)-start_in)'+1;
length_input=end_ego_inter-start_ego_inter;
max_length= max(length_input+pos_in); 

Matrix_all=nan(max_length+1,length(start_ego_inter));
%The time interval isn't consistent, so each turn has to have it's own time
%vector where t = 0 is at the decision moment (RT)
t_all = Matrix_all;
 
%Storage of the action dynamics, aligned with the decision moment (RT)
M_all = Matrix_all; 


for i = 1: length(start_ego_inter)
    t_all(pos_in(i):length_input(i)+pos_in(i),i)= Log_all(start_ego_inter(i):end_ego_inter(i),7) - Log_all(start_ego_inter(i)+start_in(i),7); 
    if dynamic == 1 
        M_all(pos_in(i):length_input(i)+pos_in(i),i)= sqrt((Log_all(start_ego_inter(i):end_ego_inter(i),17)).^2+(Log_all(start_ego_inter(i):end_ego_inter(i),18)).^2); 
    elseif dynamic == 2 
        M_all(pos_in(i):length_input(i)+pos_in(i),i) = Log_all(start_ego_inter(i): end_ego_inter(i),8);
    elseif dynamic == 3 
        M_all(pos_in(i):length_input(i)+pos_in(i),i) = Log_all(start_ego_inter(i): end_ego_inter(i),29);
    end 
end 

%% Max/ min/mean value (within time frame)
for i = 1 : length(start_ego_inter)
    if find(i == pos_go)
        t_lim_dec = t_lim_go; 
    elseif find(i == pos_wait)
        t_lim_dec = t_lim_wait; 
    else 
        disp("No position found")
    end 
      
      
    
    if isempty(find(t_all(:,i) > t_lim_dec(1),1))
       pos_t_min(i) = find(t_all(:,i) == min(t_all(:,i))); 
    else 
        pos_t_min(i)= find(t_all(:,i) > t_lim_dec(1),1); 
    end 
    if isempty(find(t_all(:,i) > t_lim_dec(2),1))
       pos_t_max(i) = find(t_all(:,i) == max(t_all(:,i))); 
    else 
        pos_t_max(i)= find(t_all(:,i) > t_lim_dec(2),1)-1;
    end 
 
    max_value_decision(i) = max(M_all(pos_t_min(i):pos_t_max(i),i));   
    min_value_decision(i) = min(M_all(pos_t_min(i):pos_t_max(i),i));
    mean_value_decision(i) = mean(M_all(pos_t_min(i):pos_t_max(i),i),'omitnan');   
     
end 

table_all.max_value = max_value_decision';
table_all.min_value = min_value_decision'; 


[corr_max_go_decision, p_max_go_decision] = corrcoef(table_all.conf(pos_go),max_value_decision(pos_go)); 
[corr_max_wait_decision, p_max_wait_decision] = corrcoef(table_all.conf(pos_wait),max_value_decision(pos_wait));

[corr_min_go_decision, p_min_go_decision] = corrcoef(table_all.conf(pos_go),min_value_decision(pos_go)); 
[corr_min_wait_decision, p_min_wait_decision] = corrcoef(table_all.conf(pos_wait),min_value_decision(pos_wait));

[corr_mean_go_decision, p_mean_go_decision] = corrcoef(table_all.conf(pos_go),mean_value_decision(pos_go)); 
[corr_mean_wait_decision, p_mean_wait_decision] = corrcoef(table_all.conf(pos_wait),mean_value_decision(pos_wait));


%% Difference w.r.t. average 
diff_mean = Matrix_all; 
diff_mean(:,pos_go)=M_all(:,pos_go)-mean(M_all(:,pos_go)','omitnan')'; 
diff_mean(:,pos_wait)=M_all(:,pos_wait)-mean(M_all(:,pos_wait)','omitnan')'; 


%% Difference with respect to mean value of an individual 
difference_mean_indiv = Matrix_all; 
nbr = unique(table_all.ID);

for i = 1:length(nbr)
   j_g = find(table_all.ID== nbr(i) & table_all.decision == 'Go');
   mean_nbr_g = mean(M_all(:,j_g)','omitnan')';
   difference_mean_indiv(:,j_g)= M_all(:,j_g)- mean_nbr_g; 

   j_w = find(table_all.ID== nbr(i) & table_all.decision == 'Wait');
   mean_nbr_w = mean(M_all(:,j_w)','omitnan')';
   difference_mean_indiv(:,j_w)= M_all(:,j_w) - mean_nbr_w; 
end 


%% Mean error (individual and general) and RMSD
for i = 1: length(start_ego_inter)
    pos_min = pos_t_min(i); 
    pos_max = pos_t_max(i); 
    
    if find(i == pos_go) 
        DM_go=sum(difference_mean_indiv(pos_min:pos_max,pos_go), 'omitnan')./length(pos_min:pos_max);
        RMSD_go=sqrt(sum(difference_mean_indiv(pos_min:pos_max,pos_go).^2, 'omitnan')./length(pos_min:pos_max));
        DM_mean_go = sum(diff_mean(pos_min:pos_max,pos_go), 'omitnan')./length(pos_min:pos_max); 
    else 
        DM_wait=sum(difference_mean_indiv(pos_min:pos_max,pos_wait), 'omitnan')./length(pos_min:pos_max);
        RMSD_wait=sqrt(sum(difference_mean_indiv(pos_min:pos_max,pos_wait).^2, 'omitnan')./length(pos_min:pos_max));
        DM_mean_wait = sum(diff_mean(pos_min:pos_max,pos_wait), 'omitnan')./length(pos_min:pos_max); 
    end 
end 
 
table_all.DM_indiv(pos_go) = DM_go'; 
table_all.DM_indiv(pos_wait) = DM_wait'; 
table_all.DM_all(pos_go) = DM_mean_go'; 
table_all.DM_all(pos_wait) = DM_mean_wait'; 
table_all.RMSD(pos_go) = RMSD_go'; 
table_all.RMSD(pos_wait) = RMSD_wait';

[DM_corr_go, DM_p_go] = corrcoef(table_all.conf(pos_go),DM_go);
[RMSD_corr_go, RMSD_p_go] = corrcoef(table_all.conf(pos_go),RMSD_go);
[DM_mean_corr_go, DM_mean_p_go] = corrcoef(table_all.conf(pos_go),DM_mean_go); 

[DM_corr_wait, DM_p_wait] = corrcoef(table_all.conf(pos_wait),DM_wait);
[RMSD_corr_wait, RMSD_p_wait] = corrcoef(table_all.conf(pos_wait),RMSD_wait);
[DM_mean_corr_wait, DM_mean_p_wait] = corrcoef(table_all.conf(pos_wait), DM_mean_wait);  

%% Correlations
names = {'max_value', 'p_max_value', 'min_value', 'p_min_value', ... 
    'mean_value','p_mean_value', 'DM_indiv','p_DM_indiv', 'RMSD','p_RMSD','DM_mean', 'p_DM_mean'};

table_time_decision = table(names',[corr_max_go_decision(1,2), p_max_go_decision(1,2),...
    corr_min_go_decision(1,2), p_min_go_decision(1,2),...
    corr_mean_go_decision(1,2), p_mean_go_decision(1,2),...
    DM_corr_go(1,2), DM_p_go(1,2), ...
    RMSD_corr_go(1,2), RMSD_p_go(1,2),...
    DM_mean_corr_go(1,2), DM_mean_p_go(1,2)]',...
    [corr_max_wait_decision(1,2), p_max_wait_decision(1,2),...
    corr_min_wait_decision(1,2), p_min_wait_decision(1,2),...
    corr_mean_wait_decision(1,2), p_mean_wait_decision(1,2),...
    DM_corr_wait(1,2), DM_p_wait(1,2), ...
    RMSD_corr_wait(1,2), RMSD_p_wait(1,2),...
    DM_mean_corr_wait(1,2), DM_mean_p_wait(1,2)]','VariableNames',{'Variable','Go','Wait'});


%% Statistical analysis 
lme_max=fitlme(table_all, 'max_value ~ conf * decision + (1|ID)','CheckHessian', true)
lme_min=fitlme(table_all, 'min_value ~ conf * decision + (1|ID)','CheckHessian', true)
lme_DM_indiv = fitlme(table_all, 'DM_indiv ~ conf*decision + (1|ID)','CheckHessian', true)
lme_DM_all=fitlme(table_all, 'DM_all ~ conf*decision+ (1|ID)','CheckHessian', true)
lme_RMSD=fitlme(table_all, 'RMSD ~ conf*decision+ (1|ID)','CheckHessian', true)


[table_decision_max, table_metric_max, F_crit1, F_crit2]=Post_hoc_analysis('max_value', table_all) 
[table_decision_min, table_metric_min, ~, ~]=Post_hoc_analysis('min_value', table_all) 
[table_decision_DM_indiv, table_metric_DM_indiv, ~, ~]=Post_hoc_analysis('DM_indiv', table_all) 
[table_decision_DM_all, table_metric_DM_all, ~, ~]=Post_hoc_analysis('DM_all', table_all) 
[table_decision_RMSD, table_metric_RMSD, ~, ~]=Post_hoc_analysis('RMSD', table_all) 

%% Figure - mean values 

[conf_sort_go, I_sort_go]=sort(table_all.conf(pos_go));
[conf_sort_wait, I_sort_wait]=sort(table_all.conf(pos_wait));
int_go = round(length(conf_sort_go)/3)*[1,2]; 
int_wait = round(length(conf_sort_wait)/3)*[1,2]; 

Conf_low_pos_go = pos_go(I_sort_go(1:int_go(1))); 
Conf_mean_pos_go = pos_go(I_sort_go(int_go(1)+1:int_go(2)));
Conf_high_pos_go = pos_go(I_sort_go(int_go(2)+1:end));
Conf_low_pos_wait = pos_wait(I_sort_wait(1:int_wait(1))); 
Conf_mean_pos_wait = pos_wait(I_sort_wait(int_wait(1)+1:int_wait(2)));
Conf_high_pos_wait = pos_wait(I_sort_wait(int_wait(2)+1:end));


M_mean_go = M_all(:,Conf_mean_pos_go); 
M_low_go = M_all(:,Conf_low_pos_go); 
M_high_go = M_all(:,Conf_high_pos_go);

t_mean_go = t_all(:,Conf_mean_pos_go); 
t_low_go = t_all(:,Conf_low_pos_go); 
t_high_go = t_all(:,Conf_high_pos_go);

CI_mean_go = nan(1,max_length+1); 
CI_low_go = nan(1,max_length+1); 
CI_high_go = nan(1,max_length+1); 

for i = 1:(max_length+1)
    if ~isempty(find(~isnan(M_mean_go(i,:)), 1))
        CI_mean_go(i) = 1.96*std(M_mean_go(i,:),'omitnan')./sqrt(length(find(~isnan(M_mean_go(i,:))))); 
    end 
    if ~isempty(find(~isnan(M_low_go(i,:)), 1))
        CI_low_go(i) = 1.96*std(M_low_go(i,:),'omitnan')./sqrt(length(find(~isnan(M_low_go(i,:)))));
    end 
    if ~isempty(find(~isnan(M_high_go(i,:)), 1))
        CI_high_go(i) = 1.96*std(M_high_go(i,:),'omitnan')./sqrt(length(find(~isnan(M_high_go(i,:)))));
    end 
end 


M_mean_wait = M_all(:,Conf_mean_pos_wait); 
M_low_wait = M_all(:,Conf_low_pos_wait); 
M_high_wait = M_all(:,Conf_high_pos_wait);

t_mean_wait = t_all(:,Conf_mean_pos_wait); 
t_low_wait = t_all(:,Conf_low_pos_wait); 
t_high_wait = t_all(:,Conf_high_pos_wait);

CI_mean_wait = nan(1,max_length+1); 
CI_low_wait = nan(1,max_length+1); 
CI_high_wait = nan(1,max_length+1); 

for i = 1:(max_length+1)
    if ~isempty(find(~isnan(M_mean_wait(i,:)), 1))
        CI_mean_wait(i) = 1.96*std(M_mean_wait(i,:),'omitnan')./sqrt(length(find(~isnan(M_mean_wait(i,:))))); 
    end 
    if ~isempty(find(~isnan(M_low_wait(i,:)), 1))
        CI_low_wait(i) = 1.96*std(M_low_wait(i,:),'omitnan')./sqrt(length(find(~isnan(M_low_wait(i,:)))));
    end 
    if ~isempty(find(~isnan(M_high_wait(i,:)), 1))
        CI_high_wait(i) = 1.96*std(M_high_wait(i,:),'omitnan')./sqrt(length(find(~isnan(M_high_wait(i,:)))));
    end 
end 

mean_fig =figure; 
subplot(1,2,1)
plot(mean(t_low_go','omitnan'),mean(M_low_go','omitnan'),'-','color',[0.8,0.8,0.8],'LineWidth', 1.3);
hold on; 
plot(mean(t_mean_go','omitnan'),mean(M_mean_go','omitnan'),'-','color',[0.6,0.6,0.6],'LineWidth', 1.2);
plot(mean(t_high_go','omitnan'),mean(M_high_go','omitnan'),'-k','LineWidth', 1.2);
% plot(mean(t_low_go','omitnan'), mean(M_low_go','omitnan') + CI_low_go,'--','color',[0.8,0.8,0.8],'LineWidth', 1.3);
% plot(mean(t_mean_go','omitnan'),mean(M_mean_go','omitnan')+ CI_mean_go,'-','color',[0.6,0.6,0.6],'LineWidth', 1.2);
% plot(mean(t_high_go','omitnan'),mean(M_high_go','omitnan') + CI_high_go,'-k','LineWidth', 1.2);
% plot(mean(t_low_go','omitnan'), mean(M_low_go','omitnan') - CI_low_go,'--','color',[0.8,0.8,0.8],'LineWidth', 1.3);
% plot(mean(t_mean_go','omitnan'),mean(M_mean_go','omitnan')- CI_mean_go,'--','color',[0.6,0.6,0.6],'LineWidth', 1.2);
% plot(mean(t_high_go','omitnan'),mean(M_high_go','omitnan') - CI_high_go,'--k','LineWidth', 1.2);
ax = gca;          ax.FontSize = 13;
title('Go', 'FontSize', 14); box off
ylabel([text unit])
xlabel('time, s')
xlim(t_lim_go)
ylim(y_lim_dynamic)


subplot(1,2,2)
plot(mean(t_low_wait','omitnan'),mean(M_low_wait','omitnan'),'-','color',[0.8,0.8,0.8],'LineWidth',1.3);
hold on;
plot(mean(t_mean_wait','omitnan'),mean(M_mean_wait','omitnan'),'-','color',[0.6,0.6,0.6],'LineWidth', 1.2);
plot(mean(t_high_wait','omitnan'),mean(M_high_wait','omitnan'),'-k','LineWidth',1.2);
ax = gca;          ax.FontSize = 12;
legend('Lowest conf', 'Intermediate conf', 'Highest conf', 'Fontsize', 12)
legend boxoff; box off
title('Wait', 'FontSize', 14)
ylabel([text unit])
xlabel('time, s')
xlim(t_lim_wait)
ylim(y_lim_dynamic)



%% Output 
fname_output = '.\Data analysis\general data analysis\action dynamics'; 
fname_figure = '.\Data analysis\figures\Action Dynamics';

if time_moment == 1
    txt_time = 'Throttle'; 
else 
    txt_time = 'Button'; 
end 

% analysis 
writetable(table_time_decision, fullfile(fname_output, ['Correlations ', text1, condition, txt_time, '.csv']))
writetable(dataset2table(lme_max.Coefficients), fullfile(fname_output, ['lme max ', text1, condition, txt_time, '.csv']))
writetable(dataset2table(lme_min.Coefficients), fullfile(fname_output, ['lme min ', text1, condition, txt_time, '.csv']))
writetable(dataset2table(lme_DM_indiv.Coefficients), fullfile(fname_output, ['lme DM indiv ', text1, condition, txt_time, '.csv']))
writetable(dataset2table(lme_DM_all.Coefficients), fullfile(fname_output, ['lme DM general ', text1, condition, txt_time, '.csv']))
writetable(dataset2table(lme_RMSD.Coefficients), fullfile(fname_output, ['lme RMSD ', text1, condition, txt_time, '.csv']))

writetable(table_decision_max, fullfile(fname_output, ['Ftest Max ', text1, condition, txt_time, '.csv']))
writetable(table_decision_min, fullfile(fname_output, ['Ftest Min ', text1, condition, txt_time, '.csv']))
writetable(table_decision_DM_indiv, fullfile(fname_output, ['Ftest DM indiv ', text1, condition, txt_time, '.csv']))
writetable(table_decision_DM_all, fullfile(fname_output, ['Ftest DM general ', text1, condition, txt_time, '.csv']))
writetable(table_decision_RMSD, fullfile(fname_output, ['Ftest DM RMSD ', text1, condition, txt_time, '.csv']))


% figures
saveas(mean_fig, fullfile(fname_figure, ['Mean ', text1, 'during turn', condition, txt_time, '.jpg']))
saveas(mean_fig, fullfile(fname_figure, ['Mean ', text1, 'during turn', condition, txt_time, '.pdf']))


