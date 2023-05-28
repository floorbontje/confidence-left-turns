function [random_effect_inter,random_effect_wait] = function_lme_results (txt, lme_model, pos_go, pos_wait, mean_go_exp, CI_mean_go_exp, ...
    mean_wait_exp, CI_mean_wait_exp)


% Fit model 
fit_lme  = fitted(lme_model);
con = ['c1';'c2';'c3';'c4'];

for i = 1:4 
    Conf_mean_c_go_fit(i) = mean(fit_lme(pos_go.(con(i,:))));
    CI_mean_go_fit(i) = 1.96 * std(fit_lme(pos_go.(con(i,:))))/length(pos_go.(con(i,:))); 
    if ~isempty(mean_wait_exp)
        Conf_mean_c_wait_fit(i) = mean(fit_lme(pos_wait.(con(i,:))));
        CI_mean_wait_fit(i) = 1.96 * std(fit_lme(pos_wait.(con(i,:))))/length(pos_wait.(con(i,:))); 
    end
end 

% Random effects - Confidence  (RT)
[~, ~, random_coef]=randomEffects(lme_model);
if ~isempty(mean_wait_exp)
    pos_effect          = find(random_coef.pValue <0.05);           %Positions significant effect 
    pos_wait_effect     = rem(pos_effect,2).*pos_effect;            %Positions: decision
    pos_inter_effect    = ~rem(pos_effect,2).*pos_effect;           %Positions: intercept
    pos_wait_effect(pos_wait_effect==0)     = [];
    pos_inter_effect(pos_inter_effect ==0)  = []; 

    random_inter = random_coef.Estimate(pos_inter_effect);          %Value effect intercept
    random_wait  = random_coef.Estimate(pos_wait_effect);           %Value effect decision  

    random_effect_inter = random_coef(pos_inter_effect,[4,8]);      %Table intercept: effect + pValue
    random_effect_wait  = random_coef(pos_wait_effect,[4,8]);       %Table decision: effect + pValue
else 
    [~, ~, random_coef]=randomEffects(lme_model);

    pos_inter_effect = find(random_coef.pValue <0.05);                  %Positions significant effect 
    random_inter = random_coef.Estimate(pos_inter_effect);              %Value effect intercept
    random_effect_inter = random_coef(pos_inter_effect,[4,8]);          %Table intercept: effect + pValue
    random_effect_wait = []; 
    random_wait = []; 
end 


% Figure model 

y_min = round(min([mean_wait_exp,mean_go_exp]),1)-0.5;
y_max = round(max([mean_wait_exp,mean_go_exp]),1)+0.5;

figure_model = figure;
if ~isempty(mean_wait_exp)
    subplot(1,2,1)
    errorbar([5.5;6.5],Conf_mean_c_go_fit(1:2),CI_mean_go_fit(1:2),'.-b','LineWidth',1.3,'MarkerSize',15)
    hold on
    errorbar([5.5;6.5],mean_go_exp(1:2),CI_mean_go_exp(1:2) ,'.--','color', [.73 .73 1],'LineWidth',1.3,'MarkerSize',15)
    errorbar([5.5;6.5],Conf_mean_c_go_fit(3:4),CI_mean_go_fit(3:4),'.-r','LineWidth',1.3,'MarkerSize',15)
    errorbar([5.5;6.5],mean_go_exp(3:4),CI_mean_go_exp(3:4) ,'.--','color', [1 .73 .73],'LineWidth',1.3, 'MarkerSize',15)
    box off;  title ('Go');
    if Conf_mean_c_wait_fit(4)>Conf_mean_c_go_fit(4)
        legend ('Model 70 m','Data 70m','Model 90 m', 'Data 90m','Location','northwest')
        legend boxoff
    end
    ylabel (txt)
    xlabel ('Time-to-arrival (TTA), s')
    xlim([5.3 6.7])
    ylim([y_min,y_max])

    subplot(1,2,2)
    errorbar([5.5;6.5],Conf_mean_c_wait_fit(1:2),CI_mean_wait_fit(1:2),'.-b','LineWidth',1.3,'MarkerSize',15)
    hold on
    errorbar([5.5;6.5],Conf_mean_c_wait_fit(3:4),CI_mean_wait_fit(3:4),'.-r','LineWidth',1.3,'MarkerSize',15)
    errorbar([5.5;6.5],mean_wait_exp(1:2),CI_mean_wait_exp(1:2) ,'.--','color', [.73 .73 1],'LineWidth',1.3,'MarkerSize',15)
    errorbar([5.5;6.5],mean_wait_exp(3:4),CI_mean_wait_exp(3:4) ,'.--','color', [1 .73 .73], 'LineWidth',1.3,'MarkerSize',15)
    if Conf_mean_c_wait_fit(4)<Conf_mean_c_go_fit(4)
        legend ('Model 70 m','Data 70m','Model 90 m', 'Data 90m','Location','northeast')
        legend boxoff
    end
    title ('Wait'); box off; 
    ylabel (txt)
    xlabel ('Time-to-arrival (TTA), s')
    xlim([5.3 6.7])
    ylim([y_min,y_max])
else 
    hold on; box off
    plot([5.5;6.5], Conf_mean_c_go_fit(1:2),'.-b','LineWidth',1.3,'MarkerSize', 20)
    plot([5.5;6.5], mean_go_exp(1:2),'.--','color', [.73 .73 1],'LineWidth',1.3,'MarkerSize', 20)
    plot([5.5;6.5], Conf_mean_c_go_fit(3:4),'.-r', 'LineWidth',1.3,'MarkerSize', 20)
    plot([5.5;6.5], mean_go_exp (3:4),'.--', 'color', [1 .73 .73],'LineWidth',1.3,'MarkerSize', 20)
    xlim([5.3,6.7]); ylim([0,1])
    xlabel('Time-to-arrival (TTA), s'); ylabel(txt)
    legend('Model 70 m', 'Data 70 m', 'Model 90m', 'Data 90m','Location', 'northwest')
    legend boxoff
end 


% Figure random effects 
figure_random_lme = figure; 
if ~isempty(mean_wait_exp)
    subplot(1,2,1)
    histogram(random_inter,'Normalization','probability')
    ylim([0,1])
    xlabel('Estimate: Intercept'); ylabel('Probability')
    subplot(1,2,2)
    histogram(random_wait,'Normalization','probability')
    sgtitle('Random Effects Coefficients') 
    ylim([0,1])
    xlabel('Estimate: Wait Decision'); ylabel('Probability')
else 
    histogram(random_inter,'Normalization','probability')
    ylim([0,1])
    xlabel('Random Intercept'); ylabel('Porportion')
end 


%% save 
fname_lme_analysis  = '.\Data anlysis\general data analysis\lme models'; 
fname_lme  = '.\Data analysis\figures\analysis'; 

% Fixed coefficients 
writetable(dataset2table(lme_model.Coefficients), fullfile(fname_lme_analysis, ['Fixed coefficients lme ', txt, '.csv'])); 

% Random coefficients 
writetable(dataset2table(random_effect_inter), fullfile(fname_lme_analysis, ['Random Intercept lme ', txt, '.csv'])); 
if ~isempty(mean_wait_exp)
    writetable(dataset2table(random_effect_wait), fullfile(fname_lme_analysis, ['Random wait lme ', txt, '.csv'])); 
end 

% Figures 
saveas(figure_model, fullfile(fname_lme, ['lme ', txt, '.jpg']))
saveas(figure_model, fullfile(fname_lme, ['lme ', txt, '.pdf']))
saveas(figure_random_lme, fullfile(fname_lme, ['random effects ',txt, '.jpg']))






