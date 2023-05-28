function [table_decision, table_metric, F_crit_1dec, F_crit_2dec]=Post_hoc_analysis(metric, table_all)

lme = fitlme(table_all, [metric,' ~ conf*decision+ (1|ID)'],'CheckHessian', true);


H_conf_go = [eye(2,2), zeros(2,2)]; 
H_conf_wait = [eye(2,2), eye(2,2)]; 

for i = 1:2
[pVal_go(i), F_go(i), DF1_1dec, DF2_1dec] = coefTest(lme, H_conf_go(i,:));
[pVal_wait(i), F_wait(i), ~, ~] = coefTest(lme, H_conf_wait(i,:));
[pVal_go_wait(i), F_go_wait(i), DF1_2dec, DF2_2dec] = coefTest(lme, [H_conf_go(i,:);H_conf_wait(i,:)]);
end

F_crit_1dec = finv(0.95, DF1_1dec, DF2_1dec);
F_crit_2dec = finv(0.95, DF2_2dec, DF2_2dec);

chategory = {'go', 'wait', 'go : wait'}'; 
table_decision = table(chategory, [F_go(1); F_wait(1); F_go_wait(1)],[pVal_go(1);pVal_wait(1);pVal_go_wait(1)],...
    'VariableNames', {'Decision', 'F', 'pValue'});
table_metric = table(chategory, [F_go(2); F_wait(2); F_go_wait(2)],[pVal_go(2);pVal_wait(2);pVal_go_wait(2)],...
    'VariableNames', {'Decision', 'F', 'pValue'});
end 