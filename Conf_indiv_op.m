function [pos_c, pos_log, Go_c, Wait_c]=Conf_indiv_op(Log, Conf, Go, Wait, change_of_mind, no_button, ii)

nbr_c=Log(1,1); 
nbrtxt=num2str(nbr_c); 
con = ['c1'; 'c2';'c3'; 'c4'];

tta=[5.5, 6.5, 5.5, 6.5]; 
dist=[70,70, 90, 90];


for i = 1:4 
    % Define positions of the conditions
    [pos, log] = positions(Conf,Log, tta(i), dist(i),change_of_mind, no_button); 
    pos_c.(con(i,:))=pos; 
    pos_log.(con(i,:))=log; 
    
    % Decisions and conditions 
    [Go_con, Wait_con, Go_pr, Wait_pr]=decision(Go, Wait, pos);

    Go_c.(con(i,:))= Go_con; 
    Wait_c.(con(i,:))= Wait_con;
    Go_c_pr.(con(i,:)) = Go_pr; 
    Wait_c_pr.(con(i,:)) = Wait_pr; 

end 

pos_tot=1:length(Conf(:,1)); 
[~, ~, Go_tot_pr, ~]=decision(Go, Wait, pos_tot);

% Confidence judgements: All decisions 
[mean_conf_c, mean_conf, std_conf_c, std_conf] = ...
    conf_anlysis(pos_c.c1, pos_c.c2, pos_c.c3, pos_c.c4, Conf);

% Confidence judgements: Go - decision 
[mean_conf_go_c, mean_conf_go, std_conf_go_c, std_conf_go] = ...
    conf_anlysis(Go_c.c1, Go_c.c2, Go_c.c3, Go_c.c4, Conf);


% Confidence judgements: Wait - decision 
[mean_conf_wait_c, mean_conf_wait, std_conf_wait_c, std_conf_wait] = ...
    conf_anlysis(Wait_c.c1, Wait_c.c2, Wait_c.c3, Wait_c.c4, Conf);


%% Save important data 

if ii == 1 
    nbr = [nbr_c]; 
    GoPr = [Go_tot_pr];
    GoC1Pr= [Go_c_pr.c1]; 
    GoC2Pr= [Go_c_pr.c2];
    GoC3Pr= [Go_c_pr.c3];
    GoC4Pr= [Go_c_pr.c4];

    conf_mean = [mean_conf]; 
    conf_mean_C = [mean_conf_c]; 
    conf_std = [std_conf];
    conf_std_C = [std_conf_c]; 

    conf_mean_go = [mean_conf_go]; 
    conf_mean_go_C = [mean_conf_go_c]; 
    conf_std_go = [std_conf_go];
    conf_std_go_C = [std_conf_go_c];

    conf_mean_wait = [mean_conf_wait]; 
    conf_mean_wait_C = [mean_conf_wait_c]; 
    conf_std_wait = [std_conf_wait];
    conf_std_wait_C = [std_conf_wait_c];
    save('Conf_data_indiv.mat','nbr')
else 
    load Conf_data_indiv.mat 
    nbr = [nbr; nbr_c];
    GoPr = [GoPr ; Go_tot_pr]; 
    GoC1Pr= [GoC1Pr; Go_c_pr.c1]; 
    GoC2Pr= [GoC2Pr; Go_c_pr.c2];
    GoC3Pr= [GoC3Pr; Go_c_pr.c3];
    GoC4Pr= [GoC4Pr; Go_c_pr.c4];

    conf_mean = [conf_mean; mean_conf]; 
    conf_mean_C = [conf_mean_C; mean_conf_c]; 
    conf_std = [conf_std; std_conf];
    conf_std_C = [conf_std_C; std_conf_c]; 

    conf_mean_go = [conf_mean_go; mean_conf_go]; 
    conf_mean_go_C = [conf_mean_go_C; mean_conf_go_c]; 
    conf_std_go = [conf_std_go; std_conf_go];
    conf_std_go_C = [conf_std_go_C; std_conf_go_c];

    conf_mean_wait = [conf_mean_wait; mean_conf_wait]; 
    conf_mean_wait_C = [conf_mean_wait_C; mean_conf_wait_c]; 
    conf_std_wait = [conf_std_wait; std_conf_wait];
    conf_std_wait_C = [conf_std_wait_C; std_conf_wait_c];
    save('Conf_data_indiv.mat','nbr','-append')
end 


save('Conf_data_indiv.mat','GoPr','GoC1Pr','GoC2Pr','GoC3Pr','GoC4Pr','-append')
save('Conf_data_indiv.mat','conf_mean','conf_mean_C','conf_std','conf_std_C','-append')
save('Conf_data_indiv.mat','conf_mean_go','conf_mean_go_C','conf_std_go','conf_std_go_C','-append')
save('Conf_data_indiv.mat','conf_mean_wait','conf_mean_wait_C','conf_std_wait','conf_std_wait_C','-append')

output = table(nbr_c,Go_tot_pr, Go_c_pr.c1,Go_c_pr.c2,Go_c_pr.c3,Go_c_pr.c4,...
mean_conf, mean_conf_c, std_conf, std_conf_c, mean_conf_go, mean_conf_go_c, std_conf_go, std_conf_go_c,...
mean_conf_wait, mean_conf_wait_c, std_conf_wait, std_conf_wait_c); 

fname='.\Data analysis\indiv data analysis';
writetable(output, fullfile(fname,[nbrtxt, '_conf_analysis.txt']),'Delimiter','tab')




end 