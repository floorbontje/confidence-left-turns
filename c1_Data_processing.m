% 1.Data processing 
% 7-9-2020
% By Floor Bontje 

clear all
close all 
clc


%%  Initialize 
Conf_all = [];  
Log_all  = [];
nbr_change_of_mind =[]; 
prs_changes_of_mind = []; 
nbr_no_button = []; 
prs_no_button = [];     
nbr_inters = []; 
prs_go =[]; prs_wait = [];

mean_delta_t1 = [];         mean_delta_t2 = [];         mean_delta_t3 = []; 
mean_delta_t1_Go = [];      mean_delta_t2_Go = [];      mean_delta_t3_Go = []; 
mean_delta_t1_Wait = [];    mean_delta_t2_Wait = [];    mean_delta_t3_Wait = []; 

CI_t1_go = [];      CI_t2_go = [];      CI_t3_go = [];
CI_t1_wait = [];    CI_t2_wait = [];    CI_t3_wait = [];
time_total = []; 


%% load data 
txt_dir = '.\Data\';

list=dir(txt_dir);
nbr_conf=[]; 
nbr_log =[];
for i = 1: length(list)
        if contains(list(i).name,'conf')
           nbr_conf(end+1) = i; 
        end 
        
        if contains(list(i).name,'log')
            nbr_log(end+1)=i; 
        end 
end 


for i= 1 :length(nbr_conf) 
    filename1=list(nbr_conf(i)).name;
    filename2=list(nbr_log(i)).name; 
    delimiterIn = '	';
    headerlinesIn = 1;

    Import_data_conf = importdata([txt_dir,filename1],delimiterIn,headerlinesIn);
    Conf=Import_data_conf.data; 
    Conf_all = [Conf_all;Conf];

    Import_data_log = importdata([txt_dir,filename2],delimiterIn,headerlinesIn);
    Log = Import_data_log.data;
    Log_all = [Log_all;Log];


    Conf_head = Import_data_conf.colheaders;
    Log_head  = Import_data_log.colheaders;

 
    %% All information for all intersections of one participant 
    close all 
    [Go, Wait, change_of_mind,no_button,time, ...
        ~, ~,change_of_mind_kind]=All_inter(Log, Conf, 1); 
    
    % Go:   indicates the positions of executed 'Go' in the Conf file
    % Wait: indicates the positions of executed 'Wait' in the Conf file  
    % change_of_mind: indicates the positions of changes of mind in Conf file 
    % no_button:    indicates the positions where no button was pressed in
        
    % time includes: 
    % delta_1: time between start decision moment and button press 
    % delta_2: time between button press and pause of simulation 
    % delta_3: time between pause of simulation and confidence judgement 
    % total  : time it took to drive all the routes during the
   
    nbr_change_of_mind(i)  = length(change_of_mind); 
    prs_changes_of_mind(i) = length(change_of_mind)/length(Conf(:,1))*100; 
    nbr_no_button(i)       = length(no_button); 
    prs_no_button(i)       = length(no_button)/length(Conf(:,1))*100;
    nbr_inters(i)          = length(Conf(:,1))-nbr_change_of_mind(end) - nbr_no_button(end); 
    
    prs_go(i) = length(Go)/length(Conf(:,1))*100; 
    prs_wait(i) = length(Wait)/length(Conf(:,1))*100; 
    
    excl_changes_no=[change_of_mind,no_button]; 
    pos=1:length(Conf(:,1));
    pos(excl_changes_no)=[]; 
    
    mean_delta_t1(i) = mean(time.delta_1(pos)); 
    mean_delta_t2(i) = mean(time.delta_2(pos)); 
    mean_delta_t3(i) = mean(time.delta_3(pos)); 
    time_total(i)= time.total/60; 
    
   
    
    %% Confidence judgements and Conditions 
    close all 
    [~, ~, Go_c, Wait_c] = Conf_indiv_op(Log, Conf, Go, Wait, change_of_mind, no_button, i);

    % Go_c: str with the positions go decision - condition in Conf matrix 
    % Wait_c: str with the positions of wait decision - condition in Conf matrix 
    con = ['c1'; 'c2';'c3'; 'c4'];
    for ii = 1:4 
        mean_delta_t1_Go(i,ii) = mean(time.delta_1(Go_c.(con(ii,:)))); 
        mean_delta_t2_Go(i,ii) = mean(time.delta_2(Go_c.(con(ii,:))));
        mean_delta_t3_Go(i,ii) = mean(time.delta_3(Go_c.(con(ii,:))));
    
        mean_delta_t1_Wait(i,ii) = mean(time.delta_1(Wait_c.(con(ii,:)))); 
        mean_delta_t2_Wait(i,ii) = mean(time.delta_2(Wait_c.(con(ii,:))));
        mean_delta_t3_Wait(i,ii) = mean(time.delta_3(Wait_c.(con(ii,:))));
        
        std_t1_go(i,ii)= std(time.delta_1(Go_c.(con(ii,:))));
        std_t2_go(i,ii)= std(time.delta_2(Go_c.(con(ii,:))));
        std_t3_go(i,ii)= std(time.delta_3(Go_c.(con(ii,:))));
        
        std_t1_wait(i,ii)= std(time.delta_1(Wait_c.(con(ii,:))));
        std_t2_wait(i,ii)= std(time.delta_2(Wait_c.(con(ii,:))));
        std_t3_wait(i,ii)= std(time.delta_3(Wait_c.(con(ii,:))));
    end 
        
   
    length_go = [length(Go_c.c1),length(Go_c.c2),length(Go_c.c3),length(Go_c.c4)]; 
  
    CI_t1_go(i,:) = 1.96*std_t1_go(i,:)./sqrt(length_go);
    CI_t2_go(i,:) = 1.96*std_t2_go(i,:)./sqrt(length_go);
    CI_t3_go(i,:) = 1.96*std_t3_go(i,:)./sqrt(length_go);
  
    length_wait = [length(Wait_c.c1),length(Wait_c.c2),length(Wait_c.c3),length(Wait_c.c4)]; 
       
    CI_t1_wait(i,:) = 1.96*std_t1_wait(i,:)./sqrt(length_wait);
    CI_t2_wait(i,:) = 1.96*std_t2_wait(i,:)./sqrt(length_wait);
    CI_t3_wait(i,:) = 1.96*std_t3_wait(i,:)./sqrt(length_wait);
        
    
end 

%% Save confindence matrices and log matrices of all participants
close all 
save('Conf_all.mat','Conf_all')
save('Log_all.mat','Log_all')


%% Behaviour information of all participants 
nbr= unique(Conf_all(:,1))';

output1 = table(nbr',prs_go',prs_wait',prs_changes_of_mind',prs_no_button',...
    mean_delta_t1',mean_delta_t2',mean_delta_t3',time_total', nbr_inters', 'VariableNames',...
    {'ID','prs_go','prs_wait','prs_change_of_mind','prs_no_button',...
    'mean_delta_t1','mean_delta_t2','mean_delta_t3','total_time','nbr_inters'});


%% confidence information of all participants 
load Conf_data_indiv.mat 
output2 = table(nbr,GoPr, GoC1Pr,GoC2Pr,GoC3Pr,GoC4Pr,...
conf_mean, conf_mean_C, conf_std, conf_std_C, conf_mean_go, conf_mean_go_C, ...
conf_std_go, conf_std_go_C, conf_mean_wait, conf_mean_wait_C, ...
conf_std_wait, conf_std_wait_C, mean_delta_t1_Go,mean_delta_t2_Go,mean_delta_t3_Go, CI_t1_go, CI_t2_go, CI_t3_go,...
mean_delta_t1_Wait,mean_delta_t2_Wait,mean_delta_t3_Wait, CI_t1_wait, CI_t2_wait, CI_t3_wait); 


%% Save 
fname='.\Data analysis\general data analysis';
writetable(output1, fullfile(fname,'behaviour_analysis.txt'),'Delimiter','tab')
writetable(output2, fullfile(fname,'conf_analysis.txt'),'Delimiter','tab')










