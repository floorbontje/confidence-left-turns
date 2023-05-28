function [Go, Wait, change_of_mind,no_button,time, start_ego_inter, end_ego_inter,change_of_mind_kind]=All_inter(Log, Conf, individual)

if individual == 1 
    nbr=Log(1,1); 
    total_time = sum(Conf(:,10)); 
end 

%% Ego vehicle and subject vehicle pressent at intersection
% Boundary condition bot (subject vehicle) is present 
dis_ego = Log(:,8);
dis_bot = Log(:,9);
pos_9999 = find(dis_bot == 99999); 
dis_bot(pos_9999) = NaN;
dis_ego(pos_9999) = NaN;

% Only left turns  
Non_left=find(Log(:,6)<1);
dis_bot(Non_left) = NaN; 
dis_ego(Non_left) = NaN; 


% Moment that the bot enters the intersection and drives by the ego vehicle 
bot_inter = find(dis_bot < dis_ego); 
start_bot_inter = bot_inter([1;find(bot_inter(2:end)-bot_inter(1:end-1)>1)+1]);
end_bot_inter   = bot_inter([find(bot_inter(2:end)-bot_inter(1:end-1)>1);end]);

% Start_ego_inter: Ego vehicle has stopped (v<1m/s and dist <10m) and left the intersection
% End_ego_inter: End of the presence of the bot vehicle 
pos_value       = find(~isnan(dis_ego));
start_ego_inter = pos_value([1;find(pos_value(2:end)-pos_value(1:end-1)>1)+1]);
end_ego_inter   = pos_value([find(pos_value(2:end)-pos_value(1:end-1)>1);end]);



%% Executed behaviour 
% Wait and Go store the executed behaviour related to the position in the
% Conf. function. 
Go=[];
Wait=[];

for i = 1: length(start_ego_inter)
    for ii = 1: length(end_bot_inter)
            if find((end_ego_inter(i) > start_bot_inter(ii)) & ...
                    (start_ego_inter(i) < end_bot_inter(ii)) & ...
                    (end_ego_inter(i) > end_bot_inter(ii)))
                Wait(end+1)=i;
            end 
    end 
    
    if (isempty(Wait)) || (Wait(end)<i)
        Go(end+1)=i; 
    end 

end 

Wait = unique(Wait); 
Go = unique(Go); 
perf_Go   = zeros(length(Conf),1); 
perf_Wait =  zeros(length(Conf),1); 
perf_Go(Go)    = 1; 
perf_Wait(Wait)= 1; 

%% Changes of mind and no decision moment 
change_of_mind=[];
change_of_mind_kind=NaN(length(Conf),1);
changes_of_mind=NaN(length(Conf),1);
no_button=[];
no_button_vec=NaN(length(Conf),1); 

for i = 1: length(Conf(:,1))
    if perf_Go(i)==0 && Conf(i,7)==1 && perf_Wait(i)==1 && Conf(i,8)==0
        % 1: did not indicate a change of mind, present change is: go -> wait
        change_of_mind(end+1)=i; 
        change_of_mind_kind(i)=1; 
        
    elseif perf_Go(i)==0 && Conf(i,8)== 2
        % 2: indicated a change of mind, go -> wait, and executed this 
        change_of_mind(end+1)=i; 
        change_of_mind_kind(i)= 2; 
        
    elseif perf_Go(i)==1 && Conf(i,8)== 2
        % 3: indicated a change of mind, go -> wait, and did not execute this 
        change_of_mind(end+1)=i; 
        change_of_mind_kind(i)= 3; 
        
    elseif perf_Go(i)==1 && Conf(i,7)==0 && perf_Wait(i)==0 && Conf(i,8)==1 
        % -1: did not indicate a change of mind, present change is: wait -> go
        change_of_mind(end+1)=i; 
        change_of_mind_kind(i)= -1; 
        
    elseif perf_Go(i)==1 && Conf(i,7)== 2
        % -2: indicated a change of mind, wait -> go, and executed this 
        change_of_mind(end+1)=i; 
        change_of_mind_kind(i)= -2; 
        
    elseif perf_Go(i)==0 && Conf(i,7)== 2
        % -3: indicated a change of mind, wait -> go, and did not execute this 
        change_of_mind(end+1)=i; 
        change_of_mind_kind(i)= -3;    
     
    elseif Conf(i,7) == 0 && Conf(i,8) == 0 
        % did not press a button to indicate the decision 
        no_button(end+1)=i; 
        no_button_vec(i)=1; 
    end 
end

changes_of_mind(change_of_mind)=1; 


%% Time intervals   
% time interval 1: moment subject vehicle visible - moment indication of
%                  decision (delta_t1)
% time interval 2: moment indication decision - moment pause(delta_t2)
% time interval 3: moment pause - moment confidence rating 


end_decision  = NaN(size(start_ego_inter));
t_paused = NaN(size(start_ego_inter)); 

for i=1:length(start_ego_inter)
    ii=start_ego_inter(i);
    end_decision_pos= min([ii+find(Log(ii:end_ego_inter(i),33)>0,1),ii+find(Log(ii:end_ego_inter(i),32)>0,1)]); 
    if ~isempty(end_decision_pos)
        end_decision(i)=end_decision_pos; 
    else 
        end_decision(i) = ii; 
    end 
end 


delta_t1= Log(end_decision,7)-Log(start_ego_inter,7);
delta_t1(delta_t1==0)= NaN; 
delta_t2= Log(end_ego_inter,7) - Log(end_decision,7);
delta_t3= Conf(:,10)- Log(end_ego_inter,7); 

incorrect_button_use=find(isnan(delta_t1)); 
for i = 1: length(incorrect_button_use)
    if isempty(find(incorrect_button_use(i)==no_button,1))
        no_button(end+1)=(incorrect_button_use(i));
        no_button_vec(incorrect_button_use(i))= 0; 
    end
end

time.delta_1=delta_t1; 
time.delta_2=delta_t2;
time.delta_3=delta_t3;

if individual == 1 
    time.total=total_time; 
end 

%% Save 
    if individual == 1
        nbr = Conf(:,1); 
        tta_con = Conf(:,5); 
        d_con = Conf(:,6); 
        conf = Conf(:,4); 
        Go_ind= Conf(:,7); 
        Wait_ind=Conf(:,8);

        data = table(nbr,tta_con, d_con, conf, Go_ind, Wait_ind, ...
            perf_Go, perf_Wait,delta_t1, delta_t2,delta_t3, changes_of_mind, change_of_mind_kind, no_button_vec); 

        nbrtxt=num2str(Conf(1,1)); 
        fname= '.\Data analysis\indiv data analysis';
        writetable(data, fullfile(fname,[nbrtxt, '_all_int_indiv.txt']),'Delimiter','tab')
    end 
end 
