function [pos_c, pos_log, Go_c, Wait_c, Go_c_pr, Wait_c_pr]=Conf_global_op(Log, Conf, Go, Wait, change_of_mind, no_button)

% Conditions
con = ['c1'; 'c2';'c3'; 'c4'];
tta=[5.5, 6.5, 5.5, 6.5]; 
dist=[70,70, 90, 90];


for i = 1:4 
    % Define positions of the conditions
    [pos, log] = positions(Conf,Log, tta(i), dist(i),change_of_mind, no_button); 
    pos_c.(con(i,:))= pos; 
    pos_log.(con(i,:))=log; 
    
    % Decisions and conditions 
    [Go_con, Wait_con, Go_pr, Wait_pr]=decision(Go, Wait, pos);

    Go_c.(con(i,:))= Go_con; 
    Wait_c.(con(i,:))= Wait_con;
    Go_c_pr.(con(i,:)) = Go_pr; 
    Wait_c_pr.(con(i,:)) = Wait_pr; 

end 


% Save outputs 
save('.\Data analysis\general data analysis\Decision_behaviour_all.mat', 'Go_c_pr', 'Wait_c_pr')

end 