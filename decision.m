function [Go_c, Wait_c, Go_c_pr, Wait_c_pr]=decision(Go, Wait, pos_c)
Go_c = []; 
Wait_c = []; 

for i= 1: length(pos_c)
    if ~isempty(find(pos_c(i) == Go, 1))
        Go_c(end+1)=pos_c(i);
    elseif ~isempty(find(pos_c(i) == Wait,1))
        Wait_c(end+1)=pos_c(i); 
    end 
end 

Go_c_pr = length(Go_c)/length(pos_c); 
Wait_c_pr = length(Wait_c)/length(pos_c); 

end 