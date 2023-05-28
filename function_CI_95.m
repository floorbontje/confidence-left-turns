function [CI] = function_CI_95(data, pos, nbr)

% Function to calculate 95% confidence interval (CI) 
% Data: input data
% Pos: positions of data of a condition 
% nbr: nbr of conditions 

for i = 1: nbr
    CI(i) = 1.96*std(data(pos(i,:)))/sqrt(length(data(pos(i,:))));
end 