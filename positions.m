function [pos_c, pos_log] = positions(Conf,Log, tta, dist,change_of_mind, no_button)

    % Position of the defined condition in the Conf matrix  
    pos_c = find(Conf(:,5) == tta & Conf(:,6) == dist);
    
    % All changes of mind are excluded from the analysis
    pos_excl=[change_of_mind,no_button]; 
    
    for i = 1:length(pos_c)
        if ~isempty(find(pos_c(i)==pos_excl,1))
            pos_c(i)=NaN; 
        end 
    end
    
    pos_c(isnan(pos_c))=[];
    
    % Route & intersection combination the condition 
    RI = Conf(pos_c,[2,3]); 
   
    % Positions of left turns with the defined condition in the Log matrix 
    pos1_log=[]; 

    for i = 1:length(RI)
    new = find(Log(:,2) == RI(i,1) & Log(:,3) == RI(i,2) & Log(:,6) == 1); 
    pos_log = [pos1_log;new];
    end 
    
 
end 
