function [conf_c_mean, conf_mean, conf_c_std, conf_std] = ...
    conf_anlysis(pos1_c, pos2_c, pos3_c, pos4_c, Conf)
close all 

conf_c1 = Conf(pos1_c,4); 
conf_c2 = Conf(pos2_c,4); 
conf_c3 = Conf(pos3_c,4); 
conf_c4 = Conf(pos4_c,4); 

mean_conf_c1=mean(conf_c1,'omitnan');
mean_conf_c2=mean(conf_c2,'omitnan');
mean_conf_c3=mean(conf_c3,'omitnan');
mean_conf_c4=mean(conf_c4,'omitnan');

conf_c_mean =[mean_conf_c1, mean_conf_c2, mean_conf_c3, mean_conf_c4]; 
conf_mean= mean([conf_c1;conf_c2;conf_c3;conf_c4]); 

std_conf_c1=std(conf_c1,'omitnan');
std_conf_c2=std(conf_c2,'omitnan');
std_conf_c3=std(conf_c3,'omitnan');
std_conf_c4=std(conf_c4,'omitnan');

conf_c_std =[std_conf_c1, std_conf_c2, std_conf_c3, std_conf_c4];
conf_std = std([conf_c1;conf_c2;conf_c3;conf_c4]); 


end 