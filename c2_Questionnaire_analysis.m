% Questionnaire analysis 
% By Floor Bontje 

% No conclustions about correlation can be drawn from this analysis due to
% small number of participants (the power of the analysis). However, this
% script gives an indication of the diversity of the different participants
% who participated in the research. 

clear all 
clc 

Questionnaire_output = readtable('Data\Questionnaire_output.txt','Format', '%f%f%C%f%C%f%f%C%C%C%f%C'); 
load Conf_data_indiv.mat


%% Input: information obtained from the questionnaire 
[~,ID] =        sort(Questionnaire_output.ID); 

age =           Questionnaire_output.Age(ID); 
general_conf =  Questionnaire_output.Confidence_general(ID); 
driving_area =  Questionnaire_output.Driving_area(ID); 
frequency_driving = Questionnaire_output.Frequency_driving(ID); 
frequency_driving = categorical(frequency_driving,{'Almost never (< 1 day a week)',...
    'Sometimes (1 day a week)',' Regularly (2-3 days a week)','Often (4-5 days a week)','Very often (> 6 days a week)'}); 
gaming_experience = Questionnaire_output.Gaming_experience(ID); 
gender =        Questionnaire_output.Gender(ID); 
risk_taking =   Questionnaire_output.Risk_taking_driving(ID); 
driving_license =   Questionnaire_output.Time_driving_licence(ID); 
year_km =       Questionnaire_output.year_km(ID); 


mean_age = mean(age); 
std_age = std(age); 

%% Demographic information 
figure1 = figure; 
subplot(2,2,1)
histogram(age, 'Normalization','probability')
title('Age') 
xlabel('Age [years]')
ylabel('Proportion') 

subplot(2,2,2)
histogram(gender, 'Normalization', 'probability')
title('Gender')
ylabel('Proportion')

subplot(2,2,3)
histogram(frequency_driving, 'Normalization', 'probability') 
title('Frequency of driving')
ylabel('Proportion')

subplot(2,2,4)
histogram(driving_license, 'Normalization', 'probability')
title("Possession of driver's license")
xlabel('years')
ylabel('Proportion')
sgtitle('General information')


figure2=figure; 
subplot(1,2,1)
histogram(general_conf, 'Normalization', 'probability')
title('Global confidence') 
xlabel('Confidence on a scale [1,5]')
ylabel('Proportion') 

subplot(1,2,2)
histogram(risk_taking, 'Normalization', 'probability')
title('Risk taking driving behaviour') 
xlabel('Amount of risk on a scale [1,5]')
ylabel('Proportion')
sgtitle('Self-reflecion on driving style')


%% Poweranalysis
% mean_M=mean(conf_mean_go(gender == 'Male',1)); 
% std_M=std(conf_mean_go(gender == 'Male',1));
% mean_F=mean(conf_mean_go(gender == 'Female',1), 'omitnan');
% std_F=std(conf_mean_go(gender == 'Female',1), 'omitnan');
% 
% n_gen=sampsizepwr('t',[mean_M, std_M],[mean_F]) ;


%% Save results 
fname = '.\Data analysis\figures\Questionnaire analysis';
saveas(figure1, fullfile(fname, 'Info participants.jpg'))
saveas(figure1, fullfile(fname, 'Global Confidence.jpg'))  




