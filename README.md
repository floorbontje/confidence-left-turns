# confidence-left-turns
Are you sure? Modelling the Confidence of a Driver in Left-Turn Gap Acceptance Decisions

In this study, several different MATLAB scripts were used to process the data from the experiment (Data_processing.m, Questionnaire_analysis.m), to perform analysis on the decision behaviour and confidence judgements (Data_analysis_confidence.m) and to create left-turn gap acceptance/rejection decision models (Training_Decision.m) and confidence models (Modelling_Confidence.m). In order for the code to function properly, all MATLAB scripts should be stored in the same folder, and in addition the following subfolders should be present in this folder: 
-	Data
-	Data analysis  
    -	figures
        - Action Dynamics
        - analysis
        - data distribution
        - Questionnaire analysis
    - general data analysis
	    - action dynamics 
      -	lme models
 
    - indiv data analysis
-	Modelling
    -	Confidence
        - Model1
        - Model2
        -	Model3
        -	Model4
     - data model feeting 
     - decision 
        -	figures
        -	Parameters
     -	Parameters Initial
  
Moreover, for the driving simulator experiment Carla and Unreal Engine4 was used in combination with Python. See folder “Data collection code > experiment”. In addition, the initial parameters of the drift-diffusion decision model were obtained with the use of the code python (“03_fit_model.py”)  presented in the study “Should I stay, or should I go? Evidence accumulation drives decision making in human drivers” by Zgonnikov et al. 2020 (see folder “Data collection code > model feeting > data_and_model_analysis”). 

1.	c1_Data_processing.m 
This script is used to load and read the data obtained from the experiment which is stored in the folder “Data”. Moreover, the script creates “Conf_all.mat” and “Log_all.mat”, the two data files used for further analysis. 

2. c2_Questionnaire_analysis.m
This script gives an indication about the diversity of participants who took part in the experiment. 
 
3. c3_Data_analysis_confidence.m
The script “Data_analysis_confidence.m” is used to perform several analyses on the obtained data of all participants on the left-turn gap acceptance decision behaviour and confidence judgements. Furthermore, this code is also used for the visualisation of the obtained data on the decision behaviour, response time, initial throttle operation moment and confidence judgements. 

The script is used to create and to visualise the linear mixed effects models (lme) and associated post-hoc hypothesis analysis presented in the paper. The following lme models are created: 
  -	Decision behaviour ~ decision outcome, time-to-arrival, and distance gap.  
  -	Response time ~ decision outcome, time-to-arrival, and distance gap. 
  -	Confidence ~ response time, decision outcome, time-to-arrival, and distance gap.
  -	Confidence ~ initial throttle operation moment, decision outcome, time-to-arrival, and distance gap. 

The additional function “function_lme_results.m” is used for the visualisation of the lme models and analysis of the random effects due to individual differences between participants. 

Moreover, it computes correlations between: 
  -	Confidence and the response time. 
  -	Confidence and the initial throttle operation moment.
  -	Response time and the initial throttle operation moment. 

The relation between confidence and the action dynamics measures, the velocity profile and distance to the centre of the intersection, is investigated. For this analysis the additional function “analysis_action_dynamics.m” is used. This function performs a linear regression analysis between confidence and the metrics of the action dynamics and creates a visualisation of the effect of confidence on the action dynamic measure over time. The metrics of the action dynamic measures that are investigated constitute the RMSD, individual MD, group MD, minimum and maximum value. 
Figures of the group and individual mean of the decision behaviour, the response times and the confidence judgements are also made in this script. Also, the experiment data is exported (in the file “data_output_RT.csv”) so that it can be used for the parameter tuning of the decision and confidence models. 

4. c4_Training_Decision.m
This script is used to create and train the two decision models presented in the paper, the drift-diffusion and the race model. These models are based on a dynamic drift-diffusion left-turn gap acceptance decision model described by Zgonnikov et al. 2020. This script searches for the optimal model parameters, using the “fmincon” function of MATLAB. The function “function_Decision_Model.m” describes the decision model based on the defined model parameters. 

5. c5_Modelling_Confidence.m
This script is used to create different potential confidence models which are built on the basis of a decision model. First, the script evaluates the different left-turn gap acceptance decision models, after which two confidence models that do not allow for additional evidence accumulation are created. For the drift-diffusion based confidence model, the “function_Conf_DDM.m” function is used and for the race model-based confidence model the “function_Conf_Race_model.m” is used. 
The script also investigates the effect of additional evidence accumulation on the performance of the  confidence models. For the additional evidence accumulation, we defined a maximum duration of 2.5 seconds.
The mean performance of the different confidence models over 10 models is computed and stored as “RMSD_models_mean.csv”.  
