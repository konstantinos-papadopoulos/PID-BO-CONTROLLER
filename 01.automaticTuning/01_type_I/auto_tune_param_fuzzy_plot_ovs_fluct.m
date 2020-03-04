% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_fuzzy_plot_ovs_fluct.m																			   																		  	
%  Project:     Automatic tuning of the parameters for PI,PID controllers
%  
%  Purpose:     main script for the automatic tuning based on the three methods																		   																		
%  Author :     kostas g. papadopoulos																	   																		
% 																										   																		
%  History:     Date: 07.07.2008  date last modified
% 																										  																		
%  Contact:     kostas g. papadopoulos,    nikos mitrakis       leonidas droukas
%               kpapadop@eng.auth.gr  ,    nmitr@auth.gr        leon_drouk@yahoo.gr
% 																										  																		
%  Place:	    Aristotle University of Thessaloniki, Thessaloniki, Greece							   																		
% 
% ************************************************************************************************************
function  [] = auto_tune_param_fuzzy_plot_ovs_fluct(auto_tune_param_fuzzy_plot_ovs_fluctStruct)

% storing back the variables from main
% -------------------------------------------------------------------------
figureIndexLocal = auto_tune_param_fuzzy_plot_ovs_fluctStruct.figureIndex                                 ;

log_ovs_while_fuzzy_pi_tuning = auto_tune_param_fuzzy_plot_ovs_fluctStruct.log_ovs_while_fuzzy_pi_tuning  ;
log_ovs_while_fuzzy_pid_tuning = auto_tune_param_fuzzy_plot_ovs_fluctStruct.log_ovs_while_fuzzy_pid_tuning;

step_Conv_fuzzyTuning_complexZeros_PI_Controller = auto_tune_param_fuzzy_plot_ovs_fluctStruct.step_Conv_fuzzyTuning_complexZeros_PI_Controller;
step_Conv_fuzzyTuning_complexZeros_PID_Controller = auto_tune_param_fuzzy_plot_ovs_fluctStruct.step_Conv_fuzzyTuning_complexZeros_PID_Controller;

ovs_PI_fuzzy_tune = auto_tune_param_fuzzy_plot_ovs_fluctStruct.ovs_PI_fuzzy_tune               ;
ovs_PID_fuzzy_tune = auto_tune_param_fuzzy_plot_ovs_fluctStruct.ovs_PID_fuzzy_tune             ;

figure(figureIndexLocal)
% string parameters for the LEGEND
% -------------------------------------------------------------------------
lowerlimitOvrst_PI_autotunestr  = num2str(log_ovs_while_fuzzy_pi_tuning(1))               ; % log_ovs_while_auto_tuning(1) = lowerlimitOvrst;
upperlimitOvrst_PI_autotunestr  = num2str(log_ovs_while_fuzzy_pi_tuning(2))               ; % log_ovs_while_auto_tuning(2) = upperlimitOvrst;
referencelimitOvrst_PI_autotunestr  = num2str(log_ovs_while_fuzzy_pi_tuning(4))           ; % log_ovs_while_auto_tuning(4) = desired...output from the pi fuzzy estimator;
optimalPIovs = num2str(log_ovs_while_fuzzy_pi_tuning(5))                                  ;

lowerlimitOvrst_PID_autotunestr  = num2str(log_ovs_while_fuzzy_pid_tuning(1))             ; % log_ovs_while_auto_tuning(1) = lowerlimitOvrst;
upperlimitOvrst_PID_autotunestr  = num2str(log_ovs_while_fuzzy_pid_tuning(2))             ; % log_ovs_while_auto_tuning(2) = upperlimitOvrst;
referencelimitOvrst_PID_autotunestr  = num2str(log_ovs_while_fuzzy_pid_tuning(4))         ; % log_ovs_while_auto_tuning(4) = desired...output from the pi fuzzy estimator;
optimalPIDovs = num2str(log_ovs_while_fuzzy_pid_tuning(5))                                ;

lowerlimitLegend_PI = strcat('lower Limit OVS is:',lowerlimitOvrst_PI_autotunestr)        ;
upperlimitLegend_PI = strcat('upper Limit OVS is:',upperlimitOvrst_PI_autotunestr)        ;
referencelimitLegend_PI = strcat('reference OVS is:',referencelimitOvrst_PI_autotunestr)  ;

lowerlimitLegend_PID = strcat('lower Limit OVS is:',lowerlimitOvrst_PID_autotunestr)      ;
upperlimitLegend_PID = strcat('upper Limit OVS is:',upperlimitOvrst_PID_autotunestr)      ;
referencelimitLegend_PID = strcat('reference OVS is:',referencelimitOvrst_PID_autotunestr);

% LINEs for upper, lower and reference - I Controller
% -------------------------------------------------------------------------
lowerLimit_PI_auto_tune = log_ovs_while_fuzzy_pi_tuning(1)*ones(step_Conv_fuzzyTuning_complexZeros_PI_Controller,1)        ;
upperLimit_PI_auto_tune = log_ovs_while_fuzzy_pi_tuning(2)*ones(step_Conv_fuzzyTuning_complexZeros_PI_Controller,1)        ;
referenceLimit_PI_auto_tune = log_ovs_while_fuzzy_pi_tuning(4)*ones(step_Conv_fuzzyTuning_complexZeros_PI_Controller,1)    ;
optimalReference_PI_auto_tune = log_ovs_while_fuzzy_pi_tuning(5)*ones(step_Conv_fuzzyTuning_complexZeros_PI_Controller,1)    ;

% LINEs for upper, lower and reference - PID Controller
% -------------------------------------------------------------------------
lowerLimit_PID_auto_tune = log_ovs_while_fuzzy_pid_tuning(1)*ones(step_Conv_fuzzyTuning_complexZeros_PID_Controller,1)    ;
upperLimit_PID_auto_tune = log_ovs_while_fuzzy_pid_tuning(2)*ones(step_Conv_fuzzyTuning_complexZeros_PID_Controller,1)    ;
referenceLimit_PID_auto_tune = log_ovs_while_fuzzy_pid_tuning(4)*ones(step_Conv_fuzzyTuning_complexZeros_PID_Controller,1);
optimalReference_PID_auto_tune = log_ovs_while_fuzzy_pid_tuning(5)*ones(step_Conv_fuzzyTuning_complexZeros_PID_Controller,1)    ;

% Final Values after tuning
% -------------------------------------------------------------------------
ovs_PI_fuzzy_auto_tuneFinal = ovs_PI_fuzzy_tune(length(ovs_PI_fuzzy_tune))      ;
ovs_PID_fuzzy_auto_tuneFinal = ovs_PID_fuzzy_tune(length(ovs_PID_fuzzy_tune))   ;

ovs_PI_auto_tuneFinal_str = num2str(ovs_PI_fuzzy_auto_tuneFinal)      ; 
ovs_PID_auto_tuneFinal_str = num2str(ovs_PID_fuzzy_auto_tuneFinal)    ; 



subplot(2,1,1) %PI Autotuned Controller
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
line([1:length(lowerLimit_PI_auto_tune)],...
     lowerLimit_PI_auto_tune,'Color','k',...
     'LineStyle','--','LineWidth',0.5)
hold on

line([1:length(upperLimit_PI_auto_tune)],...
     upperLimit_PI_auto_tune,'Color','k',...
     'LineStyle','--','LineWidth',0.5)
hold on
line([1:length(referenceLimit_PI_auto_tune)],...
    referenceLimit_PI_auto_tune,'Color','r',...
    'LineStyle','--','LineWidth',0.5)
hold on

line([1:length(optimalReference_PI_auto_tune)],...
    optimalReference_PI_auto_tune,'Color','g',...
    'LineStyle','--','LineWidth',2)
hold on

plot(ovs_PI_fuzzy_tune,'-b*')
grid on
finalovs_str1 = strcat('Overshoot Final Value:',ovs_PI_auto_tuneFinal_str);
finalovs_str2 = strcat(' - Overshoot Optimal:',optimalPIovs);
finalovs_str3 = strcat(' - tuning converged in:',num2str(step_Conv_fuzzyTuning_complexZeros_PI_Controller),' steps');
title({'PI Controller - Overshoot fluctuation while tuning';strcat(finalovs_str1,finalovs_str2,finalovs_str3)})
legend(lowerlimitLegend_PI,upperlimitLegend_PI,referencelimitLegend_PI,0)
% -------------------------------------------------------------------------

subplot(2,1,2) % PID - Cotnroller - Real Zeros
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
line([1:length(lowerLimit_PID_auto_tune)],...
     lowerLimit_PID_auto_tune,'Color','k',...
     'LineStyle','--','LineWidth',0.5)
hold on

line([1:length(upperLimit_PID_auto_tune)],...
     upperLimit_PID_auto_tune,'Color','k',...
     'LineStyle','--','LineWidth',0.5)
hold on
line([1:length(referenceLimit_PID_auto_tune)],...
    referenceLimit_PID_auto_tune,'Color','r',...
    'LineStyle','--','LineWidth',0.5)
hold on

line([1:length(optimalReference_PID_auto_tune)],...
    optimalReference_PID_auto_tune,'Color','g',...
    'LineStyle','--','LineWidth',2)
hold on

plot(ovs_PID_fuzzy_tune,'-b*')
grid on
finalovs_str1 = strcat('Overshoot Final Value:',ovs_PID_auto_tuneFinal_str);
finalovs_str2 = strcat(' - Overshoot Optimal:',optimalPIDovs);
finalovs_str3 = strcat(' - tuning converged in:',num2str(step_Conv_fuzzyTuning_complexZeros_PID_Controller),' steps');

title({'PID Controller - Overshoot fluctuation while tuning';...
       strcat(finalovs_str1,finalovs_str2,finalovs_str3)})
legend(lowerlimitLegend_PID,upperlimitLegend_PID,referencelimitLegend_PID,0)
% -------------------------------------------------------------------------

