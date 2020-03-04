% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_plot_ovs_fluct.m																			   																		  	
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
function  [] = auto_tune_param_plot_ovs_fluct(auto_tune_param_plot_ovs_fluctStructLoc)

% storing back the variables from main
% -------------------------------------------------------------------------
figureIndexLocal = auto_tune_param_plot_ovs_fluctStructLoc.figureIndex                            ;
log_ovs_while_auto_tuning = auto_tune_param_plot_ovs_fluctStructLoc.log_ovs_while_auto_tuning;
step_Conv_autoTuning_realZeros_I_Controller = auto_tune_param_plot_ovs_fluctStructLoc.step_Conv_autoTuning_realZeros_I_Controller;
step_Conv_autoTuning_realZeros_PI_Controller = auto_tune_param_plot_ovs_fluctStructLoc.step_Conv_autoTuning_realZeros_PI_Controller;
step_Conv_autoTuning_realZeros_PID_Controller = auto_tune_param_plot_ovs_fluctStructLoc.step_Conv_autoTuning_realZeros_PID_Controller;
step_Conv_autoTuning_complexZeros_PID_Controller = auto_tune_param_plot_ovs_fluctStructLoc.step_Conv_autoTuning_complexZeros_PID_Controller;




ovs_I_auto_tune = auto_tune_param_plot_ovs_fluctStructLoc.ovs_I_auto_tune                 ;
ovs_PI_auto_tune = auto_tune_param_plot_ovs_fluctStructLoc.ovs_PI_auto_tune               ;
ovs_PID_auto_tune = auto_tune_param_plot_ovs_fluctStructLoc.ovs_PID_auto_tune             ;
ovs_PID_auto_tune_ypred = auto_tune_param_plot_ovs_fluctStructLoc.ovs_PID_auto_tune_ypred ;

figure(figureIndexLocal)
% string parameters for the LEGEND
% -------------------------------------------------------------------------
lowerlimitOvrst_I_autotunestr  = num2str(log_ovs_while_auto_tuning(1))              ; % log_ovs_while_auto_tuning(1) = lowerlimitOvrst;
upperlimitOvrst_I_autotunestr  = num2str(log_ovs_while_auto_tuning(2))              ; % log_ovs_while_auto_tuning(2) = upperlimitOvrst;
referencelimitOvrst_I_autotunestr  = num2str(log_ovs_while_auto_tuning(3))          ; % log_ovs_while_auto_tuning(3) = referenceOvrst ;

lowerlimitLegend = strcat('lower Limit OVS is:',lowerlimitOvrst_I_autotunestr)      ;
upperlimitLegend = strcat('upper Limit OVS is:',upperlimitOvrst_I_autotunestr)      ;
referencelimitLegend = strcat('reference OVS is:',referencelimitOvrst_I_autotunestr);

% LINEs for upper, lower and reference - I Controller
% -------------------------------------------------------------------------
lowerLimit_I_auto_tune = log_ovs_while_auto_tuning(1)*ones(step_Conv_autoTuning_realZeros_I_Controller,1)        ;
upperLimit_I_auto_tune = log_ovs_while_auto_tuning(2)*ones(step_Conv_autoTuning_realZeros_I_Controller,1)        ;
referenceLimit_I_auto_tune = log_ovs_while_auto_tuning(3)*ones(step_Conv_autoTuning_realZeros_I_Controller,1)    ;


% LINEs for upper, lower and reference - PI Controller
% -------------------------------------------------------------------------
lowerLimit_PI_auto_tune = log_ovs_while_auto_tuning(1)*ones(step_Conv_autoTuning_realZeros_PI_Controller,1)      ;
upperLimit_PI_auto_tune = log_ovs_while_auto_tuning(2)*ones(step_Conv_autoTuning_realZeros_PI_Controller,1)      ;
referenceLimit_PI_auto_tune = log_ovs_while_auto_tuning(3)*ones(step_Conv_autoTuning_realZeros_PI_Controller,1)  ;

% LINEs for upper, lower and reference - PID Controller
% -------------------------------------------------------------------------
lowerLimit_PID_auto_tune = log_ovs_while_auto_tuning(1)*ones(step_Conv_autoTuning_realZeros_PID_Controller,1)    ;
upperLimit_PID_auto_tune = log_ovs_while_auto_tuning(2)*ones(step_Conv_autoTuning_realZeros_PID_Controller,1)    ;
referenceLimit_PID_auto_tune = log_ovs_while_auto_tuning(3)*ones(step_Conv_autoTuning_realZeros_PID_Controller,1);

% LINEs for upper, lower and reference - PID Controller
% -------------------------------------------------------------------------
lowerLimit_PID_auto_tune_complex = log_ovs_while_auto_tuning(1)*ones(step_Conv_autoTuning_complexZeros_PID_Controller,1);
upperLimit_PID_auto_tune_complex = log_ovs_while_auto_tuning(2)*ones(step_Conv_autoTuning_complexZeros_PID_Controller,1);
referenceLimit_PID_auto_tune_complex = log_ovs_while_auto_tuning(3)*ones(step_Conv_autoTuning_complexZeros_PID_Controller,1);


% Final Values after tuning
% -------------------------------------------------------------------------
ovs_I_auto_tuneFinal = ovs_I_auto_tune(length(ovs_I_auto_tune))                        ;
ovs_PI_auto_tuneFinal = ovs_PI_auto_tune(length(ovs_PI_auto_tune))                     ;
ovs_PID_auto_tuneFinal = ovs_PID_auto_tune(length(ovs_PID_auto_tune))                  ;
ovs_PID_auto_tune_ypredFinal = ovs_PID_auto_tune_ypred(length(ovs_PID_auto_tune_ypred));

ovs_I_auto_tuneFinal_str = num2str(ovs_I_auto_tuneFinal)        ; 
ovs_PI_auto_tuneFinal_str = num2str(ovs_PI_auto_tuneFinal)       ; 
ovs_PID_auto_tuneFinal_str = num2str(ovs_PID_auto_tuneFinal)      ; 
ovs_PID_auto_tune_ypredFinal_str = num2str(ovs_PID_auto_tune_ypredFinal); 


subplot(2,2,1)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
line([1:length(lowerLimit_I_auto_tune)],...
     lowerLimit_I_auto_tune,'Color','k',...
     'LineStyle','--','LineWidth',0.5)
hold on

line([1:length(upperLimit_I_auto_tune)],...
     upperLimit_I_auto_tune,'Color','k',...
     'LineStyle','--','LineWidth',0.5)
hold on
line([1:length(referenceLimit_I_auto_tune)],...
    referenceLimit_I_auto_tune,'Color','r',...
    'LineStyle','--','LineWidth',0.5)
hold on

plot(ovs_I_auto_tune,'-b*')
grid on
finalovs_str1 = strcat('Overshoot Final Value:',ovs_I_auto_tuneFinal_str);
finalovs_str2 = strcat(' - tuning converged in: ',num2str(step_Conv_autoTuning_realZeros_I_Controller),' steps');

title({'I Controller - Overshoot fluctuation while tuning';strcat(finalovs_str1,finalovs_str2)})
legend(lowerlimitLegend,upperlimitLegend,referencelimitLegend,0)
% -------------------------------------------------------------------------

subplot(2,2,2) %PI Autotuned Controller
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

plot(ovs_PI_auto_tune,'-b*')
grid on
finalovs_str1 = strcat('Overshoot Final Value:',ovs_PI_auto_tuneFinal_str);
finalovs_str2 = strcat(' - tuning converged in: ',num2str(step_Conv_autoTuning_realZeros_PI_Controller),' steps');
title({'PI Controller - Overshoot fluctuation while tuning';strcat(finalovs_str1,finalovs_str2)})
legend(lowerlimitLegend,upperlimitLegend,referencelimitLegend,0)
% -------------------------------------------------------------------------



subplot(2,2,3) % PID - Cotnroller - Real Zeros
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

plot(ovs_PID_auto_tune,'-b*')
grid on
finalovs_str1 = strcat('Overshoot Final Value:',ovs_PID_auto_tuneFinal_str);
finalovs_str2 = strcat(' - tuning converged in: ',num2str(step_Conv_autoTuning_realZeros_PID_Controller),' steps');
title({'PID Controller - Overshoot fluctuation while tuning';...
       'real zeros controller';strcat(finalovs_str1,finalovs_str2)})
legend(lowerlimitLegend,upperlimitLegend,referencelimitLegend,0)
% -------------------------------------------------------------------------



subplot(2,2,4) % PID Controller - Complex Zeros
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
line([1:length(lowerLimit_PID_auto_tune_complex)],...
     lowerLimit_PID_auto_tune_complex,'Color','k',...
     'LineStyle','--','LineWidth',0.5)
hold on

line([1:length(upperLimit_PID_auto_tune_complex)],...
     upperLimit_PID_auto_tune_complex,'Color','k',...
     'LineStyle','--','LineWidth',0.5)
hold on
line([1:length(referenceLimit_PID_auto_tune_complex)],...
    referenceLimit_PID_auto_tune_complex,'Color','r',...
    'LineStyle','--','LineWidth',0.5)
hold on

plot(ovs_PID_auto_tune_ypred,'-b*')
grid on
finalovs_str1 = strcat('Overshoot Final Value:',ovs_PID_auto_tune_ypredFinal_str);
finalovs_str2 = strcat(' - tuning converged in: ',num2str(step_Conv_autoTuning_complexZeros_PID_Controller),' steps');
title({'PID Controller - Overshoot fluctuation while tuning';
       'complex zeros controller';strcat(finalovs_str1,finalovs_str2)})
legend(lowerlimitLegend,upperlimitLegend,referencelimitLegend,0)
% -------------------------------------------------------------------------





