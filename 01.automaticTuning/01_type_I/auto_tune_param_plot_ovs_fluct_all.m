% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_plot_ovs_fluct_all.m																			   																		  	
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
function  [] = auto_tune_param_plot_ovs_fluct_all(auto_tune_param_all_plot_ovs_fluctStruct)

% restoring the matrices
% *************************************************************************
ovs_Final_auto_tune = auto_tune_param_all_plot_ovs_fluctStruct.ovs_Final_auto_tune                  ;
ovs_Final_auto_tune_margaris = auto_tune_param_all_plot_ovs_fluctStruct.ovs_Final_auto_tune_margaris;
ovs_Final_fuzzy_tune = auto_tune_param_all_plot_ovs_fluctStruct.ovs_Final_fuzzy_tune                ;

length_ovs_Final_auto_tune = length(ovs_Final_auto_tune);
length_ovs_Final_auto_tune_margaris = length(ovs_Final_auto_tune_margaris);
length_ovs_Final_fuzzy_tune = length(ovs_Final_fuzzy_tune);
lengthMatr = [length_ovs_Final_auto_tune length_ovs_Final_auto_tune_margaris length_ovs_Final_fuzzy_tune];
x_min_fig = 0;
x_max_fig = max(lengthMatr);

y_max_fig_ovs_Final_auto_tune = max(ovs_Final_auto_tune);
y_max_fig_ovs_Final_auto_tune_margaris = max(ovs_Final_auto_tune_margaris);
y_max_fig_ovs_Final_fuzzy_tune = max(ovs_Final_fuzzy_tune);
y_max_fig = max([y_max_fig_ovs_Final_auto_tune y_max_fig_ovs_Final_auto_tune_margaris y_max_fig_ovs_Final_fuzzy_tune]);

y_min_fig_ovs_Final_auto_tune = min(ovs_Final_auto_tune);
y_min_fig_ovs_Final_auto_tune_margaris = min(ovs_Final_auto_tune_margaris);
y_min_fig_ovs_Final_fuzzy_tune = min(ovs_Final_fuzzy_tune);
y_min_fig = min([y_min_fig_ovs_Final_auto_tune y_min_fig_ovs_Final_auto_tune_margaris y_min_fig_ovs_Final_fuzzy_tune]);

y_min_fig = 0.5*y_min_fig;
y_max_fig = 1.1*y_max_fig;

figureIndex = auto_tune_param_all_plot_ovs_fluctStruct.figureIndex;

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
figure(figureIndex)
plot(ovs_Final_auto_tune,'b','LineWidth',1)
grid
title({'Automatic Tuning based on Overshoot:4.47%',...
       'Overshoot fluctuation while tuning - Controller with Real Zeros'})
axis([x_min_fig x_max_fig y_min_fig y_max_fig])

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
figureIndex = figureIndex + 1;
figure(figureIndex)
plot(ovs_Final_auto_tune_margaris,'b','LineWidth',1)
grid
title({'Automatic Tuning based on Overshoot: 4.47%',...
      'Overshoot fluctuation while tuning - Controller with Complex Zeros'})
axis([x_min_fig x_max_fig y_min_fig y_max_fig])

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
figureIndex = figureIndex + 1;
figure(figureIndex)
plot(ovs_Final_fuzzy_tune,'b','LineWidth',1)
grid on 
title({'Automatic Tuning based on Overshoot Estimation via FUZZY Model',...
      'Overshoot fluctuation while tuning - Controller with Complex Zeros'}) 
axis([x_min_fig x_max_fig y_min_fig y_max_fig])



