% ************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% **************************************************************************************************
%  Title:       plotResults.m																			   																		  	
%  Project:     Automatic tuning of the parameters for PI,PID controllers
%  
%  Purpose:     main script for the automatic tuning based on the three methods																		   																		
%  Author :     kostas g. papadopoulos																	   																		
% 																										   																		
%  History:     Date: 25.06.2008  date last modified
% 																										  																		
%  Contact:     kostas g. papadopoulos,    nikos mitrakis       leonidas droukas
%               kpapadop@eng.auth.gr  ,    nmitr@auth.gr        leon_drouk@yahoo.gr
% 																										  																		
%  Place:	    Aristotle University of Thessaloniki, Thessaloniki, Greece							   																		
% 
% **************************************************************************************************
ltiview(Fcl_MO_PID_optimal,'k',Fcl_auto_PID,'--r',Fcl_PID_auto_marg,'-.b')

% ++++++++++++++++++++++++++++++ Figure 7 +++++++++++++++++++++++++++++++++
x_min_fig7 = 0; 
x_max_fig7 = 3*FinalTime;
y_min_fig7 = 0;
sortMAX_I = [y_optimal_Si_I y_autotuned_Si_I y_fuzzy_autotuned_Si_I];
sortMAX_I = [y_optimal_Si_I y_autotuned_Si_I];
chooseMaxItemp = max(max(sortMAX_I));
chooseMaxI = max(2,1 + chooseMaxItemp);
if chooseMaxI == 2
    y_max_fig7 = 2.1;
else
    y_max_fig7 = 1.05*chooseMaxI;
end


figure(figureIndex)
plot(t1,y_optimal_I,'k','LineWidth',2.5)
hold on
plot(t1,y_autotuned_I,'--r','LineWidth',2.5)
hold on
plot(t1,y_fuzzy_autotuned_I,'g')
axis([x_min_fig7 x_max_fig7 y_min_fig7 y_max_fig7])
grid
legend('optimal','autotuned [based on ovs]','fuzzy autotuned [based on ovs]',0)
title('I Controller')

% ++++++++++++++++++++++++++++++ Figure 8 +++++++++++++++++++++++++++++++++
x_min_fig8 = 0; 
x_max_fig8 = 3*FinalTime;
y_min_fig8 = 0;
sortMAX_PI = [y_optimal_Si_PI y_autotuned_Si_PI y_fuzzy_autotuned_Si_PI];
chooseMaxPItemp = max(max(sortMAX_PI));
chooseMaxPI = max(2,1 + chooseMaxPItemp);
if chooseMaxPI == 2
    y_max_fig8 = 2.1;
else
    y_max_fig8 = 1.05*chooseMaxPI;
end
    
figureIndex = figureIndex + 1;
figure(figureIndex)
plot(t1,y_optimal_PI,'k','LineWidth',2.5)
hold on
plot(t1,y_autotuned_PI,'--r','LineWidth',2.5)
hold on
plot(t1,y_fuzzy_autotuned_PI,'-.g','LineWidth',2.5)


axis([x_min_fig8 x_max_fig8 y_min_fig8 y_max_fig8])
grid
legend('optimal','autotuned [based on ovs]','fuzzy autotuned',0)
title('PI Controller')


% ++++++++++++++++++++++++++++++ Figure 9 +++++++++++++++++++++++++++++++++
x_min_fig9 = 0; 
x_max_fig9 = 3*FinalTime;
y_min_fig9 = 0;
% sortMAX_PID = [y_optimal_Si_PID y_autotuned_Si_PID y_fuzzy_autotuned_Si_PID y_auto_autotuned_marg_Si_PID];
sortMAX_PID = [y_optimal_Si_PID y_autotuned_Si_PID y_auto_autotuned_marg_Si_PID y_fuzzy_autotuned_Si_PID];
chooseMaxPIDtemp = max(max(sortMAX_PID));
chooseMaxPID = max(2,1 + chooseMaxPIDtemp);


if chooseMaxPID == 2
   y_max_fig9 = 2.1;
else
   y_max_fig9 = 1.05*chooseMaxPID;
end

figureIndex = figureIndex + 1;
figure(figureIndex)

plot(t1,y_optimal_PID,'k','LineWidth',2.5)
hold on

plot(t1,y_autotuned_PID,'--r','LineWidth',2.5)
hold on

plot(t1,y_fuzzy_autotuned_PID,'-.g','LineWidth',2.5)
hold on

plot(t1,y_auto_marg_PID,'-.b','LineWidth',2.5)

axis([x_min_fig9 x_max_fig9 y_min_fig9 y_max_fig9])
grid
legend('optimal','autotuned [based on ovs]','fuzzy autotuned','autotuned [based on ovs]: "y" estimation',0)
title('PID Controller')

% PhD Figures
% ------------------------------------------------------------------
% Closed Loop Control system Output
% =================================
figureIndex = figureIndex + 1;
figure(figureIndex)
plot(t1,y_optimal_Fcl_PID,'k',t1,y_autotuned_Fcl_PID,'r',t1,y_fuzzy_autotuned_Fcl_PID,'b','LineWidth',2.5)
grid on
title('PID controller')
legend('PID optimal','PID autotuned',0)


figureIndex = figureIndex + 1;
figure(figureIndex)
plot(t2,y_optimal_So_PID,'k',t2,y_autotuned_So_PID,'r',t2,y_fuzzy_autotuned_So_PID,'LineWidth',2.5)
grid on
title('PID controller - Output Disturbances')
legend('PID optimal','PID autotuned',0)

figureIndex = figureIndex + 1;
figure(figureIndex)
plot(t3,y_optimal_Si_PID,'k',t3,y_autotuned_Si_PID,'r',t3,y_fuzzy_autotuned_Si_PID,'b','LineWidth',2.5)
grid on
title('PID controller - Input Disturbances')
legend('PID optimal','PID autotuned',0)


% Closed Loop Control system Output
% PI Controller
% =================================
figureIndex = figureIndex + 1;
figure(figureIndex)
plot(t1,y_optimal_Fcl_PI,'k',t1,y_autotuned_Fcl_PI,'r',t1,y_fuzzy_autotuned_Fcl_PI,'b','LineWidth',2.5)
grid on
title('PI controller')
legend('PI optimal','PI autotuned',0)

figureIndex = figureIndex + 1;
figure(figureIndex)
plot(t2,y_optimal_So_PI,'k',t2,y_autotuned_So_PI,'r',t2,y_fuzzy_autotuned_So_PI,'b','LineWidth',2.5)
grid on
title('PI controller - Output Disturbances')
legend('PI optimal','PI autotuned',0)

figureIndex = figureIndex + 1;
figure(figureIndex)
plot(t3,y_optimal_Si_PI,'k',t3,y_autotuned_Si_PI,'r',t3,y_fuzzy_autotuned_Si_PI,'b','LineWidth',2.5)
grid on
title('PI controller - Input Disturbances')
legend('PI optimal','PI autotuned',0)


ovs_PI_optimal = optimal_PI_struct.Overshoot  ;
ovs_PID_optimal = optimal_PID_struct.Overshoot;
ovs_PI_auto = autotune_PI_struct.Overshoot    ;
ovs_PID_auto = autotune_PID_struct.Overshoot  ;

tss_PI_optimal = optimal_PI_struct.SettlingTime  ;
tss_PID_optimal = optimal_PID_struct.SettlingTime;
tss_PI_auto = autotune_PI_struct.SettlingTime    ;
tss_PID_auto = autotune_PID_struct.SettlingTime  ;


fprintf('ovsPI_optimal: %1.5f,  ovsPI_auto: %1.5f\n',ovs_PI_optimal,ovs_PI_auto)
fprintf('ovsPID_optimal: %1.5f, ovsPID_auto: %1.5f\n',ovs_PID_optimal,ovs_PID_auto)

fprintf('tssPI_optimal : %1.5f,  tssPI_auto: %1.5f\n',tss_PI_optimal,tss_PI_auto)
fprintf('tssPID_optimal: %1.5f, tssPID_auto: %1.5f\n',tss_PID_optimal,tss_PID_auto)
 
ltiview(Si_MO_PI_optimal,'k',Si_MO_PI_autotuned,'--r',Si_MO_PI_fuzzy_autotuned,'-.b')
ltiview(Si_MO_PID_optimal,'k',Si_MO_PID_autotuned,'--r',Si_MO_PID_fuzzy_autotuned,'-.b')
ltiview(Fcl_MO_PI_optimal,'k',Fcl_auto_PI,'r',Fcl_PI_fuzzy,'b')
ltiview(Fcl_MO_PID_optimal,'k',Fcl_auto_PID,'r',Fcl_PID_fuzzy,'b')







