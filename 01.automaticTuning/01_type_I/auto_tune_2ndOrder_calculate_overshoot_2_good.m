% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_2ndOrder_calculate_overshoot_2.m																			   																		  	
%  Project:     Automatic tuning of the parameters for PI,PID controllers
%  
%  Purpose:     Ti2_tuning for the second order system																	   																		
%  Author :     kostas g. papadopoulos																	   																		
% 																										   																		
%  History:     Date: 07.07.2008  date last modified
% 																										  																		
%  Contact:     leonidas droukas   ,       kostas g. papadopoulos    
%               leon_drouk@yahoo.gr,       kpapadop@eng.auth.gr      
% 																										  																		
%  Place:	    Aristotle University of Thessaloniki, Thessaloniki, Greece							   																		
% 
% ************************************************************************************************************

function [ovrst Fcl Gp] = auto_tune_2ndOrder_calculate_overshoot_2_good(plant_loc)

kh = plant_loc.kh            ;
Ti_2_good = plant_loc.Ti_2_good        ;

T = plant_loc.T_arx          ;     
zeta = plant_loc.zeta_arx    ;       
  
kp_est_good = plant_loc.kp_est_good    ;
zeta_est_good = plant_loc.zeta_est_good;
T_est_good = plant_loc.T_est_good      ;
Ts_est_good = plant_loc.Ts_est_good    ;

num_Gp = kp_est_good                      ;
den_Gp = [T_est_good^2 2*zeta_est_good*T_est_good 1];
Gp = tf(num_Gp,den_Gp)               ;

% Fcl implementation 
% --------------------------------------------------------------------------------------------------
numFcl = [((T_est_good^2)*kp_est_good) (2*zeta_est_good*T_est_good*kp_est_good) kp_est_good];
denFcl = [(Ti_2_good*Ts_est_good*(T^2)) (Ti_2_good*(T^2) + 2*zeta*T*Ti_2_good*Ts_est_good) (2*zeta*T*Ti_2_good +...
           Ti_2_good*Ts_est_good + kh*kp_est_good*(T_est_good^2)) (Ti_2_good + kh*kp_est_good*2*zeta_est_good*T_est_good) (kh*kp_est_good)];
Fcl = tf(numFcl,denFcl);

S = stepinfo(Fcl);
ovrst = S.Overshoot;
% --------------------------------------------------------------------------------------------------
% EOF:auto_tune_2ndOrder_calculate_overshoot_2
