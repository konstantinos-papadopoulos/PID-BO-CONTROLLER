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

function [ovrst Fcl Gp] = auto_tune_2ndOrder_calculate_overshoot_2(plant_loc)

kh = plant_loc.kh            ;
Ti_2 = plant_loc.Ti_2        ;

T = plant_loc.T_arx          ;     
zeta = plant_loc.zeta_arx    ;       
  
kp_est = plant_loc.kp_est    ;
zeta_est = plant_loc.zeta_est;
T_est = plant_loc.T_est      ;
Ts_est = plant_loc.Ts_est    ;

num_Gp = kp_est                      ;
den_Gp = [T_est^2 2*zeta_est*T_est 1];
Gp = tf(num_Gp,den_Gp)               ;

% Fcl implementation 
% --------------------------------------------------------------------------------------------------
numFcl = [((T_est^2)*kp_est) (2*zeta_est*T_est*kp_est) kp_est];
denFcl = [(Ti_2*Ts_est*(T^2)) (Ti_2*(T^2) + 2*zeta*T*Ti_2*Ts_est) (2*zeta*T*Ti_2 + Ti_2*Ts_est+kh*kp_est*(T_est^2)) (Ti_2 + kh*kp_est*2*zeta_est*T_est) (kh*kp_est)];
Fcl = tf(numFcl,denFcl);

S = stepinfo(Fcl);
ovrst = S.Overshoot;
% --------------------------------------------------------------------------------------------------
% EOF:auto_tune_2ndOrder_calculate_overshoot_2
