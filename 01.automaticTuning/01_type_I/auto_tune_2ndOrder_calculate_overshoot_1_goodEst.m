% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_2ndOrder_calculate_overshoot_1_goodEst.m																			   																		  	
%  Project:     Automatic tuning of the parameters for PI,PID controllers
%  
%  Purpose:     Ti1_tuning for the second order system																	   																		
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

function [ovrst_good Fcl_good Ti_good] = auto_tune_2ndOrder_calculate_overshoot_1_goodEst(plant_loc)

kh = plant_loc.kh    ;
Tsx = plant_loc.Tsx  ; 

% restore the estimation without error
% ---------------------------------------
T_good = plant_loc.T_est_good      ;    
zeta_good = plant_loc.zeta_est_good;       
kp_good = plant_loc.kp_est_good    ;  

Ti_good = 2*kp_good*kh*(2*zeta_good*T_good + Tsx);

% Closed Loop Transfer Function after the estimation and the open loop experiment
% --------------------------------------------------------------------------------------------------
% Estimation without error
% ---------------------------------
numFcl_good = kp_good;
denFcl_good = [((T_good^2)*Ti_good*Tsx) (T_good*Ti_good*(T_good + 2*zeta_good*Tsx)) ((Tsx + 2*zeta_good*T_good)*Ti_good) Ti_good (kh*kp_good)];
Fcl_good = tf(numFcl_good,denFcl_good);

S = stepinfo(Fcl_good)       ;
S_good = stepinfo(Fcl_good)  ;
ovrst_good = S_good.Overshoot;

% --------------------------------------------------------------------------------------------------
% EOF:auto_tune_2ndOrder_calculate_overshoot_1_goodEst
