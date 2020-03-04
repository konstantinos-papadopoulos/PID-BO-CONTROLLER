% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_2ndOrder_calculate_overshoot_1.m																			   																		  	
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

function [ovrst Fcl Ti] = auto_tune_2ndOrder_calculate_overshoot_1(plant_loc)

kh = plant_loc.kh    ;
Tsx = plant_loc.Tsx  ; 

% restore the estimation with error
% ---------------------------------------
T = plant_loc.T_est      ;    
zeta = plant_loc.zeta_est;       
kp = plant_loc.kp_est    ;  

Ti = 2*kp*kh*(2*zeta*T + Tsx);

% Closed Loop Transfer Function after the estimation and the open loop experiment
% --------------------------------------------------------------------------------------------------
% Estimation with error
% ---------------------------------
numFcl = kp;
denFcl = [((T^2)*Ti*Tsx) (T*Ti*(T+2*zeta*Tsx)) ((Tsx+2*zeta*T)*Ti) Ti (kh*kp)];
Fcl = tf(numFcl,denFcl);

S = stepinfo(Fcl)            ;
ovrst = S.Overshoot          ;

% --------------------------------------------------------------------------------------------------
% EOF:auto_tune_2ndOrder_calculate_overshoot_1
