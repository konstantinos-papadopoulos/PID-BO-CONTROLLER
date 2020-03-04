% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_main.m																			   																		  	
%  Project:     Automatic tuning of the parameters for PI,PID controllers
%  
%  Purpose:     main script for the automatic tuning based on the three methods																		   																		
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
% H sinartisi dexetai os eisodo to sistima Gp kai ipologizei mia
% proseggistiki timi gia ta kp , J kai T
% -------------------------------------------------------------------------------------------------
function [kp_estim zeta_estim T_estim kp_estim_good zeta_estim_good T_estim_good] = auto_tune_2ndOrder_kp_zeta_T_estimation(plant_loc)

T = plant_loc.T_arx      ;     
Ts = plant_loc.Ts_arx    ; 
zeta = plant_loc.zeta_arx;    kp = plant_loc.kp_arx;

% assume that the second order term is conceived in the open loop experiment
% --------------------------------------------------------------------------------------------------
numGp = kp                  ;
denGp = [(T^2) (2*zeta*T) 1];
Gp = tf(numGp,denGp)        ;


% Estimating kp from the open loop step response
% --------------------------------------------------------------------------------------------------
syms s;
Gp_s = kp / (1 + 2*zeta*T*s + (T^2)*(s^2));
r = limit(Gp_s,s,0)                       ;
kp_estim = double(r)                      ;
kp_estim_good = kp_estim                  ;

% zeta estimation
% --------------------------------------------------------------------------------------------------
[y_r,t] = step(Gp)                                       ;
M = ( (max(y_r)) / kp_estim ) - 1                        ;
zeta_estim = sqrt( ((log(M))^2) / ((pi^2)+((log(M))^2)) );
zeta_estim_good = zeta_estim                             ;

% T estimation
% --------------------------------------------------------------------------------------------------
F = stepinfo(Gp)                            ;
tss_e = F.SettlingTime                      ;
T_estim_good = (tss_e / 4) * zeta_estim_good;


flagError_kp = 0  ;
flagError_zeta = -0.15;
flagError_T = 0   ;

kp_estim   = (1 + flagError_kp)*plant_loc.kp_arx    ;
zeta_estim = (1 + flagError_zeta)*plant_loc.zeta_arx;
T_estim    = (tss_e / 4) * zeta_estim               ;

