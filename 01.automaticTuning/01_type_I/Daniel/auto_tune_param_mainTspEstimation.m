% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_mainTspEstimation.m																			   																		  	
%  Project:     Automatic tuning of the paramters for PI,PID controllers
%  
%  Purpose:     main script for the automatic tuning based on the three methods																		   																		
%  Author :     kostas g. papadopoulos																	   																		
% 																										   																		
%  History:     Date: 07.07.2008  date last modified
% 																										  																		
%  Contact:     kostas g. papadopoulos,    nikos mitrakis,       leonidas droukas
%               kpapadop@eng.auth.gr  ,    nmitr@auth.gr ,       leon_drouk@yahoo.gr
% 																										  																		
%  Place:	    Aristotle University of Thessaloniki, Thessaloniki, Greece							   																		
% 
% ************************************************************************************************************

function [tsp kp_e tss Gp] = auto_tune_param_mainTspEstimation(plant_loc)

% Assuming that the plant is a First Order Plant
% *************************************************************************
tp1 = plant_loc.tp1;     
tp2 = plant_loc.tp2;    
tp3 = plant_loc.tp3;    
tp4 = plant_loc.tp4;    
tp5 = plant_loc.tp5;   
td = plant_loc.td  ;      kp_r = plant_loc.kp;      

% Summing Up the Dominant Time constant
% *************************************************************************
pGp1 = tp1 + tp2 + tp3 + tp4 + tp5;
pGp2 = tp1*tp2 + tp1*tp3 + tp2*tp3;
pGp3 = tp1*tp2*tp3;
pGp4 = 0;
pGp5 = 0;

q4 = 0;      
q3 = 0;    
q2 = 0;     
q1 = 0;

a = 0 ;
kp_r = (1 + a)*kp_r   ;
num_Gp_s = kp_r*[q4 q3 q2 q1 1]            ;
den_Gp_s = [pGp5 pGp4 pGp3 pGp2 pGp1 pGp1] ;
Gp_s     = tf(num_Gp_s,den_Gp_s)           ;
Gp1 = Gp_s;
% define time delay constant at the output of the process
% *************************************************************************
if (td ~= 0)                          
    N = 20 ;
    [num_d,den_d] = pade(td,N) ;
    Gd = tf(num_d,den_d);
else
    Gd = 1 ;
end
Gp = Gp1 * Gd;
kp_e = kp_r;
S = stepinfo(Gp)    ;
tss = S.SettlingTime;
tsp = tss / 4       ;           % approximate calculation of tsp
                                % based on the settling time of the open
                                % loop experiment
% *************************************************************************
