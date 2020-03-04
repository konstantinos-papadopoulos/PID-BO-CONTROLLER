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
pGp2 = 0;
pGp3 = 0;
pGp4 = 0;
pGp5 = 0;

q4 = 0;      
q3 = 0;    
q2 = 0;     
q1 = 0;

numGp = kp_r;
denGp = [pGp1 1];
Gp1 = tf(numGp,denGp);

% orismos tis vathmidas xronikis kathisterisis
% *************************************************************************
if (td ~= 0)                          
    N = 20 ;
    [num_d,den_d] = pade(td,N) ;
    Gd = tf(num_d,den_d);
else
    Gd = 1 ;
end
Gp = Gp1 * Gd;

S = stepinfo(Gp)    ;
tss = S.SettlingTime;
tsp = tss / 4       ;           %ipologizetai proseggistika to tsp, me vasi to settling time tis apokrisis

% *************************************************************************
syms s;
Gp_s = (kp_r * ((q4*(s^4)) + (q3*(s^3)) + (q2*(s^2)) + (q1*s) + 1)) /...
       ((pGp5*(s^5)) + (pGp4*(s^4)) + (pGp3*(s^3)) + (pGp2*(s^2)) + (pGp1*s) + 1);
        
   r = limit(Gp_s,s,0);
   kp_e = double(r);
