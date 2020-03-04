% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_main_SettlingTimePlant.m																			   																		  	
%  Project:     Automatic tuning of the paramters for PI,PID controllers
%  
%  Purpose:     calculate "i" controller parameters																		   																		
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

function stepinformation = auto_tune_param_main_SettlingTimePlant(plant_loc)

tp1 = plant_loc.tp1;    tz1 = plant_loc.tz1; 
tp2 = plant_loc.tp2;    tz2 = plant_loc.tz2;
tp3 = plant_loc.tp3;    tz3 = plant_loc.tz3;
tp4 = plant_loc.tp4;    tz4 = plant_loc.tz4;
tp5 = plant_loc.tp5;   
td = plant_loc.td  ;    kp = plant_loc.kp  ; 
% -------------------------------------------------------------------------
% Plant Implementation
% -------------------------------------------------------------------------
q0 = 1                                                         ;
q1 = tz1 + tz2 + tz3 + tz4                                     ;
q2 = tz1*tz2 + tz1*tz3 + tz1*tz4 + tz2*tz3 + tz2*tz4 + tz3*tz4 ;
q3 = tz1*tz2*tz3 + tz1*tz2*tz4 + tz2*tz3*tz4 + tz1*tz3*tz4     ;
q4 = tz1*tz2*tz3*tz4                                           ;
   
pGp0 = 1 ;
pGp1 = tp1 + tp2 + tp3 + tp4 + tp5                                               ;

pGp2 = tp1*tp2 +tp1*tp3 +tp1*tp4 +tp1*tp5 +tp2*tp3 + ...
       tp2*tp4 +tp2*tp5  +tp3*tp4 +tp3*tp5 +tp4*tp5                              ;

pGp3 = tp1*tp2*tp3 + tp1*tp2*tp4 + tp1*tp2*tp5 + tp1*tp3*tp4 + tp1*tp3*tp5 +...
       tp1*tp4*tp5 + tp2*tp3*tp4 + tp2*tp3*tp5 + tp2*tp4*tp5 + tp3*tp4*tp5       ;

pGp4 = tp1*tp2*tp3*tp4 + tp1*tp2*tp3*tp5 + tp1*tp2*tp4*tp5 + tp1*tp3*tp4*tp5 +...
       tp2*tp3*tp4*tp5                                                           ;
   
pGp5 = tp1*tp2*tp3*tp4*tp5                                                       ;
        
numGp = kp * [q4 q3 q2 q1 q0]               ;
denGp  = [pGp5 pGp4 pGp3 pGp2 pGp1 pGp0]    ;
Gp1 = tf(numGp,denGp)                       ;

% -------------------------------------------------------------------------
% Pade Approximation 
% -------------------------------------------------------------------------  
if td ~= 0
    N = 30 ;
    [num_d,den_d] = pade(td,N) ;
    Gd = tf(num_d,den_d);      
    Gp = Gp1 * Gd ;      % Transfer Function of the Process with Time Delay
    n1 = td ;
    n2 = td^2 / 2 ;
else
    Gd = 1 ;
    Gp = Gp1 ;        % Transfer Function of the Process without Time Delay
    n1 = 0 ;
    n2 = 0 ;
end
  
stepinformation = stepinfo(Gp);

