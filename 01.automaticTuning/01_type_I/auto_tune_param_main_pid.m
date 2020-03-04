% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_main_pid.m																			   																		  	
%  Project:     Automatic tuning of the paramters for PI,PID controllers
%  
%  Purpose:     calculate "pid" controller parameters																		   																		
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

function [Gc Gp stepinformation ti_MO x_MO y_MO stepinformation_Gp] = auto_tune_param_main_pid(plant_loc)
% Coefficients for Pade Approximation
% -------------------------------------------------------------------------
 N = 15                     ;
% -------------------------------------------------------------------------

tp1 = plant_loc.tp1;    tz1 = plant_loc.tz1; 
tp2 = plant_loc.tp2;    tz2 = plant_loc.tz2;
tp3 = plant_loc.tp3;    tz3 = plant_loc.tz3;
tp4 = plant_loc.tp4;    tz4 = plant_loc.tz4;
tp5 = plant_loc.tp5;   
tp6 = plant_loc.tp6;   

td = plant_loc.td  ;    kp = plant_loc.kp  ;    kh = plant_loc.kh  ;


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


% -------------------------------------------------------------------------
% Controller Parameter Calculation
% -------------------------------------------------------------------------
p1 = tp1 + tp2 + tp3 + tp4 + tp5 + tp6                                    ;

p2 = tp1*tp2 + tp1*tp3 + tp1*tp4 + tp1*tp5 + tp1*tp6  + ...
     tp2*tp3 + tp2*tp4 + tp2*tp5 + tp2*tp6 + tp3*tp4  + ...
     tp3*tp5 + tp3*tp6 + tp4*tp5 + tp4*tp6 + tp5*tp6                      ;
 
p3 = tp1*tp2*tp3 + tp1*tp2*tp4 + tp1*tp2*tp5 + tp1*tp2*tp6 + ...
     tp1*tp3*tp4 + tp1*tp3*tp5 + tp1*tp3*tp6 + tp1*tp4*tp5 + ...
     tp1*tp4*tp6 + tp1*tp5*tp6 + tp2*tp3*tp4 + tp2*tp3*tp5 + ...
     tp2*tp3*tp6 + tp2*tp4*tp5 + tp2*tp4*tp6 + tp2*tp5*tp6 + ...
     tp3*tp4*tp5 + tp3*tp4*tp6 + tp3*tp5*tp6 + tp4*tp5*tp6                ;
 
p4 = tp1*tp2*tp3*tp4 + tp1*tp2*tp3*tp5 + tp1*tp2*tp3*tp6 + ...
     tp1*tp2*tp4*tp5 + tp1*tp2*tp4*tp6 + tp1*tp2*tp5*tp6 + ... 
     tp1*tp3*tp4*tp5 + tp1*tp3*tp4*tp6 + tp1*tp3*tp5*tp6 + ...
     tp1*tp4*tp5*tp6 + tp2*tp3*tp4*tp5 + tp2*tp3*tp4*tp6 + ...
     tp2*tp3*tp5*tp6 + tp2*tp4*tp5*tp6 + tp3*tp4*tp5*tp6                  ;
 
p5 = tp1*tp2*tp3*tp4*tp5 + tp1*tp2*tp3*tp4*tp6 + tp1*tp2*tp3*tp5*tp6 +...
     tp1*tp2*tp4*tp5*tp6 + tp1*tp3*tp4*tp5*tp6 + tp2*tp3*tp4*tp5*tp6      ;
 
p6 = tp1*tp2*tp3*tp4*tp5*tp6                                              ;
    
w1 = n1 + p1            ;
w2 = n2 + n1*p1 +p2     ;
w3 = n2*p1 + n1*p2 + p3 ;
w4 = n2*p2 + n1*p3 + p4 ;
w5 = n2*p3 + n1*p4 + p5 ;
w6 = n2*p4 + n1*p5 + p6 ;
w7 = n2*p5 + n1*p6      ;
w8 = n2*p6              ;



A1 = w1^2- w2 - w1*q1 + q2                         ;
B1 = q1 -  w1                                      ;
C1 = w1^3 - w1^2*q1 - 2*w1*w2 + w2*q1 + w3 + w1*q2 - q3 ;

A2 = w2^2 + w1*q3 + w4 + w3*q1 - w2*q2 - 2*w1*w3 - q4 ;
B2 = w3 + w1*q2 - w2*q1  - q3 ;
C2 = (w1 - q1)*(w2^2 - 2*w1*w3 + 2*w4) - (w1*q4 + w3*q2 + w5 - w4*q1 - w2*q3);


C = [C1;C2]       ;
D = [A1 B1;A2 B2] ;
E = inv(D)        ;
F = E*C           ; 

x_MO  = F(1,1)               ;
y_MO  = F(2,1)               ;  
ti_MO = 2*kp*(w1 - q1 - x_MO);
% alternate calculation of x_MO

num0Gc = 1                       ;
num1Gc = x_MO                    ;
num2Gc = y_MO                    ;
numGc_MO = [num2Gc num1Gc num0Gc];
    
den0Gc = 0                       ;
den1Gc = ti_MO                   ;
den2Gc = ti_MO*tp6               ;
denGc_MO = [den2Gc den1Gc den0Gc];
Gc_MO = tf(numGc_MO,denGc_MO)    ; 
Gc = Gc_MO                       ;
% -------------------------------------------------------------------------
% Forward Path and Closed Loop Implementation
% -------------------------------------------------------------------------  
Ffp_MO = Gc_MO*Gp            ;
Fcl_MO = feedback(Ffp_MO,kh) ;
stepinformation = stepinfo(Fcl_MO);
stepinformation_Gp = stepinfo(Gp);
% figure(2)
% step(Fcl_MO);
% hold on
% grid        

% -------------------------------------------------------------------------
% Complementary Sensitivities
% -------------------------------------------------------------------------  
So_MO =  1 - Fcl_MO             ; 
Si_MO =  series(So_MO,Gp)       ;
Su_MO = -kh*series(So_MO,Gc_MO) ;
Sr_MO =  series(So_MO,Gc_MO)    ;
% step(Fcl_MO,'r')
% hold on
