% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       autotune_param_typeII_calculateOvs.m																			   																		  	
%  Project:     Automatic tuning of the parameters for PI,PID controllers - TYPE II systems
%  
%  Purpose:     main script for the automatic tuning type II systems																		   																		
%  Author :     konstantina mermikli, kostas g. papadopoulos																	   																		
% 																										   																		
%  History:     Date: 25.06.2008  date last modified
% 																										  																		
%  Contact:     konstantina i. mermikli,    kostas g. papadopoulos,    
%               kmermikl@auth.gr       ,    kpapadop@eng.auth.gr
% 																										  																		
%  Place:	    Aristotle University of Thessaloniki, Thessaloniki, Greece							   																		
% 
% ************************************************************************************************************

function [Fcl_MO Fcl_optimal Gp stepinformation ti_MO x_MO y_MO] = autotune_param_typeII_calculateOvs(plant_loc)

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
tn = plant_loc.tn  ;

% -------------------------------------------------------------------------
% Plant Implementation
% -------------------------------------------------------------------------
q0 = 1                                                         ;
q1 = tz1 + tz2 + tz3 + tz4                                     ;
q2 = tz1*tz2 + tz1*tz3 + tz1*tz4 + tz2*tz3 + tz2*tz4 + tz3*tz4 ;
q3 = tz1*tz2*tz3 + tz1*tz2*tz4 + tz2*tz3*tz4 + tz1*tz3*tz4     ;
q4 = tz1*tz2*tz3*tz4                                           ;
   
pGp0 = 1                                                       ;
pGp1 = tp1 + tp2 + tp3 + tp4 + tp5                             ;

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
    Gd = tf(num_d,den_d)       ;      
    Gp = Gp1 * Gd              ;      % Transfer Function of the Process with Time Delay
    n1 = td                    ;
    n2 = td^2 / 2              ;
else
    Gd = 1   ;
    Gp = Gp1 ;        % Transfer Function of the Process without Time Delay
    n1 = 0   ;
    n2 = 0   ;
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

% 1st Method (...leading to 4 roots): eq.47
% -------------------------------------------------------------------------
% tn^4 : 1
% tn^3 : 4*(z1 - q1)
% tn^2 : 6*z1^2 - 8*z2 + 4*q2 - 4*q1*z1
% tn^1 : 4*(z1 - q1)*(z1^2 - 2*z2)
% tn^0 : (z1^2 + 2*z2 + 4*q2 - 4*q1*z1)*(z1^2 - 2*z2) + 4*(z2^2 - 2*z1*z3 + 2z4)

% z0 = q0;
% z1 = q1;
% z2 = q2;
% z3 = q3;
% z4 = q4;


c_meth_1(1) = 1                                                                       ;   
c_meth_1(2) = 4*(q1 - w1)                                                             ;
c_meth_1(3) = 6*q1^2 - 8*q2 + 4*w2 - 4*w1*q1                                          ;
c_meth_1(4) = 4*(q1 - w1)*(q1^2 - 2*q2)                                               ;
c_meth_1(5) = (q1^2 + 2*q2 + 4*w2 - 4*w1*q1)*(q1^2 - 2*q2) + 4*(q2^2 - 2*q1*q3 + 2*q4);

tn_meth1 = roots(c_meth_1);

% collecting candidate solutions
% -------------------------------------------------------------------------
tn_meth1_candidate = -15*ones(4,1);
m = 1;
for i = 1:4
    if isreal(tn_meth1(i)) == 1
        tn_meth1_candidate(m) = tn_meth1(i);
        m = m + 1                          ;
    end
end
tn_meth1_sol = max(tn_meth1_candidate);

if tn_meth1_sol == 0
   tn_meth1_sol = -100;
end

% 2nd Method (...leading to 2 roots): eq.54
% -------------------------------------------------------------------------
% tn^2 : 1
% tn^1 : 4*(z1 - q1)
% tn^0 : 4*q2 + 2*z2 - 4*q1*z1 + z1^2

c_meth_2(1) = 1                               ;   
c_meth_2(2) = 4*(q1 - w1)                     ;
c_meth_2(3) = 4*w2 + 2*q2 - 4*w1*q1 + q1^2    ;
tn_meth2 = roots(c_meth_2)                    ;                    

if isreal(tn_meth2) == 1
    tn_meth2_1 = tn_meth2(1)                  ;
    tn_meth2_2 = tn_meth2(2)                  ;
    tn_meth2_sol = max(tn_meth2_1,tn_meth2_2) ;
else
    disp('...complex solution for tn: METHOD 2 - 2:roots')
end

% Choose the optimal controller
% -------------------------------------------------------------------------
tn_optimal = tn_meth1_sol  ;
% -------------------------------------------------------------------------

x_optimal  = tn_optimal            ;
y_optimal  = 0                     ;  
ti_optimal = sqrt((1/2)*kp*kh*(tn_optimal^2 + q1^2 - 2*q2));
ti_optimal_square = ti_optimal^2   ;


num0Gc_opt = 1                            ;
num1Gc_opt = x_optimal                    ;
num2Gc_opt = y_optimal                    ;
numGc_optimal = [num2Gc_opt num1Gc_opt num0Gc_opt];
    
den0Gc_opt = 0                       ;
den1Gc_opt = 0                       ;
den2Gc_opt = ti_optimal_square            ;
den3Gc_opt = ti_optimal_square*tp6        ;

denGc_optimal = [den3Gc_opt den2Gc_opt den1Gc_opt den0Gc_opt];
Gc_optimal = tf(numGc_optimal,denGc_optimal)    ; 
Gc = Gc_optimal                       ;
% -------------------------------------------------------------------------
% Forward Path and Closed Loop Implementation
% -------------------------------------------------------------------------  
Ffp_optimal = Gc_optimal*Gp            ;
Fcl_optimal = feedback(Ffp_optimal,kh) ;

% Automatic tuning
% -------------------------------------------------------------------------  
tn_MO = tn;

x_MO  = tn_MO                      ;
y_MO  = 0                          ;  
% ti_MO = sqrt((1/2)*kp*kh*(tn_MO^2 + q1^2 - 2*q2));
ti_MO = sqrt((1/2)*kp*kh*(tn_MO^2));
ti_MO_square = ti_MO^2             ;


num0Gc = 1                       ;
num1Gc = x_MO                    ;
num2Gc = y_MO                    ;
numGc_MO = [num2Gc num1Gc num0Gc];
    
den0Gc = 0                       ;
den1Gc = 0                       ;
den2Gc = ti_MO_square            ;
den3Gc = ti_MO_square*tp6        ;

denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
Gc_MO = tf(numGc_MO,denGc_MO)    ; 
Gc = Gc_MO                       ;
% -------------------------------------------------------------------------
% Forward Path and Closed Loop Implementation
% -------------------------------------------------------------------------  
Ffp_MO = Gc_MO*Gp            ;
Fcl_MO = feedback(Ffp_MO,kh) ;
stepinformation = stepinfo(Fcl_MO);
% hold on
% grid        
% -------------------------------------------------------------------------
% Complementary Sensitivities
% -------------------------------------------------------------------------  
% So_MO =  1 - Fcl_MO             ; 
% Si_MO =  series(So_MO,Gp)       ;
% Su_MO = -kh*series(So_MO,Gc_MO) ;
% Sr_MO =  series(So_MO,Gc_MO)    ;
