% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_main_i.m																			   																		  	
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

function [Gc Gp stepinformation ti_MO x_MO] = auto_tune_param_main_pi_dig(plant_loc)

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
z0 = 1                                                         ;
z1 = tz1 + tz2 + tz3 + tz4                                     ;
z2 = tz1*tz2 + tz1*tz3 + tz1*tz4 + tz2*tz3 + tz2*tz4 + tz3*tz4 ;
z3 = tz1*tz2*tz3 + tz1*tz2*tz4 + tz2*tz3*tz4 + tz1*tz3*tz4     ;
z4 = tz1*tz2*tz3*tz4                                           ;
   
pGp0 = 1 ;
pGp1 = tp1 + tp2 + tp3 + tp4 + tp5                                               ;

pGp2 = tp1*tp2 +tp1*tp3 +tp1*tp4 +tp1*tp5 +tp2*tp3 + ...
     tp2*tp4 +tp2*tp5  +tp3*tp4 +tp3*tp5 +tp4*tp5                              ;

pGp3 = tp1*tp2*tp3 + tp1*tp2*tp4 + tp1*tp2*tp5 + tp1*tp3*tp4 + tp1*tp3*tp5 +...
     tp1*tp4*tp5 + tp2*tp3*tp4 + tp2*tp3*tp5 + tp2*tp4*tp5 + tp3*tp4*tp5       ;

pGp4 = tp1*tp2*tp3*tp4 + tp1*tp2*tp3*tp5 + tp1*tp2*tp4*tp5 + tp1*tp3*tp4*tp5 +...
     tp2*tp3*tp4*tp5                                                           ;
   
pGp5 = tp1*tp2*tp3*tp4*tp5                                                       ;
        
numGp = kp * [z4 z3 z2 z1 z0]               ;
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
    

% ------------------------------------------
n1 = factorial(1)  ;
n2 = 1/factorial(2);
n3 = 1/factorial(3);
n4 = 1/factorial(4);
% -----------------------------------------
d = td;
q0 = 1;
q1 = p1 + d;
q2 = p2 + p1*d + n2*(d^2);
q3 = p3 + p2*d + n2*(d^2)*p1 + n3*(d^3);
q4 = p4 + p3*d + n2*(d^2)*p2 + n3*(d^3)*p1 + n4*(d^4);
% -----------------------------------------
k1 = n1;
k2 = n2 + q1*n1;
k3 = n3 + q1*n2 + q2*n1;
k4 = q1*n3 + q2*n2 + q3*n1;
k5 = q2*n3 + q3*n2 + q4*n1;
% -----------------------------------------
A1 = k2/k1            ;
A2 = k2 - 2*(k1/k2)*k3;
A3 = k3/k2            ;
A4 = k4/k2            ;
% -----------------------------------------
x1 = n1                ;
x2 = z1*n1 + n2        ;
x3 = n3 + z1*n2 + z2*n1;
% -----------------------------------------

num_tk = A1*A2 + A4 - (A2 + A3)*z1 + z2 - A1^(-1)*z3;
den_tk = x1*(A2 + A3) - x2 + A1^(- 1)*x3            ;

tk = num_tk/den_tk;
tn = tk - 1       ;


x_MO  = tn        ;
y_MO  = 0         ;  
ti_MO = 2*kp*(p1 + td - z1 - x_MO - (1/2)) ;


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

% Closing the Loop
% ---------------------------------------------------------     
Ffp_MO = Gc_MO*Gp            ;
Fcl_MO = feedback(Ffp_MO,kh) ;

stepinformation = stepinfo(Fcl_MO);

% hold on
% grid
% -------------------------------------------------------------------------
% Complementary Sensitivities
% -------------------------------------------------------------------------  
So_MO =  1 - Fcl_MO             ; 
Si_MO =  series(So_MO,Gp)       ;
Su_MO = -kh*series(So_MO,Gc_MO) ;
Sr_MO =  series(So_MO,Gc_MO)    ;
