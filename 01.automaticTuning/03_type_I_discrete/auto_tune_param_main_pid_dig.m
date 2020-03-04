% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_main_pid_dig.m																			   																		  	
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

function [Gc Gp stepinformation ti_MO x_MO y_MO] = auto_tune_param_main_pid_dig(plant_loc)

% Coefficients for Pade Approximation
% -------------------------------------------------------------------------
 N = 15            ;
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
z5 = 0 ;

pGp0 = 1 ;
pGp1 = tp1 + tp2 + tp3 + tp4 + tp5                                             ;

pGp2 = tp1*tp2 +tp1*tp3 +tp1*tp4 +tp1*tp5 +tp2*tp3 + ...
     tp2*tp4 +tp2*tp5  +tp3*tp4 +tp3*tp5 +tp4*tp5                              ;

pGp3 = tp1*tp2*tp3 + tp1*tp2*tp4 + tp1*tp2*tp5 + tp1*tp3*tp4 + tp1*tp3*tp5 +...
     tp1*tp4*tp5 + tp2*tp3*tp4 + tp2*tp3*tp5 + tp2*tp4*tp5 + tp3*tp4*tp5       ;

pGp4 = tp1*tp2*tp3*tp4 + tp1*tp2*tp3*tp5 + tp1*tp2*tp4*tp5 + tp1*tp3*tp4*tp5 +...
     tp2*tp3*tp4*tp5                                                           ;
   
pGp5 = tp1*tp2*tp3*tp4*tp5                                                     ;
        
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
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
r1 = 2*n1            ;
r2 = (n1^2) + 2*n2   ;   
r3 = 2*(n3 + n1*n2)  ;
r4 = (n2^2) + 2*n1*n3;
r5 = 2*n2*n3         ;
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
a11 = r1           ;    a12 = -(n1)           ;
a21 = z0*r2 + z1*r1;    a22 = -(z1*n1 + z0*n2);
a31 = z1*r2 + z2*r1;    a32 = -(z2*n1 + z1*n2);
a41 = z2*r2 + z3*r1;    a42 = -(z3*n1 + z2*n2);
a51 = z3*r2 + z4*r1;    a52 = -(z4*n1 + z3*n2);

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
u1 = r1 - n1;
u2 = r2 - n2;
u3 = r3 - n3;
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
w1 = u1                   ;
w2 = u2    + u1*q1        ;
w3 = u3    + u2*q1 + u1*q2;
w4 = u3*q1 + u2*q2 + u1*q3;
w5 = u3*q2 + u2*q3 + u1*q4;
w6 = u3*q3 + u2*q4        ;
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
A = (w2/w1^2) - (z1/w1)         ;
B = w2^2 - 2*w1*w3              ;
C = (w3^2) - 2*( w2*w4 - w1*w5 );

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
h1 = - a12/w1;
h2 = - a11/w1;

A11 = w1*a32 + w3*a12 - w2*a22 - h1*B ;
B11 = w1*a31 + w3*a11 - w2*a21 - h2*B ;
C1  = A*B    +   w4    + w2*z2  - w3*z1 - w1*z3 ;

A22 = w2*a42 + w4*a22 - w1*a52 - w5*a12 - w3*a32 - h1*C;
B22 = w2*a41 + w4*a21 - w1*a51 - w5*a11 - w3*a31 - h2*C;
C2  = A*C    + w3*z3  - w6     + w5*z1  + w1*z5 - z4*w2 - w4*z2;

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Controller parameters
% --------------------------------
num_y_MO_dot = C2*A11  - A22*C1       ;
den_y_MO_dot = B22*A11 - A22*B11      ;
y_MO_dot  = num_y_MO_dot/ den_y_MO_dot;  

x_MO_dot  = (C1 - B11*y_MO_dot)/A11   ;

ti_MO = 2*kp*( A + h1*x_MO_dot + h2*y_MO_dot ); 

x_MO = y_MO_dot - 1;
y_MO = (1/2)*(1 + x_MO_dot - y_MO_dot);

num0Gc =        1                    ;
num1Gc =  r1*y_MO_dot - n1*x_MO_dot  ;
num2Gc =  r2*y_MO_dot - n2*x_MO_dot  ;
numGc_MO = [num2Gc num1Gc num0Gc]    ;
    
den0Gc = 0        ;
den1Gc = ti_MO*u1 ;
den2Gc = ti_MO*u2 ;
den3Gc = ti_MO*u3 ; 

denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
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
