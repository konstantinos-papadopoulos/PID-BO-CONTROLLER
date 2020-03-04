% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_main_pid_typeII.m																			   																		  	
%  Project:     Automatic tuning of the paramters for PI,PID controllers - 
%  
%  Purpose:     calculate "pi" controller parameters type II systems
%  
%  Author :     kostas g. papadopoulos																	   																		
% 																										   																		
%  History:     Date: 10.07.2008  date last modified
% 																										  																		
%  Contact:     kostas g. papadopoulos,    konstantina mermikli
%               kpapadop@eng.auth.gr  ,    kmermikl@auth.gr ,      
% 																										  																		
%  Place:	    Aristotle University of Thessaloniki, Thessaloniki, Greece							   																		
% 
% ************************************************************************************************************

function [Gc_meth_2 Gp stepinformation_meth_2 ti_sq_meth_2 x_meth2_sol y_meth2_sol ...
          Gc_meth_1 stepinformation_meth_1 ti_sq_meth_1 x_meth1_sol y_meth1_sol ] = auto_tune_param_main_pid_typeII(plant_loc)

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
q5 = 0;
q6 = 0; 
% z0 = q0;
% z1 = q1;
% z2 = q2;
% z3 = q3;
% z4 = q4;
% -------------------------------------------------------------------------

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
    Gd = tf(num_d,den_d)       ;      
    Gp = Gp1 * Gd              ;      % Transfer Function of the Process with Time Delay
    n1 = td                    ;
    n2 = td^2 / 2              ;
    n3 = td^3 / 6              ;   
else
    Gd = 1                     ;
    Gp = Gp1                   ;        % Transfer Function of the Process without Time Delay
    n1 = 0                     ;
    n2 = 0                     ;
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
w2 = n2 + n1*p1 + p2    ;
w3 = n2*p1 + n1*p2 + p3 ;
w4 = n2*p2 + n1*p3 + p4 ;
w5 = n2*p3 + n1*p4 + p5 ;
w6 = n2*p4 + n1*p5 + p6 ;
w7 = n2*p5 + n1*p6      ;
w8 = n2*p6              ;

%dikes mou prakseis NEW
Q0 = w1^2 - 2*w2;
Q1 = w1 - q1                        ;
Q2 = w2 + q2 - w1*q1                ;
Q3 = w1*q2 - w2*q1 + w3 - q3        ;
Q4 = w1*q3 + w3*q1 - w2*q2 - w4 - q4;

Z1 = q1^2 - 2*q2                    ;
Z2 = 0                              ;
Z3 = q2^2 + 2*q4 - 2*q1*q3          ;
Z4 = q3^2 + 2*q1*q5 - 2*q2*q4 - 2*q6;
Z5 = 0                              ;

W1 = 4*Z3 - 3*(Z1^2)                ;
W2 = 2*Q4 - Z1*Q2 - (Z1^2) + 2*Z3   ;
W3 = Z1^3 + 4*Z4 - 4*Z1*Z3          ;
W4 = 0                              ;



% 1st Method (...leading to 4 roots): eq.47
% -------------------------------------------------------------------------
% x^8 : C(8) = 4*(Q0+2*Q2+Z1);
% x^7 : C(7) = -32*Q1*(Z1+Q2);
% x^6 : C(6) = 16*(Z1*(Q0+2*Q2+Z1)+2*((Z1+Q2)^2+2*Q1*Q3)+4*(Q1^2)*Z1);
% x^5 : C(5) = -32*(2*(Q2+Z1)*(Z1*Q1+Q3)-Q1*W2+4*Z1*Q1*(Z1+Q2));
% x^4 : C(4) = 8*((Q0+2*Q2+Z1)*(2*(Z1^2)+W1)+4*(2*Z1*(Z1+Q2)^2+4*Q1*Q3*Z1-W2*(Z1+Q2))+8*(Z1*(Z1+Q2)^2+2*(Z1^2)*(Q1^2)));
% x^3 : C(3) = -32*((Z1+Q2)*(Q1*W1+4*Q3*Z1)-2*Q1*Z1*W2+8*Q1*(Z1^2)*(Z1+Q2));
% x^2 : C(2) = 16*(Z1*W1*(Q0+2*Q2+Z1)+2*(W1*(Z1+Q2)^2+2*Q1*Q3*W1-2*Z1*W2*(Z1+Q2))+4*(2*(Z1^2)*(Z1+Q2)^2+W3*(Q1^2)));
% x^1 : C(1) = -32*(2*Q3*W1*(Z1+Q2)-Q1*W1*W2+4*W3*Q1*(Z1+Q2));
% x^0 : C(0) = 4*((W1^2)*(Q0+2*Q2+Z1)-8*W1*W2*(Z1+Q2)+16*W3*(Z1+Q2)^2);

A8 = 4*(Q0+2*Q2+Z1);
A7 = -32*Q1*(Z1+Q2);
A6 = 16*(Z1*(Q0+2*Q2+Z1)+2*((Z1+Q2)^2+2*Q1*Q3)+4*(Q1^2)*Z1);
A5 = -32*(2*(Q2+Z1)*(Z1*Q1+Q3)-Q1*W2+4*Z1*Q1*(Z1+Q2));
A4 = 8*((Q0+2*Q2+Z1)*(2*(Z1^2)+W1)+4*(2*Z1*(Z1+Q2)^2+4*Q1*Q3*Z1-W2*(Z1+Q2))+8*(Z1*(Z1+Q2)^2+2*(Z1^2)*(Q1^2)));
A3 = -32*((Z1+Q2)*(Q1*W1+4*Q3*Z1)-2*Q1*Z1*W2+8*Q1*(Z1^2)*(Z1+Q2));
A2 = 16*(Z1*W1*(Q0+2*Q2+Z1)+2*(W1*(Z1+Q2)^2+2*Q1*Q3*W1-2*Z1*W2*(Z1+Q2))+4*(2*(Z1^2)*(Z1+Q2)^2+W3*(Q1^2)));
A1 = -32*(2*Q3*W1*(Z1+Q2)-Q1*W1*W2+4*W3*Q1*(Z1+Q2));
A0 = 4*((W1^2)*(Q0+2*Q2+Z1)-8*W1*W2*(Z1+Q2)+16*W3*(Z1+Q2)^2);


% c_meth_1(1) = 4*(Q0 + 2*Q2 + Z1)                                                  ;  %x^8 
% c_meth_1(2) = -32*Q1*(Z1 + Q2)                                                    ;  %x^7 
% c_meth_1(3) = 16*(Z1*(Q0 + 2*Q2 + Z1) + 2*((Z1 + Q2)^2 + 2*Q1*Q3) + 4*(Q1^2)*Z1)  ;  %x^6 
% c_meth_1(4) = -32*(2*(Q2 + Z1)*(Z1*Q1 + Q3) - Q1*W2 + 4*Z1*Q1*(Z1 + Q2))          ;  %x^5 
% c_meth_1(5) = 8*((Q0+2*Q2+ Z1)*(2*(Z1^2) + W1) + 4*(2*Z1*(Z1 + Q2)^2 + ...
%               4*Q1*Q3*Z1 - W2*(Z1 + Q2)) + 8*(Z1*(Z1 + Q2)^2 + 2*(Z1^2)*(Q1^2)))  ;  %x^4 
% c_meth_1(6) = -32*((Z1+Q2)*(Q1*W1+4*Q3*Z1)-2*Q1*Z1*W2+8*Q1*(Z1^2)*(Z1+Q2))        ;  %x^3 
% c_meth_1(7) = 16*(Z1*W1*(Q0 + 2*Q2 + Z1) + 2*(W1*(Z1 + Q2)^2 + 2*Q1*Q3*W1 - ...
%               2*Z1*W2*(Z1 + Q2)) + 4*(2*(Z1^2)*(Z1 + Q2)^2 + W3*(Q1^2)))          ;  %x^2 
% c_meth_1(8) = -32*(2*Q3*W1*(Z1 + Q2) - Q1*W1*W2 + 4*W3*Q1*(Z1 + Q2))              ;  %x^1 
% c_meth_1(9) = 4*((W1^2)*(Q0+2*Q2+Z1)-8*W1*W2*(Z1+Q2)+16*W3*(Z1+Q2)^2)             ;  %x^0 
                      
% -----------------------------------------------ss
% C(1)*X^N + ... + C(N)*X + C(N+1) [roots(C)]   |
% -----------------------------------------------
% 1st Method (...leading to 8 roots): eq.54
% -------------------------------------------------------------------------
% candidate real solutions
% -------------------------------------------------------------------------
c = [A8 A7 A6 A5 A4 A3 A2 A1 A0];
p = roots(c);
% p = roots(c_meth_1);
j = 1;
% collecting candidate solutions
% **********************************
for i = 1:8
    if (isreal(p(i)) && p(i) > 0)
        k(j) = p(i);
        j = j + 1;
    end
end
if j == 1
    x_MO = -10;
else
    x_MO = max(k);
end
R = (x_MO^4 + 2*Z1*(x_MO^2) + W1)/(8*(Q1*x_MO - Z1 - Q2));
y_MO = (1/2)*(x_MO^2 + Z1)-R;
% ++++++++++++++++++++++++++++++++++++++++++++++
    x_meth1_sol = x_MO;
    y_meth1_sol = y_MO;
    ti_sq_meth_1 = 0.5*kp*(x_meth1_sol^2 - 2*y_meth1_sol + q1^2 - 2*q2) ;
% ++++++++++++++++++++++++++++++++++++++++++++++

% 2nd Method (...leading to 2 roots): eq.54
% -------------------------------------------------------------------------
% x^2 : 2*(w1*(w1-q1)-w2+q2)
% x^1 : -4*(w1^3-3*q1*w1^2+2*w1*q1^2+w1*q2+w2*q1-w3+q3-2*q1*q2);
% x^0 : (w1^2-2*w1*q1+2*q2)*(q1^2+2*q2+4*w2-4*w1*q1)+(w1^2-2*w2)*(q1^2-2*q2)+4*(w1*q3+w3*q1-w4-q4-w2*q2);

c_meth_2(1) = 2*(w1*(w1 - q1) - w2 + q2)                                       ;
c_meth_2(2) = -4*(w1^3 - 3*q1*w1^2 + 2*w1*q1^2 +...
              w1*q2 + w2*q1 - w3 + q3 - 2*q1*q2)                               ;
c_meth_2(3) = (w1^2 - 2*w1*q1 + 2*q2)*(q1^2 + 2*q2 + 4*w2 - 4*w1*q1) + ...
              (w1^2 - 2*w2)*(q1^2 - 2*q2) + 4*(w1*q3 + w3*q1 - w4 - q4 - w2*q2);
 
r1 = c_meth_2(1);
r2 = c_meth_2(2);
r3 = c_meth_2(3);

x_meth2_1 = -r2/(2*r1) + sqrt(r2^2 - 4*r1*r3)/(2*r1);          

if (isreal(x_meth2_1) && x_meth2_1>0)
    x_meth2_1 = x_meth2_1;
    x = x_meth2_1;
else
    x_meth2_1 = -10;
    x = x_meth2_1  ;
end

% x_meth2 = roots(c_meth_2)        ;                            
% if isreal(x_meth2) == 1
%     x_meth2_1 = x_meth2(1)       ;
%     x_meth2_2 = x_meth2(2)       ;
%     x = max(x_meth2_1,x_meth2_2) ;
% else
%     disp('...complex solution for tn: METHOD 2 - 2:roots')
%     x_meth2_2 = -10;
%     x = x_meth2_2  ;
% end
% ++++++++++++++++++++++++++++++++++++++++++++++
x_meth2_sol = x;
% -------------------------------------------------------------------------
% y = -0.5*(x^2) + 2*(w1 - q1)*x - 0.5*(q1^2 + 2*q2 + 4*w2 - 4*w1*q1)
% -------------------------------------------------------------------------
y_MO = -0.5*(x_meth2_sol^2) + 2*(w1 - q1)*x_meth2_sol - 0.5*(q1^2 + 2*q2 + 4*w2 - 4*w1*q1);
y_meth2_sol = y_MO;
% ++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
% ti_sq = 0.5*kp(x^2 - 2*y + z1^2 - 2*z2)
% -------------------------------------------------------------------------
 ti_sq_meth_2 = 0.5*kp*(x_meth2_sol^2 - 2*y_meth2_sol + q1^2 - 2*q2) ;

% 1st Method - Controller
% -------------------------------------------------------------------------
    
num0Gc_meth_1 = 1                       ;
num1Gc_meth_1 = x_meth1_sol             ;
num2Gc_meth_1 = y_meth1_sol             ;
numGc_MO_meth_1 = [num2Gc_meth_1 num1Gc_meth_1 num0Gc_meth_1];
    
den0Gc_meth_1 = 0                       ;
den1Gc_meth_1 = 0                       ;
den2Gc_meth_1 = ti_sq_meth_1            ;
den3Gc_meth_1 = ti_sq_meth_1*tp6        ;

denGc_MO_meth_1 = [den3Gc_meth_1 den2Gc_meth_1 den1Gc_meth_1 den0Gc_meth_1];
Gc_MO_meth_1 = tf(numGc_MO_meth_1,denGc_MO_meth_1); 
Gc_meth_1 = Gc_MO_meth_1                          ;

% -------------------------------------------------------------------------
% Forward Path and Closed Loop Implementation
% -------------------------------------------------------------------------  
Ffp_MO_meth1 = Gc_meth_1*Gp                    ;
Fcl_MO_meth1 = feedback(Ffp_MO_meth1,kh)       ;
stepinformation_meth_1 = stepinfo(Fcl_MO_meth1);

% 2nd Method - Controller
% -------------------------------------------------------------------------
    
num0Gc_meth_2 = 1                       ;
num1Gc_meth_2 = x_meth2_sol             ;
num2Gc_meth_2 = y_meth2_sol             ;
numGc_MO_meth_2 = [num2Gc_meth_2 num1Gc_meth_2 num0Gc_meth_2];
    
den0Gc_meth_2 = 0                       ;
den1Gc_meth_2 = 0                       ;
den2Gc_meth_2 = ti_sq_meth_2            ;
den3Gc_meth_2 = ti_sq_meth_2*tp6        ;

denGc_MO_meth_2 = [den3Gc_meth_2 den2Gc_meth_2 den1Gc_meth_2 den0Gc_meth_2];
Gc_MO_meth_2 = tf(numGc_MO_meth_2,denGc_MO_meth_2); 
Gc_meth_2 = Gc_MO_meth_2                          ;

% -------------------------------------------------------------------------
% Forward Path and Closed Loop Implementation
% -------------------------------------------------------------------------  
Ffp_MO_meth2 = Gc_meth_2*Gp                    ;
Fcl_MO_meth2 = feedback(Ffp_MO_meth2,kh)       ;
stepinformation_meth_2 = stepinfo(Fcl_MO_meth2);


% figure(2)
% step(Fcl_MO);
% ltiview(Fcl_MO_meth1,'b',Fcl_MO_meth2,'r')
% step(Fcl_MO_meth2)
% ltiview(Fcl_MO_meth2,'r')
