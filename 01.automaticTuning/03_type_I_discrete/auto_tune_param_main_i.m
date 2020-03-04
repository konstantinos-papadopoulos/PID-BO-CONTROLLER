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

function [Gc Gp stepinformation ti_MO] = auto_tune_param_main_i(plant_loc)
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

% find tz_max
% -------------------------------------------------------------------------
tz_matr = [tz1 tz2 tz3 tz4];
tz_max = max(tz_matr)      ;


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

pGp2 = tp1*tp2 + tp1*tp3 + tp1*tp4 + tp1*tp5 + tp2*tp3 + ...
       tp2*tp4 + tp2*tp5 + tp3*tp4 + tp3*tp5 + tp4*tp5                           ;

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
    % Transfer Function of the Process with Time Delay
    % ---------------------------------------------------------------------  
    Gp = Gp1 * Gd              ; 
    n1 = td                    ;
    n2 = td^2 / 2              ;
    n3 = td^3 / 6              ;
else
    Gd = 1                     ;
    % Transfer Function of the Process without Time Delay
    % ---------------------------------------------------------------------  
    Gp = Gp1                   ; 
    n1 = 0                     ;
    n2 = 0                     ;
    n3 = td^3 / 6              ;
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



A1 = w1^2- w2 - w1*q1 + q2                              ;
B1 = q1 -  w1                                           ;
C1 = w1^3 - w1^2*q1 - 2*w1*w2 + w2*q1 + w3 + w1*q2 - q3 ;

A2 = w2^2 - w2*q2 + w4 + w3*q1 + w1*q3 - 2*w1*w3 - q4 ;
B2 = w3 - w2*q1 + w1*q2 - q3 ;
C2 = (w1 - q1)*(w2^2 - 2*w1*w3 + 2*w4) - (w1*q4 + w3*q2 + w5 - w4*q1 + w2*q3);

C = [C1;C2]       ;
D = [A1 B1;A2 B2] ;
E = inv(D)        ;
F = E*C           ; 

x_MO  = 0         ;
y_MO  = 0         ;  
ti_MO = 2*kp*(w1 - q1 - x_MO) ;
ti_MOInit = ti_MO ;

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
             
            
                
% Switching to I-LAg Controller...... because ti_MO < 0.
% -------------------------------------------------------------------------
flag = 0        ;
flag_optimal = 0;
if ti_MO <= 0
    flag = 1;
    
    if flag_optimal == 1
  
        disp('Switching to I-Lag Optimal Controller......')
        disp('********************************************************************')

        k1 = p1 + td         ;
        k2 = n2 + n1*p1 + p2 ;
        k3 = n2*p1 + n1*p2 + p3 + n3;

        f3 = 1;
        f2 = k1 - q1; 
        f1 = k1^2 - k1*q1 - k2 + q2;
        f0 = k1^3 - (k1^2)*q1 - 2*k1*k2 + q1*k2 + k1*q2 + k3 - q3;

        % solve the equation
        % ---------------------------------------------------------------------     
        % tx^3 + f2*tx^2 + f1*tx + f0 = 0
        Cmatr = [f3 f2 f1 f0];
        tx_MO = roots(Cmatr);

        tx_MO_1 = tx_MO(1); 
        tx_MO_2 = tx_MO(2);  
        tx_MO_3 = tx_MO(3);

        tx_MOFinal = max(tx_MO);
        tx_MO_comp = tz_max    ;
 
        tx_MO = tx_MOFinal;
        tx_MO;
        
        ti_MO = 2*kp*(tx_MO + w1 - q1)           ;
        ti_MO_comp = 2*kp*(tx_MO_comp + w1 - q1) ;
        
        % Optimal Lag Controller
        % ---------------------------------------------------------------------     
        x_MO = 0                                ;
        y_MO = 0                                ;
        num0Gc = 1                              ;
        num1Gc = x_MO                           ;
        num2Gc = y_MO                           ;
        numGc_MO = [num2Gc num1Gc num0Gc]       ;
        den0Gc = 0                              ;
        den1Gc = ti_MO                          ;
        den2Gc = ti_MO*(tp6 + tx_MO)            ;
        den3Gc = ti_MO*tp6*tx_MO                ;
        denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
        Gc_MO = tf(numGc_MO,denGc_MO)           ; 
        Gc = Gc_MO                              ;
        
        % Compensated Controller (Compensation between tx_MO and the greatest
        % zero of the plant) tx_MO_comp = tz_max    ;
        % -----------------------------------------------------------------     
        x_MO_comp = 0                                                    ;
        y_MO_comp = 0                                                    ;
        num0Gc_comp = 1                                                  ;
        num1Gc_comp = x_MO_comp                                          ;
        num2Gc_comp = y_MO_comp                                          ;
        numGc_MO_comp = [num2Gc_comp num1Gc_comp num0Gc_comp]            ;
        den0Gc_comp = 0                                                  ;
        den1Gc_comp = ti_MO_comp                                         ;
        den2Gc_comp = ti_MO_comp*(tp6 + tx_MO_comp)                      ;
        den3Gc_comp = ti_MO_comp*tp6*tx_MO_comp                          ;
        denGc_MO_comp = [den3Gc_comp den2Gc_comp   den1Gc_comp den0Gc_comp];
        Gc_MO_comp = tf(numGc_MO_comp,denGc_MO_comp)                     ;
    else
        disp('********************************************************************')
        disp('Switching to I-Lag Controller......')
        disp('********************************************************************')
        tx_MO = 0                               ;
        ti_MO = 2*kp*(tx_MO + w1 - q1)          ;
        %Controller Implementation
        % -----------------------------------------------------------------     
        x_MO = 0                                ;
        y_MO = 0                                ;
        num0Gc = 1                              ;
        num1Gc = x_MO                           ;
        num2Gc = y_MO                           ;
        numGc_MO = [num2Gc num1Gc num0Gc]       ;
        den0Gc = 0                              ;
        den1Gc = ti_MO                          ;
        den2Gc = ti_MO*(tp6 + tx_MO)            ;
        den3Gc = ti_MO*tp6*tx_MO                ;
        denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
        Gc_MO = tf(numGc_MO,denGc_MO)           ; 
        Gc = Gc_MO                              ;

        % Closing the Loop
        % -----------------------------------------------------------------     
        Ffp_MO = Gc_MO*Gp            ;
        Fcl_MO = feedback(Ffp_MO,kh) ;
         
        % measuring the overshoot
        % -----------------------------------------------------------------     
        overshoot_Fcl_MO_struct = stepinfo(Fcl_MO)          ;
        overshoot_Fcl_MO = overshoot_Fcl_MO_struct.Overshoot;
        ovrst_Lag = overshoot_Fcl_MO                        ;
        
        % Increasing tx_MO so that ti_MO increases until the ovs of the output
        % is equeal to 4.47%
        % -----------------------------------------------------------------     
        ovsUpperLimit = 4.475                        ;
        ovsLowerLimit = 4.4                          ;
        referenceOvershoot = 4.47                    ;
        errorRef = referenceOvershoot - ovrst_Lag    ;
        switchParam = 0.005                          ;
        stepNo_tx_MO = 0                             ;
        plotNo = 2                                   ;
        figureIndex = 1                              ;
        tx_MO = abs(ti_MO/2)                         ;
        
        while ((ovrst_Lag < ovsLowerLimit) || (ovrst_Lag > ovsUpperLimit)) || isnan(ovrst_Lag) == 1
            stepNo_tx_MO = stepNo_tx_MO + 1          ;
            if (errorRef > 0)
                if_index = 1;
                tx_MO = abs(tx_MO - (tx_MO*switchParam));
                % Updating ti_MO
                %****************************************
                ti_MO = 2*kp*(tx_MO + w1 - q1)          ;
                %Controller Implementation
                % ---------------------------------------------------------     
                x_MO = 0                                ;
                y_MO = 0                                ;
                num0Gc = 1                              ;
                num1Gc = x_MO                           ;
                num2Gc = y_MO                           ;
                numGc_MO = [num2Gc num1Gc num0Gc]       ;
                den0Gc = 0                              ;
                den1Gc = ti_MO                          ;
                den2Gc = ti_MO*(tp6 + tx_MO)            ;
                den3Gc = ti_MO*tp6*tx_MO                ;
                denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
                Gc_MO = tf(numGc_MO,denGc_MO)           ; 
                Gc = Gc_MO                              ;

                % Closing the Loop
                % ---------------------------------------------------------     
                Ffp_MO = Gc_MO*Gp            ;
                Fcl_MO = feedback(Ffp_MO,kh) ;

                % measuring the overshoot
                % ---------------------------------------------------------     
                overshoot_Fcl_MO_struct = stepinfo(Fcl_MO)          ;
                overshoot_Fcl_MO = overshoot_Fcl_MO_struct.Overshoot;
                ovrst_Lag = overshoot_Fcl_MO                        ;
                k = mod(stepNo_tx_MO,plotNo);
                      if (k == 0)
                           fprintf('stepNo_tx_MO: %d - tx_MO: %1.5f - ovrst_Lag: %1.5f - ti_MO: %1.5f - errorRef: %1.5f - if_index: %1.5f\n',stepNo_tx_MO,tx_MO,ovrst_Lag,ti_MO,errorRef,if_index)
                      end

                % Calculate the new errorRef
                % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                errorRef = referenceOvershoot - ovrst_Lag; 

            elseif (errorRef < 0)

                if_index = 2;
                tx_MO = abs(tx_MO + (tx_MO*switchParam)) ;
                % Updating ti_MO
                %****************************************
                ti_MO = 2*kp*(tx_MO + w1 - q1)          ;

                %Controller Implementation
                % ---------------------------------------------------------     
                x_MO = 0                                ;
                y_MO = 0                                ;
                num0Gc = 1                              ;
                num1Gc = x_MO                           ;
                num2Gc = y_MO                           ;
                numGc_MO = [num2Gc num1Gc num0Gc]       ;
                den0Gc = 0                              ;
                den1Gc = ti_MO                          ;
                den2Gc = ti_MO*(tp6 + tx_MO)            ;
                den3Gc = ti_MO*tp6*tx_MO                ;
                denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
                Gc_MO = tf(numGc_MO,denGc_MO)           ; 
                Gc = Gc_MO                              ;

                % Closing the Loop
                % ---------------------------------------------------------     
                Ffp_MO = Gc_MO*Gp            ;
                Fcl_MO = feedback(Ffp_MO,kh) ;

                % measuring the overshoot
                % ---------------------------------------------------------     
                overshoot_Fcl_MO_struct = stepinfo(Fcl_MO)          ;
                overshoot_Fcl_MO = overshoot_Fcl_MO_struct.Overshoot;
                ovrst_Lag = overshoot_Fcl_MO                        ;
                k = mod(stepNo_tx_MO,plotNo);
                      if (k == 0)
                           fprintf('stepNo_tx_MO: %d - tx_MO: %1.5f - ovrst_Lag: %1.5f - ti_MO: %1.5f - errorRef: %1.5f - if_index: %1.5f\n',stepNo_tx_MO,tx_MO,ovrst_Lag,ti_MO,errorRef,if_index)
                      end
                % Calculate the new errorRef
                % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                errorRef = referenceOvershoot - ovrst_Lag  ;      

            elseif isnan(errorRef) == 1
                if (ti_MO > 0) && (tx_MO > 0)
                    if_index = 4;
                    tx_MO = abs(tx_MO + (tx_MO*switchParam)) ;
                    % Updating ti_MO
                    %****************************************
                    ti_MO = 2*kp*(tx_MO + w1 - q1)          ;

                    %Controller Implementation
                    % ---------------------------------------------------------     
                    x_MO = 0                                ;
                    y_MO = 0                                ;
                    num0Gc = 1                              ;
                    num1Gc = x_MO                           ;
                    num2Gc = y_MO                           ;
                    numGc_MO = [num2Gc num1Gc num0Gc]       ;
                    den0Gc = 0                              ;
                    den1Gc = ti_MO                          ;
                    den2Gc = ti_MO*(tp6 + tx_MO)            ;
                    den3Gc = ti_MO*tp6*tx_MO                ;
                    denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
                    Gc_MO = tf(numGc_MO,denGc_MO)           ; 
                    Gc = Gc_MO                              ;

                    % Closing the Loop
                    % ---------------------------------------------------------     
                    Ffp_MO = Gc_MO*Gp            ;
                    Fcl_MO = feedback(Ffp_MO,kh) ;

                    % measuring the overshoot
                    % ---------------------------------------------------------     
                    overshoot_Fcl_MO_struct = stepinfo(Fcl_MO)          ;
                    overshoot_Fcl_MO = overshoot_Fcl_MO_struct.Overshoot;
                    ovrst_Lag = overshoot_Fcl_MO                        ;
                    k = mod(stepNo_tx_MO,plotNo);
                          if (k == 0)
                               fprintf('stepNo_tx_MO: %d - tx_MO: %1.5f - ovrst_Lag: %1.5f - ti_MO: %1.5f - errorRef: %1.5f - if_index: %1.5f\n',stepNo_tx_MO,tx_MO,ovrst_Lag,ti_MO,errorRef,if_index)
                          end

                    % Calculate the new errorRef
                    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    errorRef = referenceOvershoot - ovrst_Lag  ;      

                elseif (ti_MO < 0) && (tx_MO > 0)
                    if_index = 5;
                    tx_MO = abs(tx_MO + (tx_MO*switchParam)) ;
                    % Updating ti_MO
                    %****************************************
                    ti_MO = 2*kp*(tx_MO + w1 - q1)          ;

                    %Controller Implementation
                    % ---------------------------------------------------------     
                    x_MO = 0                                ;
                    y_MO = 0                                ;
                    num0Gc = 1                              ;
                    num1Gc = x_MO                           ;
                    num2Gc = y_MO                           ;
                    numGc_MO = [num2Gc num1Gc num0Gc]       ;
                    den0Gc = 0                              ;
                    den1Gc = ti_MO                          ;
                    den2Gc = ti_MO*(tp6 + tx_MO)            ;
                    den3Gc = ti_MO*tp6*tx_MO                ;
                    denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
                    Gc_MO = tf(numGc_MO,denGc_MO)           ; 
                    Gc = Gc_MO                              ;

                    % Closing the Loop
                    % ---------------------------------------------------------     
                    Ffp_MO = Gc_MO*Gp            ;
                    Fcl_MO = feedback(Ffp_MO,kh) ;

                    % measuring the overshoot
                    % ---------------------------------------------------------     
                    overshoot_Fcl_MO_struct = stepinfo(Fcl_MO)          ;
                    overshoot_Fcl_MO = overshoot_Fcl_MO_struct.Overshoot;
                    ovrst_Lag = overshoot_Fcl_MO                        ;
                    k = mod(stepNo_tx_MO,plotNo);
                          if (k == 0)
                               fprintf('stepNo_tx_MO: %d - tx_MO: %1.5f - ovrst_Lag: %1.5f - ti_MO: %1.5f - errorRef: %1.5f - if_index: %1.5f\n',stepNo_tx_MO,tx_MO,ovrst_Lag,ti_MO,errorRef,if_index)
                          end
                    % Calculate the new errorRef
                    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    errorRef = referenceOvershoot - ovrst_Lag  ;      

                else %(ti_MO < 0) && (tx_MO < 0)
                    if_index = 5;
                    tx_MO = abs(tx_MO + (tx_MO*switchParam)) ;

                    % Updating ti_MO
                    %****************************************
                    ti_MO = 2*kp*(tx_MO + w1 - q1)          ;

                    %Controller Implementation
                    % ---------------------------------------------------------     
                    x_MO = 0                                ;
                    y_MO = 0                                ;
                    num0Gc = 1                              ;
                    num1Gc = x_MO                           ;
                    num2Gc = y_MO                           ;
                    numGc_MO = [num2Gc num1Gc num0Gc]       ;
                    den0Gc = 0                              ;
                    den1Gc = ti_MO                          ;
                    den2Gc = ti_MO*(tp6 + tx_MO)            ;
                    den3Gc = ti_MO*tp6*tx_MO                ;
                    denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
                    Gc_MO = tf(numGc_MO,denGc_MO)           ; 
                    Gc = Gc_MO                              ;

                    % Closing the Loop
                    % ---------------------------------------------------------     
                    Ffp_MO = Gc_MO*Gp            ;
                    Fcl_MO = feedback(Ffp_MO,kh) ;

                    % measuring the overshoot
                    % ---------------------------------------------------------     
                    overshoot_Fcl_MO_struct = stepinfo(Fcl_MO)          ;
                    overshoot_Fcl_MO = overshoot_Fcl_MO_struct.Overshoot;
                    ovrst_Lag = overshoot_Fcl_MO                        ;
                    k = mod(stepNo_tx_MO,plotNo);
                          if (k == 0)
                               fprintf('stepNo_tx_MO: %d - tx_MO: %1.5f - ovrst_Lag: %1.5f - ti_MO: %1.5f - errorRef: %1.5f - if_index: %1.5f\n',stepNo_tx_MO,tx_MO,ovrst_Lag,ti_MO,errorRef,if_index)
                          end
                    % Calculate the new errorRef
                    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    errorRef = referenceOvershoot - ovrst_Lag  ;      
                end

            end %of IF - line:269

        end %of WHILE  - line:267
        
        
    end %of IF - line:165
    
  
end %of IF - line:162

% -------------------------------------------------------------------------
% Forward Path and Closed Loop Implementation
% -------------------------------------------------------------------------  
% if flag == 1
%     ltiview(Fcl_MO,'r');
%     ltiview(Gp,'k');
% else
%     ltiview(Fcl_MO,'b');
% end



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
