% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_calculateOvs.m																			   																		  	
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

function [ovrst Fcl ti tnx tvx] = auto_tune_param_calculateOvs(plant_loc)

tp1 = plant_loc.tp1;    tz1 = plant_loc.tz1; 
tp2 = plant_loc.tp2;    tz2 = plant_loc.tz2;
tp3 = plant_loc.tp3;    tz3 = plant_loc.tz3;
tp4 = plant_loc.tp4;    tz4 = plant_loc.tz4;
tp5 = plant_loc.tp5;   
tp6 = plant_loc.tp6;   
td = plant_loc.td  ;    kp = plant_loc.kp;      kh = plant_loc.kh;
tsx = plant_loc.tsx;    tnx = plant_loc.tnx;    tvx = plant_loc.tvx;

% Formulate plant transfer function [numerator polynomial coefficients]
% -------------------------------------------------------------------------
q1 = tz1 + tz2 + tz3 + tz4;
q2 = tz1*tz2 + tz1*tz3 + tz1*tz4 + tz2*tz3 + tz2*tz4 + tz3*tz4;
q3 = tz1*tz2*tz3 + tz1*tz2*tz4 + tz2*tz3*tz4 + tz1*tz3*tz4;
q4 = tz1*tz2*tz3*tz4;

% Formulate plant transfer function [denominator polynomial coefficients]
% -------------------------------------------------------------------------
pGp1 = tp1 + tp2 + tp3 + tp4 + tp5;
pGp2 = tp1*tp2 + tp1*tp3 + tp1*tp4 + tp1*tp5 + tp2*tp3 + ...
       tp2*tp4 + tp2*tp5 + tp3*tp4 + tp3*tp5 + tp4*tp5      ;
pGp3 = tp1*tp2*tp3 + tp1*tp2*tp4 + tp1*tp2*tp5 + tp1*tp3*tp4 + tp1*tp3*tp5 +...
       tp1*tp4*tp5 + tp2*tp3*tp4 + tp2*tp3*tp5 + tp2*tp4*tp5 + tp3*tp4*tp5;
pGp4 = tp1*tp2*tp3*tp4 + tp1*tp2*tp3*tp5 + tp1*tp2*tp4*tp5 + tp1*tp3*tp4*tp5 +...
       tp2*tp3*tp4*tp5;
pGp5 = tp1*tp2*tp3*tp4*tp5;
        
numGp = kp * [q4 q3 q2 q1 1]        ;
denGp = [pGp5 pGp4 pGp3 pGp2 pGp1 1];
Gp1 = tf(numGp,denGp)               ;

% Formulate time delay transfer function using PADE approximation
% -------------------------------------------------------------------------
if (td ~= 0)                          
    N = 20                     ;
    [num_d,den_d] = pade(td,N) ;
    Gd = tf(num_d,den_d)       ;
else
    Gd = 1   ;
end
Gp = Gp1 * Gd;

% Formulate controller transfer function with the updated tsx, tnx, tvx
% -------------------------------------------------------------------------
numGc = [(tnx*tvx) (tnx + tvx) 1];    
denGc = [(2*kp*(tsx - tnx - tvx)*tp6) (2*kp*(tsx - tnx - tvx)) 0];
Gc = tf(numGc,denGc) ;

Ffp = Gc * Gp          ;  % forward path transfer function
Fcl = feedback(Ffp,kh) ;  % closed loop trnasfer function  

S = stepinfo(Fcl)      ;
ovrst = S.Overshoot    ;
ti = 2*kp*(tsx - tnx - tvx);

% 
% 
% % Switching to I-Lag Controller
% % -------------------------------------------------------------------------
% 
% disp('********************************************************************')
% disp('Switching to I-Lag Controller......')
% disp('********************************************************************')
% tx_MO = 0                   ;
% ti_MO = 2*kp*(tx_MO)        ;
% 
% %Controller Implementation
% % -----------------------------------------------------------------     
% x_MO = 0                                ;
% y_MO = 0                                ;
% num0Gc = 1                              ;
% num1Gc = x_MO                           ;
% num2Gc = y_MO                           ;
% numGc_MO = [num2Gc num1Gc num0Gc]       ;
% den0Gc = 0                              ;
% den1Gc = ti_MO                          ;
% den2Gc = ti_MO*(tp6 + tx_MO)            ;
% den3Gc = ti_MO*tp6*tx_MO                ;
% denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
% Gc_MO = tf(numGc_MO,denGc_MO)           ; 
% Gc = Gc_MO                              ;
% 
% % Closing the Loop
% % -----------------------------------------------------------------     
% Ffp_MO = Gc_MO*Gp            ;
% Fcl_MO = feedback(Ffp_MO,kh) ;
%          
% % measuring the overshoot
% % -----------------------------------------------------------------     
% overshoot_Fcl_MO_struct = stepinfo(Fcl_MO)          ;
% overshoot_Fcl_MO = overshoot_Fcl_MO_struct.Overshoot;
% ovrst_Lag = overshoot_Fcl_MO                        ;
%         
% % Increasing tx_MO so that ti_MO increases until the ovs of the output
% % is equeal to 4.47%
% % -----------------------------------------------------------------     
%         ovsUpperLimit = 4.475                        ;
%         ovsLowerLimit = 4.4                          ;
%         referenceOvershoot = 4.47                    ;
%         errorRef = referenceOvershoot - ovrst_Lag    ;
%         switchParam = 0.005                          ;
%         stepNo_tx_MO = 0                             ;
%         plotNo = 2                                   ;
%         figureIndex = 1                              ;
% 
%         % Setting the Initial Condition for the tx_MO
%         % ...as if it started with I Control
%         % -----------------------------------------------------------------     
% 
%         tx_MO = abs(ti_MO/2)                         ;
%         
%         while ((ovrst_Lag < ovsLowerLimit) || (ovrst_Lag > ovsUpperLimit)) || isnan(ovrst_Lag) == 1
%             stepNo_tx_MO = stepNo_tx_MO + 1          ;
%             if (errorRef > 0)
%                 if_index = 1;
%                 tx_MO = abs(tx_MO - (tx_MO*switchParam));
%                 % Updating ti_MO
%                 %****************************************
%                 ti_MO = 2*kp*(tx_MO + w1 - q1)          ;
%                 %Controller Implementation
%                 % ---------------------------------------------------------     
%                 x_MO = 0                                ;
%                 y_MO = 0                                ;
%                 num0Gc = 1                              ;
%                 num1Gc = x_MO                           ;
%                 num2Gc = y_MO                           ;
%                 numGc_MO = [num2Gc num1Gc num0Gc]       ;
%                 den0Gc = 0                              ;
%                 den1Gc = ti_MO                          ;
%                 den2Gc = ti_MO*(tp6 + tx_MO)            ;
%                 den3Gc = ti_MO*tp6*tx_MO                ;
%                 denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
%                 Gc_MO = tf(numGc_MO,denGc_MO)           ; 
%                 Gc = Gc_MO                              ;
% 
%                 % Closing the Loop
%                 % ---------------------------------------------------------     
%                 Ffp_MO = Gc_MO*Gp            ;
%                 Fcl_MO = feedback(Ffp_MO,kh) ;
% 
%                 % measuring the overshoot
%                 % ---------------------------------------------------------     
%                 overshoot_Fcl_MO_struct = stepinfo(Fcl_MO)          ;
%                 overshoot_Fcl_MO = overshoot_Fcl_MO_struct.Overshoot;
%                 ovrst_Lag = overshoot_Fcl_MO                        ;
%                 k = mod(stepNo_tx_MO,plotNo);
%                       if (k == 0)
%                            fprintf('stepNo_tx_MO: %d - tx_MO: %1.5f - ovrst_Lag: %1.5f - ti_MO: %1.5f - errorRef: %1.5f - if_index: %1.5f\n',stepNo_tx_MO,tx_MO,ovrst_Lag,ti_MO,errorRef,if_index)
%                       end
% 
%                 % Calculate the new errorRef
%                 % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                 errorRef = referenceOvershoot - ovrst_Lag; 
% 
%             elseif (errorRef < 0)
% 
%                 if_index = 2;
%                 tx_MO = abs(tx_MO + (tx_MO*switchParam)) ;
%                 % Updating ti_MO
%                 %****************************************
%                 ti_MO = 2*kp*(tx_MO + w1 - q1)          ;
% 
%                 %Controller Implementation
%                 % ---------------------------------------------------------     
%                 x_MO = 0                                ;
%                 y_MO = 0                                ;
%                 num0Gc = 1                              ;
%                 num1Gc = x_MO                           ;
%                 num2Gc = y_MO                           ;
%                 numGc_MO = [num2Gc num1Gc num0Gc]       ;
%                 den0Gc = 0                              ;
%                 den1Gc = ti_MO                          ;
%                 den2Gc = ti_MO*(tp6 + tx_MO)            ;
%                 den3Gc = ti_MO*tp6*tx_MO                ;
%                 denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
%                 Gc_MO = tf(numGc_MO,denGc_MO)           ; 
%                 Gc = Gc_MO                              ;
% 
%                 % Closing the Loop
%                 % ---------------------------------------------------------     
%                 Ffp_MO = Gc_MO*Gp            ;
%                 Fcl_MO = feedback(Ffp_MO,kh) ;
% 
%                 % measuring the overshoot
%                 % ---------------------------------------------------------     
%                 overshoot_Fcl_MO_struct = stepinfo(Fcl_MO)          ;
%                 overshoot_Fcl_MO = overshoot_Fcl_MO_struct.Overshoot;
%                 ovrst_Lag = overshoot_Fcl_MO                        ;
%                 k = mod(stepNo_tx_MO,plotNo);
%                       if (k == 0)
%                            fprintf('stepNo_tx_MO: %d - tx_MO: %1.5f - ovrst_Lag: %1.5f - ti_MO: %1.5f - errorRef: %1.5f - if_index: %1.5f\n',stepNo_tx_MO,tx_MO,ovrst_Lag,ti_MO,errorRef,if_index)
%                       end
%                 % Calculate the new errorRef
%                 % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                 errorRef = referenceOvershoot - ovrst_Lag  ;      
% 
%             elseif isnan(errorRef) == 1
%                 if (ti_MO > 0) && (tx_MO > 0)
%                     if_index = 4;
%                     tx_MO = abs(tx_MO + (tx_MO*switchParam)) ;
%                     % Updating ti_MO
%                     %****************************************
%                     ti_MO = 2*kp*(tx_MO + w1 - q1)          ;
% 
%                     %Controller Implementation
%                     % ---------------------------------------------------------     
%                     x_MO = 0                                ;
%                     y_MO = 0                                ;
%                     num0Gc = 1                              ;
%                     num1Gc = x_MO                           ;
%                     num2Gc = y_MO                           ;
%                     numGc_MO = [num2Gc num1Gc num0Gc]       ;
%                     den0Gc = 0                              ;
%                     den1Gc = ti_MO                          ;
%                     den2Gc = ti_MO*(tp6 + tx_MO)            ;
%                     den3Gc = ti_MO*tp6*tx_MO                ;
%                     denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
%                     Gc_MO = tf(numGc_MO,denGc_MO)           ; 
%                     Gc = Gc_MO                              ;
% 
%                     % Closing the Loop
%                     % ---------------------------------------------------------     
%                     Ffp_MO = Gc_MO*Gp            ;
%                     Fcl_MO = feedback(Ffp_MO,kh) ;
% 
%                     % measuring the overshoot
%                     % ---------------------------------------------------------     
%                     overshoot_Fcl_MO_struct = stepinfo(Fcl_MO)          ;
%                     overshoot_Fcl_MO = overshoot_Fcl_MO_struct.Overshoot;
%                     ovrst_Lag = overshoot_Fcl_MO                        ;
%                     k = mod(stepNo_tx_MO,plotNo);
%                           if (k == 0)
%                                fprintf('stepNo_tx_MO: %d - tx_MO: %1.5f - ovrst_Lag: %1.5f - ti_MO: %1.5f - errorRef: %1.5f - if_index: %1.5f\n',stepNo_tx_MO,tx_MO,ovrst_Lag,ti_MO,errorRef,if_index)
%                           end
% 
%                     % Calculate the new errorRef
%                     % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                     errorRef = referenceOvershoot - ovrst_Lag  ;      
% 
%                 elseif (ti_MO < 0) && (tx_MO > 0)
%                     if_index = 5;
%                     tx_MO = abs(tx_MO + (tx_MO*switchParam)) ;
%                     % Updating ti_MO
%                     %****************************************
%                     ti_MO = 2*kp*(tx_MO + w1 - q1)          ;
% 
%                     %Controller Implementation
%                     % ---------------------------------------------------------     
%                     x_MO = 0                                ;
%                     y_MO = 0                                ;
%                     num0Gc = 1                              ;
%                     num1Gc = x_MO                           ;
%                     num2Gc = y_MO                           ;
%                     numGc_MO = [num2Gc num1Gc num0Gc]       ;
%                     den0Gc = 0                              ;
%                     den1Gc = ti_MO                          ;
%                     den2Gc = ti_MO*(tp6 + tx_MO)            ;
%                     den3Gc = ti_MO*tp6*tx_MO                ;
%                     denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
%                     Gc_MO = tf(numGc_MO,denGc_MO)           ; 
%                     Gc = Gc_MO                              ;
% 
%                     % Closing the Loop
%                     % ---------------------------------------------------------     
%                     Ffp_MO = Gc_MO*Gp            ;
%                     Fcl_MO = feedback(Ffp_MO,kh) ;
% 
%                     % measuring the overshoot
%                     % ---------------------------------------------------------     
%                     overshoot_Fcl_MO_struct = stepinfo(Fcl_MO)          ;
%                     overshoot_Fcl_MO = overshoot_Fcl_MO_struct.Overshoot;
%                     ovrst_Lag = overshoot_Fcl_MO                        ;
%                     k = mod(stepNo_tx_MO,plotNo);
%                           if (k == 0)
%                                fprintf('stepNo_tx_MO: %d - tx_MO: %1.5f - ovrst_Lag: %1.5f - ti_MO: %1.5f - errorRef: %1.5f - if_index: %1.5f\n',stepNo_tx_MO,tx_MO,ovrst_Lag,ti_MO,errorRef,if_index)
%                           end
%                     % Calculate the new errorRef
%                     % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                     errorRef = referenceOvershoot - ovrst_Lag  ;      
% 
%                 else %(ti_MO < 0) && (tx_MO < 0)
%                     if_index = 5;
%                     tx_MO = abs(tx_MO + (tx_MO*switchParam)) ;
% 
%                     % Updating ti_MO
%                     %****************************************
%                     ti_MO = 2*kp*(tx_MO + w1 - q1)          ;
% 
%                     %Controller Implementation
%                     % ---------------------------------------------------------     
%                     x_MO = 0                                ;
%                     y_MO = 0                                ;
%                     num0Gc = 1                              ;
%                     num1Gc = x_MO                           ;
%                     num2Gc = y_MO                           ;
%                     numGc_MO = [num2Gc num1Gc num0Gc]       ;
%                     den0Gc = 0                              ;
%                     den1Gc = ti_MO                          ;
%                     den2Gc = ti_MO*(tp6 + tx_MO)            ;
%                     den3Gc = ti_MO*tp6*tx_MO                ;
%                     denGc_MO = [den3Gc den2Gc den1Gc den0Gc];
%                     Gc_MO = tf(numGc_MO,denGc_MO)           ; 
%                     Gc = Gc_MO                              ;
% 
%                     % Closing the Loop
%                     % ---------------------------------------------------------     
%                     Ffp_MO = Gc_MO*Gp            ;
%                     Fcl_MO = feedback(Ffp_MO,kh) ;
% 
%                     % measuring the overshoot
%                     % ---------------------------------------------------------     
%                     overshoot_Fcl_MO_struct = stepinfo(Fcl_MO)          ;
%                     overshoot_Fcl_MO = overshoot_Fcl_MO_struct.Overshoot;
%                     ovrst_Lag = overshoot_Fcl_MO                        ;
%                     k = mod(stepNo_tx_MO,plotNo);
%                           if (k == 0)
%                                fprintf('stepNo_tx_MO: %d - tx_MO: %1.5f - ovrst_Lag: %1.5f - ti_MO: %1.5f - errorRef: %1.5f - if_index: %1.5f\n',stepNo_tx_MO,tx_MO,ovrst_Lag,ti_MO,errorRef,if_index)
%                           end
%                     % Calculate the new errorRef
%                     % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                     errorRef = referenceOvershoot - ovrst_Lag  ;      
%                 end
% 
%             end %of IF - line:269
% 
%         end %of WHILE  - line:267
%         
%         
%     end %of IF - line:165

