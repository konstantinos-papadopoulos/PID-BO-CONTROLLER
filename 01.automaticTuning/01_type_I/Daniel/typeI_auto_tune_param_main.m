% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       typeI_auto_tune_param_main.m																			   																		  	
%  Project:     Automatic tuning of the parameters for PI,PID controllers
%  
%  Purpose:     main script for the automatic tuning based on the three methods																		   																		
%  Author :     kostas g. papadopoulos																	   																		
% 																										   																		
%  History:     Date: 10.04.2012  date last modified
% 																										  																		
%  Contact:     kostas g. papadopoulos
%               kpapadop@eng.auth.gr     
% 																										  																		
%  Place:	    Aristotle University of Thessaloniki, Thessaloniki, Greece							   																		
% 
% ************************************************************************************************************
% *************************************************************************
% ********************    Random plant Generation  ************************
% *************************************************************************
clc
clear all
figureIndex = 0 ; % initialize FigureIndex
kp = 1          ;
kh = 1          ; % Essential Condition for TYPE-I systems
Tz1 = 0         ;
Tz2 = 0         ;
Tz3 = 0         ;
Tz4 = 0         ;

Tm = 0.0033     ; % [ms]
R_sigma = 12e-3 ; % [Ohm]
L_sigma = 477e-6; % [H]
T_k = L_sigma / R_sigma; %[s]
T_net = 0.455;    % [ms]
Tp1 = T_k    ;    %[s]
Tp2 = Tm     ;    %[s]
Tp3 = 0.05*T_k;   %[s]
Tp4 = 0;
Tp5 = 0;
Tp6 = 0.025*Tp3;
Td  = 0;
% Normalizing time constants with Tp1
% *************************************************************************
tz1 = Tz1 / Tp1;    tp1 = Tp1 / Tp1;
tz2 = Tz2 / Tp1;    tp2 = Tp2 / Tp1;
tz3 = Tz3 / Tp1;    tp3 = Tp3 / Tp1;
tz4 = Tz4 / Tp1;    tp4 = Tp4 / Tp1;
                    tp5 = Tp5 / Tp1;
                    tp6 = Tp6 / Tp1;
td = Td / Tp1  ;


% Creating the plant Structure
% *************************************************************************
plant.tp1 = tp1;    plant.tz1 = tz1;
plant.tp2 = tp2;    plant.tz2 = tz2;
plant.tp3 = tp3;    plant.tz3 = tz3;
plant.tp4 = tp4;    plant.tz4 = tz4;
plant.tp5 = tp5;
plant.tp6 = tp6;
plant.td = td  ;    plant.kp = kp;   plant.kh = kh;

[tsp_est kp_est plant_tss GpLoc]= auto_tune_param_mainTspEstimation(plant);

kp = kp_est          ;   % R_sigma
plant.kp = kp        ;
Gp = GpLoc;

% Optimal controller parameters calculated after automatic tuning
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% ++++++++++++++++++   Closed Loop Control Information ++++++++++++++
% ++++++++++++++++++          PID Controller           ++++++++++++++
% ti_auto : 0.02416 - tnx_auto: 0.99225 - tvx_auto: 0.11963
% ti_auto : 0.02416 - x_auto  : 1.11188 - y_auto  : 0.11870
% PID_auto_tuned_ovs: 4.54843%
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
ti_auto  = 0.02416;
tnx_auto = 0.99225;
tvx_auto = 0.11963;
num_Gc = [tnx_auto*tvx_auto tnx_auto + tvx_auto 1];
den_Gc = [ti_auto*tp6 ti_auto 0];
Gc = tf(num_Gc,den_Gc);
alpha = 0.5;