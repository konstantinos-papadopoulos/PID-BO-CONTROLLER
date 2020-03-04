% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_main.m																			   																		  	
%  Project:     Automatic tuning of the paramters for PI,PID controllers
%  
%  Purpose:     main script for the automatic tuning based on the three methods																		   																		
%  Author :     kostas g. papadopoulos																	   																		
% 																										   																		
%  History:     Date: 19.02.2008  date last modified
% 																										  																		
%  Contact:     kostas g. papadopoulos,    leonidas droukas
%               kpapadop@eng.auth.gr  ,    leon_drouk@yahoo.gr
% 																										  																		
%  Place:	    Aristotle University of Thessaloniki, Thessaloniki, Greece							   																		
% 
% ************************************************************************************************************

clc;
clear all;

% *************************************************************************
% ********************    Random plant Generation  ************************
% *************************************************************************
kp = 1     ;
kh = 1        ; % Essential Condition for TYPE-I systems
% *************************************************************************
disp('**************** Include zeros in your analysis  *******************')
includeZer = input('y/n:','s');
if  strcmp(includeZer,'y') == 1
    inclZer = 1;
else
    inclZer = 0;
end
disp('************  Include time delay in your analysis  *****************')
includetimeDel = input('y/n:','s');
if  strcmp(includetimeDel,'y') == 1
    inclTDel = 1;
else
    inclTDel = 0;
end


% RANDOM plant generation
% *************************************************************************
Tz1 = inclZer*rand  ;     Tp1 = 0.8236; 
Tz2 = inclZer*rand  ;     Tp2 = 0.7226;
Tz3 = inclZer*rand  ;     Tp3 = 0.7213;
Tz4 = inclZer*rand  ;     Tp4 = 0.3218;
                          Tp5 = 0.1128;
Td  = inclTDel*rand ;     Tp6 = 0.0129;

Tsum = Tp1 + Tp2 + Tp3 + Tp4 + Tp5 + Tp6;

% Normalizing time constants with Tp1
% *************************************************************************
tz1 = Tz1 / Tsum ;    tp1 = Tp1 / Tsum 
tz2 = Tz2 / Tsum ;    tp2 = Tp2 / Tsum ;
tz3 = Tz3 / Tsum ;    tp3 = Tp3 / Tsum ;
tz4 = Tz4 / Tsum ;    tp4 = Tp4 / Tsum ;
                      tp5 = Tp5 / Tsum ;
                      tp6 = Tp6 / Tsum ;
td = Td / Tsum ;


% Creating the plant Structure
% *************************************************************************
plant.tp1 = tp1;    plant.tz1 = tz1;
plant.tp2 = tp2;    plant.tz2 = tz2;
plant.tp3 = tp3;    plant.tz3 = tz3;
plant.tp4 = tp4;    plant.tz4 = tz4;
plant.tp5 = tp5;
plant.tp6 = tp6;
plant.td = td  ;    plant.kp = kp;   plant.kh = kh;

pGp0 = 1 ;
pGp1 = tp1 + tp2 + tp3 + tp4 + tp5                                               ;

pGp2 = tp1*tp2 +tp1*tp3 +tp1*tp4 +tp1*tp5 +tp2*tp3 + ...
       tp2*tp4 +tp2*tp5  +tp3*tp4 +tp3*tp5 +tp4*tp5                              ;

pGp3 = tp1*tp2*tp3 + tp1*tp2*tp4 + tp1*tp2*tp5 + tp1*tp3*tp4 + tp1*tp3*tp5 +...
       tp1*tp4*tp5 + tp2*tp3*tp4 + tp2*tp3*tp5 + tp2*tp4*tp5 + tp3*tp4*tp5       ;

pGp4 = tp1*tp2*tp3*tp4 + tp1*tp2*tp3*tp5 + tp1*tp2*tp4*tp5 + tp1*tp3*tp4*tp5 +...
       tp2*tp3*tp4*tp5                                                           ;
   
pGp5 = tp1*tp2*tp3*tp4*tp5                                                       ;


% Control of the approximate plant [I Controller]
% *************************************************************************
% we conceive the plant as a first order plant
numGp = kp;
denGp = [pGp1 1];
Gpappr = tf(numGp,denGp);

% Control of the real plant [I Controller]
% *************************************************************************
numGp = kp;
denGp_tp1_tp2 = conv([tp1 1],[tp2 1]);
denGp_tp3_tp4 = conv([tp3 1],[tp4 1]);
denGp_tp1_tp2_tp3_tp4 = conv(denGp_tp1_tp2,denGp_tp3_tp4);
denFinal = conv(denGp_tp1_tp2_tp3_tp4,[tp5 1]);
Gpreal = tf(numGp,denFinal);

x_MO  = 0               ;
y_MO  = 0               ;
w1 = (Tsum + Tp6)/Tsum  ;
ti_MO = 2*kp*w1         ;

num0Gc = 1                       ;
num1Gc = x_MO                    ;
num2Gc = y_MO                    ;
numGc_MO = [num2Gc num1Gc num0Gc];
    
den0Gc = 0                       ;
den1Gc = ti_MO                   ;
den2Gc = ti_MO*tp6               ;
denGc_MO = [den2Gc den1Gc den0Gc];
Gc_I = tf(numGc_MO,denGc_MO)     ;

Ffp_MO_appr= Gc_I*Gpappr        ;
Ffp_MO_real = Gc_I*Gpreal       ;

Fcl_MO_appr = feedback(Ffp_MO_appr,kh) ;
Fcl_MO_real = feedback(Ffp_MO_real,kh) ;
So_appr = 1 - Fcl_MO_appr;
So_real = 1 - Fcl_MO_real;

t = 0:0.01:15;
[ycl_MO_appr t] = step(Fcl_MO_appr,t);
[ycl_MO_real t] = step(Fcl_MO_real,t);

[so_MO_appr t] = step(So_appr,t);
[so_MO_real t] = step(So_real,t);

figure(1)
plot(t,ycl_MO_appr,'r--','LineWidth',2)
hold on

plot(t,ycl_MO_real,'k','LineWidth',2)
hold on

plot(t,so_MO_appr,'r--','LineWidth',2)
hold on

plot(t,so_MO_real,'k','LineWidth',2)
grid

% BodeMagnitude
w = logspace(-3,3,15000);
for i = 1:15000;
     
      s = j*w(i);
      % ******************************************************************
      % real plant
      % ******************************************************************
      numReal(i) = kp                     ;
      denReal(i) = pGp5*s^5 + pGp4*s^4 + ...
                   pGp3*s^3 + pGp2*s^2 + ...
                   pGp1*s^1 + pGp0*s^0    ;
      FclReal(i) = numReal(i)/denReal(i)  ;
      SoReal(i)  = 1 - FclReal(i)         ; 
      magReal_Fcl(i) = abs(FclReal(i))    ;
      magReal_So(i)  = abs(SoReal(i))     ;
      % ******************************************************************
      % approximate plant      
      % ******************************************************************
      numAppr(i) = kp                   ;
      denAppr(i) = pGp2*s^2 + pGp1*s^1 + pGp0*s^0  ;
      FclAppr(i) = numAppr(i)/denAppr(i);
      SoAppr(i) = 1 - FclAppr(i)        ;
      magAppr_Fcl(i) = abs(FclAppr(i))  ; 
      magAppr_So(i) = abs(SoAppr(i))    ;
end
figure(2)
loglog(w,magReal_Fcl,'k','LineWidth',3)
hold on
loglog(w,magAppr_Fcl,'r--','LineWidth',3)
hold on
loglog(w,magReal_So,'K','LineWidth',3)
hold on
loglog(w,magAppr_So,'r--','LineWidth',3)
grid on

xmin = 10^(-3);
xmax = 10^3   ;
ymin = 10^(-2);
ymax = 2      ;
axis([xmin xmax ymin ymax]);
return
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                          2nd Order System
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

