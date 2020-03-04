% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_calculateOvs_x_y.m																			   																		  	
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

%gia kathorismeno sistima kai sigekrimena tsx , tnx kai tvx ta opoia kai
%dexetai os eisodo,
%i sinartisi ipologizei to overshoot tis vimatikis apokrisis

function [ovrst Fcl ti x y] = auto_tune_param_calculateOvs_x_y(plant_loc)

y_pred = zeros(8190,3);

load fismat_ypred

tp1 = plant_loc.tp1;    tz1 = plant_loc.tz1; 
tp2 = plant_loc.tp2;    tz2 = plant_loc.tz2;
tp3 = plant_loc.tp3;    tz3 = plant_loc.tz3;
tp4 = plant_loc.tp4;    tz4 = plant_loc.tz4;
tp5 = plant_loc.tp5;   
tp6 = plant_loc.tp6;   
td = plant_loc.td  ;    kp = plant_loc.kp ;      kh = plant_loc.kh;
tsx = plant_loc.tsx;    x  = plant_loc.x  ;      
                        x2 = plant_loc.x2 ;
% Formulate plant transfer function [numerator polynomial coefficients]
% -------------------------------------------------------------------------
q1 = tz1 + tz2 + tz3 + tz4                                    ;
q2 = tz1*tz2 + tz1*tz3 + tz1*tz4 + tz2*tz3 + tz2*tz4 + tz3*tz4;
q3 = tz1*tz2*tz3 + tz1*tz2*tz4 + tz2*tz3*tz4 + tz1*tz3*tz4    ;
q4 = tz1*tz2*tz3*tz4                                          ;

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
    N = 15                     ;
    [num_d,den_d] = pade(td,N) ;
    Gd = tf(num_d,den_d)       ;
else
    Gd = 1                     ;
end
Gp = Gp1 * Gd;

% Formulate controller transfer function with the updated tsx, tnx, tvx
% -------------------------------------------------------------------------

min_y(1) = min(y_pred(:,1));
min_y(2) = min(y_pred(:,2));
max_y(1) = max(y_pred(:,1));
max_y(2) = max(y_pred(:,2));


% Check the inputs in ANFIS for extreme values 
% x from PI control
% -------------------------------------------------------------------------

if (x2 < min_y(1))
    fuzzy_y(1) = min_y(1) ;
elseif (x2 > max_y(1))
     fuzzy_y(1) = max_y(1);
else
     fuzzy_y(1) = x2      ;
end


% x from PID control
% -------------------------------------------------------------------------
if (x < min_y(2))
    fuzzy_y(2) = min_y(2);
elseif (x > max_y(2))
    fuzzy_y(2) = max_y(2);
else
    fuzzy_y(2) = x       ;
end

y = evalfis([fuzzy_y(1) fuzzy_y(2)],fismat_ypred);
% fprintf('x2: %1.5f - x3: %1.5f - y: %1.5f \n',fuzzy_y(1),fuzzy_y(2),y)


% y = x*(x - x2) ;
numGc = [y x 1];    
denGc = [(2*kp*(tsx - x)*tp6) (2*kp*(tsx - x)) 0];
Gc = tf(numGc,denGc)  ; 
Ffp = Gc * Gp         ;  % forward path transfer function
Fcl = feedback(Ffp,kh);  % closed loop trnasfer function  

S = stepinfo(Fcl)     ;
ovrst = S.Overshoot   ;
ti = 2*kp*(tsx - x)   ;
