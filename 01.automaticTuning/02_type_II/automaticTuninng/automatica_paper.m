% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       automatica_paper.m																		   																		  	
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
[Fcl_MO Fcl_optimal Gp_loc stepinformation ti_MO x_MO y_MO optimalControllerDataLocal Gc_PI_1_opt] = autotune_param_typeII_calculateOvs(plant);

% Logging the Optimal PI Controller
% *********************************************
Fcl_auto_PI = Fcl_MO                      ;
So_MO_auto_PI  = 1 - Fcl_auto_PI          ;
Si_MO_auto_PI  = series(So_MO_auto_PI,Gp) ;


[Gc_meth_2_Loc Gp_Loc stepinformation_pid_meth_2 ti_sq_meth_2_pid x_meth2_sol_pid y_meth2_sol_pid ...
   Gc_meth_1_Loc         stepinformation_pid_meth_1 ti_sq_meth_1_pid x_meth1_sol_pid y_meth1_sol_pid ] = auto_tune_param_main_pid_typeII(plant);
% Plant Error
% -------------------------------------------------------------------
e = 0.1;
Gp_error = plant_error(plant);
kh_dot = (1 + e)*kh          ;

% Construct the optimal PID transfer Function
% -------------------------------------------------------------------
Gc_meth_2_Loc; % optimal PID controller : second part no zeros 
Gc_PI_1_opt  ; % optimal PI  controller : second part no zeros
Gp_Loc       ; % plant transfer function
Gp_error     ; % plant transfer function with altered parameters

Ffp_MO_PID = Gc_meth_2_Loc*Gp_Loc    ;
Fcl_MO_PID = feedback(Ffp_MO_PID,kh) ;

Ffp_MO_PID_er = Gc_meth_2_Loc*Gp_error    ;
Fcl_MO_PID_er = feedback(Ffp_MO_PID_er,kh_dot);


% -------------------------------------------------------------------
Ffp_MO_PI_kostas = Gc_PI_1_opt*Gp_Loc    ;
Fcl_optimal_PI_kostas = feedback(Ffp_MO_PI_kostas,kh) ;
Fcl_optimal_PI = Fcl_optimal_PI_kostas;

Ffp_MO_PI_kostas_er = Gc_PI_1_opt*Gp_error    ;
Fcl_optimal_PI_kostas_er = feedback(Ffp_MO_PI_kostas_er,kh_dot) ;
Fcl_optimal_PI_er = Fcl_optimal_PI_kostas_er;

ltiview(Fcl_optimal_PI,'k',Fcl_optimal_PI_er,'k-.')
ltiview(Fcl_MO_PID,'k',Fcl_MO_PID_er,'k-.')


