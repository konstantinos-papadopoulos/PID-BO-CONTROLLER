% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       automatic_tuning_main.m																			   																		  	
%  Project:     Automatic tuning of the parameters for PI,PID controllers
%  
%  Purpose:     main script for the automatic tuning based on the three methods																		   																		
%  Author :     kostas g. papadopoulos																	   																		
% 																										   																		
%  History:     Date: 07.07.2008  date last modified
% 																										  																		
%  Contact:     kostas g. papadopoulos,    nikos mitrakis       leonidas droukas        mermikli konstantina
%               kpapadop@eng.auth.gr  ,    nmitr@auth.gr        leon_drouk@yahoo.gr     kmermikl@auth.gr
% 																										  																		
%  Place:	    Aristotle University of Thessaloniki, Thessaloniki, Greece							   																		
% 
% ************************************************************************************************************
clc
clear all
close all

disp('********************************************************************')
disp('********************************************************************')
workingDir = cd;

disp('******************** Examine Type-I systems  ***********************')
type_I = input('y/n:','s');

if  strcmp(type_I,'y') == 1
    % call the type_I main Function
    % *********************************************************************
    cd(workingDir)
    cd 01_type_I
    typeI_auto_tune_param_main()
end

disp('******************** Examine Type-II systems  **********************')
type_II = input('y/n:','s');
if  strcmp(type_II,'y') == 1
    % call the type_I main Function
    % *********************************************************************
    cd(workingDir)
    cd 02_type_II
    typeII_auto_tune_param_main()
end

disp('************** Examine Type-I-Discrete - systems  ******************')
type_I_discrete = input('y/n:','s');
if  strcmp(type_I_discrete,'y') == 1
    % call the type_I main Function
    % *********************************************************************
    cd(workingDir)
    cd 03_type_I_discrete
    typeI_auto_tune_param_main_dig()
end