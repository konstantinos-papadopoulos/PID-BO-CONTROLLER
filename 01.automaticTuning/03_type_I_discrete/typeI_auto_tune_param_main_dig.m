% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       typeI_auto_tune_param_main_dig.m																			   																		  	
%  Project:     Automatic tuning of the parameters for PI,PID controllers (Discrete Systems)
%  
%  Purpose:     main script for the automatic tuning based on the three methods																		   																		
%  Author :     kostas g. papadopoulos																	   																		
% 																										   																		
%  History:     Date: 18.02.2009  date last modified
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

figureIndex = 0 ;                  % initialize FigureIndex
kp = 2*rand     ;
kh = 1          ;                  % Essential Condition for TYPE-I systems


% *************************************************************************
disp('**************** Include zeros in your analysis  *******************')
includeZer = input('y/n:','s');
if  strcmp(includeZer,'y') == 1
    inclZer = 1;
    fprintf('....................................\n')
    fprintf('enter the number of zeros [max:4]...\n')
    NoOfZeros = input('number of zeros:');    
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

fprintf('....................................................................\n')
fprintf('..................Enter the number of poles [max:5].................\n')
NoOfPoles = input('number of poles:');

% Matrix Initialization
% *************************************************************************
NoOfPolesMatr = zeros(5,1);
NoOfZerosMatr = zeros(4,1);


% RANDOM plant generation
% *************************************************************************
disp('********************************************************************')
disp('Generate random plant or not...(y: for random)')
disp('********************************************************************')

randomPlant = input('y/n:','s');
if  strcmp(randomPlant,'y') == 1
    % Random Zeros Generation
    % ************************
    if inclZer == 1
        for i = 1:NoOfZeros
          NoOfZerosMatr(i) = rand;
        end
    end

    for j = 1:NoOfPoles
      NoOfPolesMatr(j) = rand;
    end
    Tp6 = 0.1*rand           ;
    
    Tz1 = NoOfZerosMatr(1);     Tp1 = NoOfPolesMatr(1);  
    Tz2 = NoOfZerosMatr(2);     Tp2 = NoOfPolesMatr(2);
    Tz3 = NoOfZerosMatr(3);     Tp3 = NoOfPolesMatr(3);
    Tz4 = NoOfZerosMatr(4);     Tp4 = NoOfPolesMatr(4);
                                Tp5 = NoOfPolesMatr(5);
    Td  = inclTDel*rand   ;  
    scaleTsc = 0.01        ;
    Tp6 = scaleTsc*Tp1    ;

else
fprintf('...................Run the Debug Version with a, b..................\n')
fprintf('....................................................................\n')
    debugVer = input('y/n:','s');
    if strcmp(debugVer,'y') == 0
        Tz = 'Tz';
        Tp = 'Tp';
       if strcmp(includeZer,'y') == 1
            for i = 1:NoOfZeros
                istr = num2str(i);
                space = ' ';
                inputArg = strcat('Enter',space,Tz,istr,':');
                NoOfZerosMatr(i) = input(inputArg);
            end
       end
       
        for j = 1:NoOfPoles
            space = ' ';
            jstr = num2str(j);
            inputArg = strcat('Enter',space,Tp,jstr,':');
            NoOfPolesMatr(j) = input(inputArg);
        end
        Tz1 = NoOfZerosMatr(1);     Tp1 = NoOfPolesMatr(1);  
        Tz2 = NoOfZerosMatr(2);     Tp2 = NoOfPolesMatr(2);
        Tz3 = NoOfZerosMatr(3);     Tp3 = NoOfPolesMatr(3);
        Tz4 = NoOfZerosMatr(4);     Tp4 = NoOfPolesMatr(4);
                                    Tp5 = NoOfPolesMatr(5);
        scaleTsc = 0.01      ;
        Tp6 = scaleTsc*Tp1  ;
        if inclTDel == 1
            disp('Enter Td - (0.01 <= Td <= 0.9)')
            Td = input('Td = ');
        else
            Td = 0;
        end
    else
        disp('Enter a - (0.01 <= a <= 0.9)')
        a = input('a = ');
        a_power = 0;
        for j = 1:NoOfPoles
           NoOfPolesMatr(j) = a^(a_power);
           a_power = a_power + 1;
        end
        
        if inclZer == 1
            disp('Enter b - (0.01 <= b <= 0.9)')
            b = input('b = ');
            b_power = 0;
               for i = 1:NoOfZeros
                   NoOfZerosMatr(i) = b^(b_power);
                   b_power = b_power + 1         ;
               end
        end
        if inclTDel == 1
            disp('Enter Td - (0.01 <= Td <= 0.9)')
            Td = input('Td = ');
        else
            Td = 0;
        end
       
        
        Tz1 = NoOfZerosMatr(1);     Tp1 = NoOfPolesMatr(1);  
        Tz2 = NoOfZerosMatr(2);     Tp2 = NoOfPolesMatr(2);
        Tz3 = NoOfZerosMatr(3);     Tp3 = NoOfPolesMatr(3);
        Tz4 = NoOfZerosMatr(4);     Tp4 = NoOfPolesMatr(4);
                                    Tp5 = NoOfPolesMatr(5);
        scaleTsc = 0.01       ;
        Tp6 = scaleTsc*Tp1    ;
    end
end

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                

sortMatrtemp = [Tp1 Tp2 Tp3 Tp4 Tp5]   ;
sortMatr = sort(sortMatrtemp,'descend');
Tp1 = sortMatr(1); Tp2 = sortMatr(2)   ;
Tp3 = sortMatr(3); Tp4 = sortMatr(4)   ;
                   Tp5 = sortMatr(5)   ;

                   Ts = 0.01*Tp1        ;
% Normalizing time constants with Ts
% *************************************************************************
tz1 = Tz1 / Ts;    tp1 = Tp1 / Ts;
tz2 = Tz2 / Ts;    tp2 = Tp2 / Ts;
tz3 = Tz3 / Ts;    tp3 = Tp3 / Ts;
tz4 = Tz4 / Ts;    tp4 = Tp4 / Ts;
                   tp5 = Tp5 / Ts;
                   tp6 = Tp6 / Ts;
td = Td / Ts  ;

% Creating the plant Structure
% *************************************************************************
plant.tp1 = tp1;    plant.tz1 = tz1;
plant.tp2 = tp2;    plant.tz2 = tz2;
plant.tp3 = tp3;    plant.tz3 = tz3;
plant.tp4 = tp4;    plant.tz4 = tz4;
plant.tp5 = tp5;
plant.tp6 = tp6;
plant.td = td  ;    plant.kp = kp;   plant.kh = kh;

% Calling the Optimal I Controller [Digital]
% *************************************************************************
[Gc_i_dig Gp stepInformation_i_dig ti_MO_i_dig] = auto_tune_param_main_i_dig(plant) ;
[Gc_i     Gp stepInformation_i     ti_MO_i] = auto_tune_param_main_i(plant)         ;

Ffp_MO_dig_i = Gc_i_dig*Gp          ;
Fcl_MO_i_da = feedback(Ffp_MO_dig_i,kh);
Fcl_MO_dig_i = c2d(Fcl_MO_i_da,1,'zoh');

Ffp_MO_i = Gc_i*Gp              ;
Fcl_MO_i = feedback(Ffp_MO_i,kh);


% Calling the Optimal PI Controller [Digital]
% *************************************************************************
[Gc_pi_dig Gp stepInformation_pi_dig ti_MO_pi_dig x_MO_pi_dig] = auto_tune_param_main_pi_dig(plant)   ;
[Gc_pi     Gp stepInformation_pi ti_MO_pi x_MO_pi] = auto_tune_param_main_pi(plant)   ;

Ffp_MO_dig_pi = Gc_pi_dig*Gp              ;
Fcl_MO_pi_da  = feedback(Ffp_MO_dig_pi,kh);
Fcl_MO_dig_pi = c2d(Fcl_MO_pi_da,1,'zoh') ;

Ffp_MO_pi = Gc_pi*Gp              ;
Fcl_MO_pi = feedback(Ffp_MO_pi,kh);

% Calling the Optimal PID Controller [Digital]
% *************************************************************************
[Gc_pid_dig Gp stepInformation_pid_dig ti_MO_pid_dig x_MO_pid_dig y_MO_pid_dig] = auto_tune_param_main_pid_dig(plant)   ;
[Gc_pid     Gp stepInformation_pid     ti_MO_pid x_MO_pid y_MO_pid] = auto_tune_param_main_pid(plant)   ;


Ffp_MO_dig_pid = Gc_pid_dig*Gp              ;
Fcl_MO_pid_da = feedback(Ffp_MO_dig_pid,kh) ;
Fcl_MO_dig_pid = c2d(Fcl_MO_pid_da,1,'zoh') ;


Ffp_MO_pid = Gc_pid*Gp              ;
Fcl_MO_pid = feedback(Ffp_MO_pid,kh);
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% figure(1)
% step(1-Fcl_MO_dig_i,'b',1-Fcl_MO_i,'--b',1-Fcl_MO_dig_pi,'r',1-Fcl_MO_pi,'--r',1-Fcl_MO_dig_pid,'k',1-Fcl_MO_pid,'--k')
% figure(2)
% step(Fcl_MO_dig_i,'b',Fcl_MO_i,'--b',Fcl_MO_dig_pi,'r',Fcl_MO_pi,'--r',Fcl_MO_dig_pid,'k',Fcl_MO_pid,'--k')
% grid on
figure(3)
step(Fcl_MO_i_da,'--b',Fcl_MO_i,'b',Fcl_MO_pi_da,'--r',Fcl_MO_pi,'r',Fcl_MO_pid_da,'--k',Fcl_MO_pid,'k')
grid on
ltiview(Fcl_MO_i_da,'--b',Fcl_MO_i,'b',Fcl_MO_pi_da,'--r',Fcl_MO_pi,'r',Fcl_MO_pid_da,'--k',Fcl_MO_pid,'k')

