% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       typeII_auto_tune_param_main.m																		   																		  	
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
clc
clear all
close all
% *************************************************************************
% ********************    Random plant Generation  ************************
% *************************************************************************
kp = 2*rand     ;
kh = 1          ;                  % Essential Condition for TYPE-II systems
% *************************************************************************
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
    scaleTsc = 0.1        ;
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
        scaleTsc = 0.1      ;
        Tp6 = scaleTsc*Tp1  ;
        if inclTDel == 1
            disp('Enter Td - (0.01 <= Td <= 0.9)')
            Td = input('Td = ');
        else
            Td = 0;
        end
    else
        disp('Enter a - (0.01 <= a <= 0.9)')
        kp = 1;
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
        scaleTsc = 0.1        ;
        Tp6 = scaleTsc*Tp1    ;
    end
end

disp('********************************************************************')
plotNo = input('plottimes = ');
disp('********************************************************************')


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

% Starting the Automatic Tuning
% *************************************************************************
% Initial Conditions
% *************************************************************************
figureIndex = 1         ;
ovs_ref = 45            ;
tn_init = 10            ;
plant.tn = tn_init      ;
upperlimitOvrst = 1.02*ovs_ref;
lowerlimitOvrst = 0.98*ovs_ref;

[Fcl_MO Fcl_optimal Gp_loc stepinformation ti_MO x_MO y_MO optimalControllerDataLocal] = autotune_param_typeII_calculateOvs(plant);
Gp = Gp_loc;

while isnan(stepinformation.Overshoot) == 1
     % Increasing the Initial Condition if tn_init is not apporopriate
     % and the closed loop response is unstable
     %*********************************************************************
     tn_init = tn_init + 0.25;
     plant.tn = tn_init      ;
     [Fcl_MO Fcl_optimal Gp_loc stepinformation ti_MO x_MO y_MO optimalControllerDataLocal] = autotune_param_typeII_calculateOvs(plant);
end
    plant.tn = tn_init      ;
    ovs = stepinformation.Overshoot;
    error = ovs_ref - ovs          ;
    stepNo_tn = 1                  ;
    switchParam = 0.015            ;
    relaxCounter = 1               ;

while (ovs < lowerlimitOvrst) || (ovs > upperlimitOvrst)
    stepNo_tn = stepNo_tn + 1;
    if error > 0
        tn_auto = abs(x_MO - ( x_MO*switchParam ));
        plant.tn = tn_auto                ;
        [Fcl_MO Fcl_optimal Gp_loc stepinformation ti_MO x_MO y_MO optimalControllerDataLocal] = autotune_param_typeII_calculateOvs(plant);
        ovs = stepinformation.Overshoot;
        error = ovs_ref - ovs          ;
        k = mod(stepNo_tn,plotNo)     ;
        if (k == 0)
           figure(figureIndex)
           step(Fcl_MO)
           title('Automatic Tuning [based on overshoot] - PI Controller')
           hold on
        end
    else
        tn_auto = abs(x_MO + ( x_MO*switchParam ));
        plant.tn = tn_auto                ;
        [Fcl_MO Fcl_optimal Gp_loc stepinformation ti_MO x_MO y_MO optimalControllerDataLocal] = autotune_param_typeII_calculateOvs(plant);
        ovs = stepinformation.Overshoot;
        error = ovs_ref - ovs          ;
        k = mod(stepNo_tn,plotNo)      ;
        if (k == 0)
           figure(figureIndex)
           step(Fcl_MO)
           title('Automatic Tuning [based on overshoot] - PI Controller')
           hold on
        end
    end
        n = mod(stepNo_tn,1);
        if (n == 0)
            fprintf('step_tn: %d - ti_auto: %1.5f - x_MO: %1.5f - y_MO: %1.5f - ovs: %1.5f\n',stepNo_tn,ti_MO,x_MO,y_MO,ovs)
        end
            % relaxing the band
            % ----------------------------------------------------------------------------------------------
            if mod(stepNo_tn,relaxCounter*100) == 0
               relaxCounter = relaxCounter + 1;
                % widen the band
                % ----------------------------------------------------------------------------------------------
                upperlimitOvrst = upperlimitOvrst + 0.01*upperlimitOvrst; 
                lowerlimitOvrst = lowerlimitOvrst - 0.01*lowerlimitOvrst;
                disp('----------------------------------------------------------------------');
                fprintf('Relaxing the Reference Band so that the PI controller tuning....\n')
                fprintf('converges faster...\n')
                fprintf('New Upper Limit: %2.5f - New Lower Limit: %2.5f \n',upperlimitOvrst,lowerlimitOvrst)
                disp('----------------------------------------------------------------------');
            end

end
Gp = Gp_loc;

% Logging the Optimal Controller
% *********************************************
Fcl_auto_PI = Fcl_MO                      ;
So_MO_auto_PI  = 1 - Fcl_auto_PI          ;
Si_MO_auto_PI  = series(So_MO_auto_PI,Gp) ;

% Measured Input for the Fuzzy Classifier, Estimator
% Logging the values that are going to be the inputs for the Fuzzy estimator
% ------------------------------------------------------------------------------------------------------------
riseTime_PI_auto_comp = stepinformation.RiseTime        ;
settlingTime_PI_auto_comp = stepinformation.SettlingTime;

pi_estimator_input = zeros(4,1); 

pi_estimator_input(1) = ti_MO                    ;   % ti_PI
pi_estimator_input(2) = x_MO                     ;   % x_PI
pi_estimator_input(3) = riseTime_PI_auto_comp    ;   % riseTime_PI
pi_estimator_input(4) = settlingTime_PI_auto_comp;   % setlingTime_PI



% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   PI Fuzzy Classifier   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% ######################        FUZZY TUNING   ############################
% ######################         PI Control    ############################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
disp('----------------------------------------------------------------------');
tr_choice = input('Train Fuzzy PI Estimator [ovs. estimation]?: y/n:','s');

if  strcmp(tr_choice,'y') == 1
        
    % Estimator Training
    % *********************************************************************
    % Train the estimator so that the PI controller tuning is made through
    % an optimal overshoot reference estimated by the fuzzy estimator
    % *********************************************************************
    
    % Estimator training
    % *********************************************************************
    load ov_pi_train.txt
    
    ov_pi_control = ov_pi_train;

   
    %Formulate the training data set
    index = 1;
    k = 1    ;
    Xtrain_pi_control = ov_pi_control;
    
    % Use the Substractive Clustering method to generate a new Sugeno
    % type fuzzy system
    % *********************************************************************
    
    disp('Training the PI control estimator ....')
 
    % We use the default values of Jiang for the method, fismat contains
    % the initial fuzzy system
    % *********************************************************************
    fismat = genfis2(Xtrain_pi_control(:,1:4),Xtrain_pi_control(:,5),0.2); % We use the default values 

    % Train the system using the anfis editor
    % *********************************************************************
    [fismat_control_pi,error] = anfis(Xtrain_pi_control,fismat,[100],[0 0 0 0]);

    %Calculate the ouput of ANFIS for the training data set
    % *********************************************************************
    % trn_out_fismat_control_pi = evalfis(Xtrain_pi_control(:,1:4),fismat_control_pi);
    trn_out_fismat_control_pi = evalfis(Xtrain_pi_control(:,1:4),fismat_control_pi);

else

    % Load the already trained fuzzy PI overshoot estimator
    % *********************************************************************
    load fismat_est_pi
%     load ov_pi_zeros.txt
end
    % Input for the Fuzzy Estimator (optimal overshoot)
    % *********************************************************************
    % 1. ti_pi [is coming from the tuning based on 43% ovs]
    % 2. x_pi [is coming from the tuning based on 43% ovs]
    % 3. rise_time [is coming from the tuning based on 43% ovs]
    % 4. settling_time [is coming from the tuning based on 43% ovs]

    ovPIdesired = evalfis([pi_estimator_input(1) pi_estimator_input(2)...
                           pi_estimator_input(3) pi_estimator_input(4)],...
                           fismat_control_pi);
struct_stpInfo_Fcl_optimal = stepinfo(Fcl_optimal);

% Logging Optimal Characteristics
% *********************************************************************
ovPIOptimal = struct_stpInfo_Fcl_optimal.Overshoot;
ti_PI_optimal = optimalControllerDataLocal(1)     ;
tn_PI_optimal = optimalControllerDataLocal(2)     ;
rise_time_PI_optimal = struct_stpInfo_Fcl_optimal.RiseTime ;
set_time_PI_optimal  = struct_stpInfo_Fcl_optimal.SettlingTime ;
% *********************************************************************

fprintf('auto_tune_characteristics: [1] ti_auto: %5.5f, [2] tn_auto:  %5.5f, [3] setTime_auto:  %5.5f, [2] riseTime_auto:  %5.5f,  \n',...
        ti_PI_optimal,tn_PI_optimal,rise_time_PI_optimal,set_time_PI_optimal)
fprintf('opti_tune_characteristics: [1] ti_opti: %5.5f, [2] tn_opti:  %5.5f, [3] setTime_opti:  %5.5f, [2] riseTime_opti:  %5.5f,  \n',...
         pi_estimator_input(1),pi_estimator_input(2),pi_estimator_input(3),pi_estimator_input(4))

fprintf('ovPIdesired: %d - ovPIoptimal: %5.5f \n',ovPIdesired,ovPIOptimal)


% Starting the Automatic Tuning based on the Optimal Overshoot that came
% fromt the Estimator
% *************************************************************************
% Initial Conditions
% *************************************************************************
figureIndex = figureIndex + 1 ;
ovs_ref = ovPIdesired         ;
tn_init = 7.7                 ;
plant.tn = tn_init            ;
upperlimitOvrst = ovPIdesired +  0.001*ovPIdesired;
lowerlimitOvrst = ovPIdesired -  0.001*ovPIdesired;

[Fcl_MO Fcl_optimal Gp_loc stepinformation ti_MO x_MO y_MO] = autotune_param_typeII_calculateOvs(plant);
while isnan(stepinformation.Overshoot) == 1
     % Increasing the Initial Condition if tn_init is not apporopriate
     % and the closed loop response is unstable
     %*********************************************************************
     tn_init = tn_init + 0.25;
     plant.tn = tn_init      ;
     [Fcl_MO Fcl_optimal Gp_loc stepinformation ti_MO x_MO y_MO] = autotune_param_typeII_calculateOvs(plant);
end
    plant.tn = tn_init      ;
    ovs = stepinformation.Overshoot;
    error = ovs_ref - ovs          ;
    stepNo_tn = 1                  ;
    switchParam = 0.015            ;
    relaxCounter = 1               ;

while (ovs < lowerlimitOvrst) || (ovs > upperlimitOvrst)
    stepNo_tn = stepNo_tn + 1;
    if error > 0
        tn_auto = abs(x_MO - ( x_MO*switchParam ));
        plant.tn = tn_auto                ;
        [Fcl_MO Fcl_optimal Gp_loc stepinformation ti_MO x_MO y_MO] = autotune_param_typeII_calculateOvs(plant);
        ovs = stepinformation.Overshoot;
        error = ovs_ref - ovs          ;
        k = mod(stepNo_tn,plotNo)     ;
        if (k == 0)
           figure(figureIndex)
           step(Fcl_MO)
           title('Automatic Tuning [based on optimal overshoot] - PI Controller \\ Fuzzy PI Estimator')
           hold on
        end
    else
        tn_auto = abs(x_MO + ( x_MO*switchParam ));
        plant.tn = tn_auto                ;
        [Fcl_MO Fcl_optimal Gp_loc stepinformation ti_MO x_MO y_MO] = autotune_param_typeII_calculateOvs(plant);
        ovs = stepinformation.Overshoot;
        error = ovs_ref - ovs          ;
        k = mod(stepNo_tn,plotNo)      ;
        if (k == 0)
           figure(figureIndex)
           step(Fcl_MO)
           title('Automatic Tuning [based on optimal overshoot] - PI Controller \\ Fuzzy PI Estimator')
           hold on
        end
    end
        n = mod(stepNo_tn,1);
        if (n == 0)
            fprintf('step_tn: %d - ti_auto: %1.5f - x_MO: %1.5f - y_MO: %1.5f - ovs: %1.5f\n',stepNo_tn,ti_MO,x_MO,y_MO,ovs)
        end
            % relaxing the band
            % ----------------------------------------------------------------------------------------------
            if mod(stepNo_tn,relaxCounter*100) == 0
               relaxCounter = relaxCounter + 1;
                % widen the band
                % ----------------------------------------------------------------------------------------------
                upperlimitOvrst = upperlimitOvrst + 0.01*upperlimitOvrst; 
                lowerlimitOvrst = lowerlimitOvrst - 0.01*lowerlimitOvrst;
                disp('----------------------------------------------------------------------');
                fprintf('Relaxing the Reference Band so that the PI controller tuning....\n')
                fprintf('converges faster...\n')
                fprintf('New Upper Limit: %2.5f - New Lower Limit: %2.5f \n',upperlimitOvrst,lowerlimitOvrst)
                disp('----------------------------------------------------------------------');
            end

end

% Storing the closed loop transfer function controller by PI controller (automatically tuned)
% ------------------------------------------------------------------------------------------------------------

% Print out the results of the tuning
% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
disp('----------------------------------------------------------------------');
fprintf('Overshoot Fcl(s) after Tuning: %2.5f - ti: %2.5f\n',ovs,ti_MO)       
disp('----------------------------------------------------------------------');
Gp = Gp_loc                                    ;

%  Logging the Fuzzy Estimator Controller
% *********************************************
Fcl_MO_PI_est = Fcl_MO                         ;
So_MO_PI_est  = 1 - Fcl_MO_PI_est              ;
Si_MO_PI_est  = series(So_MO_PI_est,Gp)        ;

% Logging the Optimal Controller
% *********************************************
Fcl_optimal_PI = Fcl_optimal                    ;
So_MO_optimal_PI  = 1 - Fcl_optimal_PI          ;
Si_MO_optimal_PI  = series(So_MO_optimal_PI,Gp) ;


% Plot the results
% ----------------------------------------------------------------------------------------------
Fcl_MO_auto_PI_struct = stepinfo(Fcl_auto_PI)    ;
Fcl_optimal_PI_struct = stepinfo(Fcl_optimal_PI) ;
Fcl_MO_PI_est_struct = stepinfo(Fcl_MO_PI_est)   ;

settlingtime_Fcl_MO_auto_PI_struct = Fcl_MO_auto_PI_struct.SettlingTime    ;
settlingtime_Fcl_optimal_PI_struct = Fcl_optimal_PI_struct.SettlingTime    ;
settlingtime_Fcl_optimal_PI_est_struct = Fcl_MO_PI_est_struct.SettlingTime ;


FinalTime_Fcl_MO_auto_PI = 2*settlingtime_Fcl_MO_auto_PI_struct         ;
FinalTime_Fcl_optimal_PI = 2*settlingtime_Fcl_optimal_PI_struct        ;
FinalTime_Fcl_est_PI = 2*settlingtime_Fcl_optimal_PI_est_struct        ;


sortMatr_Fcl = [FinalTime_Fcl_MO_auto_PI FinalTime_Fcl_optimal_PI FinalTime_Fcl_est_PI];
FinalTime = max(sortMatr_Fcl)                                          ;



NoofSamples = 3000                           ;
t1 = 0:(3*FinalTime/NoofSamples):3*FinalTime ;
t2 = 0:(3*FinalTime/NoofSamples):2*FinalTime ;
t3 = 0:(3*FinalTime/NoofSamples):FinalTime   ;

% ----------------------------------------------------------------------------------------------
% Automatic Tuning based on ovs 43% [PI Controller]
% ***************************************************
y_Fcl_auto_PI = step(Fcl_auto_PI,t1);
y_So_auto_PI = step(So_MO_auto_PI,t2)  ;
y_Si_auto_PI = step(Si_MO_auto_PI,t3)  ;

y_auto_PI = zeros(NoofSamples + 1,1);
for i = 1:1:NoofSamples + 1
    y_auto_PI(i) = y_Fcl_auto_PI(i);
    if i > NoofSamples/3
        y_auto_PI(i) = y_auto_PI(i) + y_So_auto_PI(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_auto_PI(i) = y_auto_PI(i) + y_Si_auto_PI(i - (2/3)*NoofSamples);
        end
    end
end
figureIndex = figureIndex + 1;
figure(figureIndex)
plot(t1,y_auto_PI,'k','LineWidth',2)
hold on
% ----------------------------------------------------------------------------------------------
% Optimal Tuning [PI Controller]
% ***************************************************
y_Fcl_optimal_PI = step(Fcl_optimal_PI,t1)   ;
y_So_optimal_PI = step(So_MO_optimal_PI,t2)  ;
y_Si_optimal_PI = step(Si_MO_optimal_PI,t3)  ;

y_optimal_PI = zeros(NoofSamples + 1,1);
for i = 1:1:NoofSamples + 1
    y_optimal_PI(i) = y_Fcl_optimal_PI(i);
    if i > NoofSamples/3
        y_optimal_PI(i) = y_optimal_PI(i) + y_So_optimal_PI(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_optimal_PI(i) = y_optimal_PI(i) + y_Si_optimal_PI(i - (2/3)*NoofSamples);
        end
    end
end
plot(t1,y_optimal_PI,'-.r','LineWidth',2)
% ----------------------------------------------------------------------------------------------
% Automatic Tuning based on Fuzzy Estimator [PI Controller]
% *********************************************************
y_Fcl_est_PI = step(Fcl_MO_PI_est,t1)   ;
y_So_est_PI = step(So_MO_PI_est,t2)  ;
y_Si_est_PI = step(Si_MO_PI_est,t3)  ;

y_est_PI = zeros(NoofSamples + 1,1);
for i = 1:1:NoofSamples + 1
    y_est_PI(i) = y_Fcl_est_PI(i);
    if i > NoofSamples/3
        y_est_PI(i) = y_est_PI(i) + y_So_est_PI(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_est_PI(i) = y_est_PI(i) + y_Si_est_PI(i - (2/3)*NoofSamples);
        end
    end
end
plot(t1,y_est_PI,'-.g','LineWidth',2)

% ----------------------------------------------------------------------------------------------
x_min_fig3 = 0; 
x_max_fig3 = 3*FinalTime;
y_min_fig3 = 0;

sortMAX = [y_Si_auto_PI y_Si_optimal_PI y_Si_est_PI];
chooseMaxtemp = max(max(sortMAX));
chooseMax = max(2,1 + chooseMaxtemp);
if chooseMax == 2
    y_max_fig3 = 2.1;
else
  y_max_fig3 = 1.05*chooseMax;
end
axis([x_min_fig3 x_max_fig3 y_min_fig3 y_max_fig3])
grid
legend('F_{cl} - automatic tuning','F_{cl} - optimal tuning','F_{cl} - fuzzy estimator tuning',0)
title('Automatic Tuning under PI Control \newline Type II Closed Loop System')

cd ..
% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   PID Fuzzy Classifier   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   PID Fuzzy Estimator    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   PID Optimal Control    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

