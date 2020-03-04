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
%  History:     Date: 07.07.2008  date last modified
% 																										  																		
%  Contact:     kostas g. papadopoulos,    nikos mitrakis       leonidas droukas
%               kpapadop@eng.auth.gr  ,    nmitr@auth.gr        leon_drouk@yahoo.gr
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
% *************************************************************************
% disp('******************** Examine 2nd ORDER system  *********************')
% secondOrderSys = input('y/n:','s');
% if  strcmp(secondOrderSys,'y') == 1
%     % call the 2nd order System Function
%     % *********************************************************************
%     auto_tune_2ndOrder_main();
%     disp('continue the automatic tuning with plant of real poles,zeros, delay...')
%     realPole = input('y/n:','s');
%     if strcmp(realPole,'n') == 1
%         return
%     end
% end


% Matrix Initialization
% *************************************************************************
NoOfPolesMatr = zeros(5,1);
NoOfZerosMatr = zeros(4,1);


% RANDOM plant generation
% *************************************************************************
% disp('********************************************************************')
% disp('Generate random plant or not...(y: for random)')
% disp('********************************************************************')
disp('********************************************************************')
fprintf('Enter "1" for random plant generation.........,\n')
fprintf('Enter "2" for the debug version...............,\n')
fprintf('Enter "3" for the ARU current controller......,\n')
disp('********************************************************************')
method = input('Enter [1,2 or 3]:');
switch method
   case 1
       disp('Method is random plant generation')
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
       
       kp = input('Enter kp = ');
       sortMatrtemp = [Tp1 Tp2 Tp3 Tp4 Tp5]   ;
       sortMatr = sort(sortMatrtemp,'descend');
       Tp1 = sortMatr(1); Tp2 = sortMatr(2)   ;
       Tp3 = sortMatr(3); Tp4 = sortMatr(4)   ;
       Tp5 = sortMatr(5)   ;
   case 2
       disp('Method is debug version')
       fprintf('...................Run the Debug Version with a, b..................\n')
       fprintf('....................................................................\n')
       disp('.............. Include zeros in your analysis  .....................')
       includeZer = input('y/n:','s');
       if  strcmp(includeZer,'y') == 1
           inclZer = 1;
           fprintf('....................................\n')
           fprintf('enter the number of zeros [max:4]...\n')
           NoOfZeros = input('number of zeros:');
       else
           inclZer = 0;
       end
       
       disp('............  Include time delay in your analysis  .................')
       includetimeDel = input('y/n:','s');
       if  strcmp(includetimeDel,'y') == 1
           inclTDel = 1;
       else
           inclTDel = 0;
       end
       
       fprintf('....................................................................\n')
       fprintf('..................Enter the number of poles [max:5].................\n')
       NoOfPoles = input('number of poles:');
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
    
    kp = input('Enter kp = ');
    sortMatrtemp = [Tp1 Tp2 Tp3 Tp4 Tp5]   ;
    sortMatr = sort(sortMatrtemp,'descend');
    Tp1 = sortMatr(1); Tp2 = sortMatr(2)   ;
    Tp3 = sortMatr(3); Tp4 = sortMatr(4)   ;
    Tp5 = sortMatr(5)   ;
    otherwise
      disp('Method is ARU current controller automatic tuning')   
      kp = 1    ;
      kh = 1    ;       % Essential Condition for TYPE-I systems
      Tz1 = 0   ; Tz2 = 0; Tz3 = 0; Tz4 = 0;
      Tm = 0.0033     ; % [s]
      R_sigma = 111e-3 ; % [Ohm]
      L_sigma = 2900e-6; % [H]
      T_k = L_sigma / R_sigma; %[s]
      T_net = 0.455;    % [ms]
      Tp1 = T_k    ;    %[s]
      Tp2 = Tm     ;    %[s]
      Tp3 = 0.0125*T_k;   %[s]
      Tp4 = 0;
      Tp5 = 0;
      Tp6 = 0.0005*Tp3;
      Td  = 0;
 end
      
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

% Calling the Optimal I Controller
% *************************************************************************
[Gc_i Gp stepInformation_i ti_MO_i] = auto_tune_param_main_i(plant)    ;

Gc_i_optimal = Gc_i                                  ;
Ffp_MO_I_optimal = Gc_i_optimal*Gp                   ;
Fcl_MO_I_optimal = feedback(Ffp_MO_I_optimal,kh)     ;
Fcl_MO_I_optimalStruct = stepinfo(Fcl_MO_I_optimal);
Fcl_MO_I_optimal_ovs = Fcl_MO_I_optimalStruct.Overshoot;

upperlimitOvrst = 4.475;
referenceOvrst = 4.47  ;
lowerlimitOvrst = 4.465;

disp('********************************************************************')
plotNo = input('plottimes = ');
disp('********************************************************************')
disp('Upper Limit overshoot is 4.47%....')
% lowerlimitOvrst = input('lowerlimitOvrst = ');

disp('********************************************************************')
while (upperlimitOvrst < lowerlimitOvrst)
    disp('-------------------------------------------------------------------------')
    disp('upperLimit < lower limit...please re-enter the limits for the ovs')
    disp('-------------------------------------------------------------------------')
    lowerlimitOvrst = input('lowerlimitOvrst = ');
end

[tsp_est kp_est plant_tss GpLoc]= auto_tune_param_mainTspEstimation(plant);

kp = kp_est          ;
plant.kp = kp        ;
tsx_estim = tsp_est  ;
tsx_auto = tsx_estim /4;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% ######################    AUTOMATIC TUNING   &###########################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%tsx_auto : tuning
% -------------------------------------------------------------------------
stepNo_tsx = 0;
for param = 1:3
    clc
    % setting tnx = 0 & tvx = 0
    % -------------------------------------------------------------------------
    tnx_auto = 0;
    tvx_auto = 0;
    plant.tsx = tsx_estim;  plant.tnx = tnx_auto;  plant.tvx = tvx_auto;

    % -------------------------------------------------------------------------
    [ovrst_r1 Fcl_auto_I ti_1_auto tnx_auto tvx_auto] = auto_tune_param_calculateOvs(plant);
    
    % fprintf('make the tuning by adding the linear PI... \n')
    % fprintf('if "n", the tuning is done with the if conditions...\n')
    % linearPI = input('y/n:','s');
    linearPI = 'y';
    disp('********************************************************************')
    relaxCounter = 1;
    % investigating three settings for the PI controller
    tn_linear_matr = [5 5 5];
    ti_linear_matr = [80 100 500];
    
    if strcmp(linearPI,'y') == 1
        % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        % ++++++++++++++++++  Linear PI Controller Implementation  ++++++++++++++++
        % *************************************************************************
        switchParam = 0.015               ; % parameter for the IF condition
        figureIndex = 0                   ;
        error1 = abs(referenceOvrst - ovrst_r1); %Initial Condition for the ODE Function
        %   x_error = error1;
        tn_linear = tn_linear_matr(param)
        ti_linear = ti_linear_matr(param)
        k = 2;
        error_Ref(1) = 0     ;
        error_Ref(2) = error1;
        x_state(1) = error1  ;
        figureIndex = figureIndex +1;
        while ((ovrst_r1 < lowerlimitOvrst) || (ovrst_r1 > upperlimitOvrst)) || isnan(ovrst_r1) == 1
            stepNo_tsx = stepNo_tsx + 1;
            x_state(k) = x_state(k - 1) + (1/ti_linear)*error_Ref(k)  ;
            tsx_auto(k) = (tn_linear/ti_linear)*error_Ref(k) + x_state(k);
            plant.tsx = abs(tsx_auto(k));
            % plant.tsx = tsx_auto(k);
            [ovrst_r1 Fcl_auto_I ti_1_auto tnx_auto tvx_auto] = auto_tune_param_calculateOvs(plant);
            ovrst_r1_temp(k) = ovrst_r1 ;
            ti_1_auto_temp(k)= ti_1_auto;
            fprintf('k_step: %1.f - errorRef: %3.5f - ovs: %3.5f - tsx: %3.5f - ti: %1.5f\n',k,error_Ref(k),ovrst_r1,tsx_auto(k),ti_1_auto_temp(k));
            k = k + 1;
            % error_Ref(k) = abs(referenceOvrst - ovrst_r1);
            error_Ref(k) = referenceOvrst - ovrst_r1;
            if  isnan(error_Ref(k)) == 1
                error_Ref(k) = 0.1*error_Ref(k-1);
                if mod(stepNo_tsx,relaxCounter*10) == 0
                    relaxCounter = relaxCounter + 1;
                    % widen the band
                    % ----------------------------------------------------------------------------------------------
                    upperlimitOvrst = upperlimitOvrst + 0.2*upperlimitOvrst;
                    lowerlimitOvrst = lowerlimitOvrst - 0.2*lowerlimitOvrst;
                    disp('----------------------------------------------------------------------');
                    fprintf('Relaxing the Reference Band so that the I controller tuning....\n')
                    fprintf('converges faster...\n')
                    fprintf('New Upper Limit: %2.5f - New Lower Limit: %2.5f \n',upperlimitOvrst,lowerlimitOvrst)
                    disp('----------------------------------------------------------------------');
                end
            end
            j = mod(stepNo_tsx,plotNo);
            if (j == 0)
                figure(figureIndex)
                step(Fcl_auto_I,'k')
                title('Automatic Tuning [based on overshoot] - I Controller')
                hold on
            end
            % Logging the ovs values
            % -------------------------------------------------------------
            ovs_I_auto_tune(stepNo_tsx) = ovrst_r1;
            step_Conv_autoTuning_realZeros_I_Controller = stepNo_tsx;
        end
        %     tsx_auto = abs(tsx_auto(end));
    else
        
        % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        % +++++++++++++  TUNING is done by using the IF condition  ++++++++++++
        % *********************************************************************
        switchParam = 0.015               ; % parameter for the IF condition
        Fcltr = cell(10);
        error1 = referenceOvrst - ovrst_r1;
        figureIndex = figureIndex + 1;
        while ((ovrst_r1 < lowerlimitOvrst) || (ovrst_r1 > upperlimitOvrst))
            stepNo_tsx = stepNo_tsx + 1;
            if (error1 > 0)
                % if overshoot < 4.47 (reference value)
                % decrease tsx
                % ---------------------------------------------------------
                tsx_auto = abs(tsx_auto - (tsx_auto*switchParam));
                plant.tsx = tsx_auto;
                [ovrst_r1 Fcl_auto_I ti_1_auto tnx_auto tvx_auto] = auto_tune_param_calculateOvs(plant);
                k = mod(stepNo_tsx,plotNo);
                if (k == 0)
                    figure(figureIndex)
                    step(Fcl_auto_I,'k')
                    title('Automatic Tuning [based on overshoot] - I Controller')
                    hold on
                end
                Fcltr{stepNo_tsx} = Fcl_auto_I;
                error1 = referenceOvrst - ovrst_r1;
            elseif (error1 < 0)
                % if overshoot > 4.47 (reference value)
                % increase tsx
                % ---------------------------------------------------------
                tsx_auto = abs(tsx_auto + (tsx_auto*switchParam));
                plant.tsx = tsx_auto;
                [ovrst_r1 Fcl_auto_I ti_1_auto tnx_auto tvx_auto] = auto_tune_param_calculateOvs(plant);
                k = mod(stepNo_tsx,plotNo);
                if (k == 0)
                    figure(figureIndex)
                    step(Fcl_auto_I,'k')
                    title('Automatic Tuning [based on overshoot] - I Controller')
                    hold on
                end
                error1 = referenceOvrst - ovrst_r1;
                Fcltr{stepNo_tsx} = Fcl_auto_I;
            end
            
            % Logging the ovs values
            % -------------------------------------------------------------
            ovs_I_auto_tune(stepNo_tsx) = ovrst_r1;
            
            n = mod(stepNo_tsx,1);
            if (n == 0)
                fprintf('step_tsx: %d - ti_auto: %1.5f - tnx_auto: %1.5f - tvx_auto: %1.5f\n',stepNo_tsx,ti_1_auto,tnx_auto,tvx_auto)
            end
            % relaxing the band
            % ----------------------------------------------------------------------------------------------
            if mod(stepNo_tsx,relaxCounter*100) == 0
                relaxCounter = relaxCounter + 1;
                % widen the band
                % ----------------------------------------------------------------------------------------------
                upperlimitOvrst = upperlimitOvrst + 0.01*upperlimitOvrst;
                lowerlimitOvrst = lowerlimitOvrst - 0.01*lowerlimitOvrst;
                disp('----------------------------------------------------------------------');
                fprintf('Relaxing the Reference Band so that the I controller tuning....\n')
                fprintf('converges faster...\n')
                fprintf('New Upper Limit: %2.5f - New Lower Limit: %2.5f \n',upperlimitOvrst,lowerlimitOvrst)
                disp('----------------------------------------------------------------------');
            end
        end
        
        if stepNo_tsx == 0
            stepNo_tsx = stepNo_tsx + 1;
            ovs_I_auto_tune(stepNo_tsx) = ovrst_r1;
        end
        % Logging Convergence Steps for I Controller - Automatic Tuning Real Zeros
        % *************************************************************************
        step_Conv_autoTuning_realZeros_I_Controller = stepNo_tsx;
        
        % adding figure properties
        % *************************************************************************
        title('Automatic Tuning [based on overshoot] - I Controller')
        grid on
        
    end    
        if param == 1
            error_Ref_matr_1 = zeros(length(error_Ref),1);
            error_Ref_matr_1(:,1) = error_Ref;
            x_ax_1 = 1:1:length(error_Ref_matr_1);
            tsx_1 = ti_1_auto/(2*kp);
        elseif param == 2
            error_Ref_matr_2 = zeros(length(error_Ref),1);
            error_Ref_matr_2(:,1) = error_Ref;
            x_ax_2 = 1:1:length(error_Ref_matr_2);
            tsx_2 = ti_1_auto/(2*kp);
        else
            error_Ref_matr_3 = zeros(length(error_Ref),1);
            error_Ref_matr_3(:,1) = error_Ref;
            x_ax_3 = 1:1:length(error_Ref_matr_3);
            tsx_3 = ti_1_auto/(2*kp);
        end
end
return
close all hidden
figure(2)
plot(x_ax_1,error_Ref_matr_1,'k-d','MarkerSize',2); grid on; hold on
plot(x_ax_2,error_Ref_matr_2,'r-+','MarkerSize',2); grid on; hold on
plot(x_ax_3,error_Ref_matr_3,'b-o','MarkerSize',2); grid on;
return
% *************************************************************************
Fcl_auto_I_struct = stepinfo(Fcl_auto_I)         ;
ovsFinal_Fcl_auto_I = Fcl_auto_I_struct.Overshoot;
tsx_auto = abs(tsx_auto(end)); % log the last value of tsx_auto during its tuning
ti_auto = 2*kp*(tsx_auto - tnx_auto - tvx_auto)  ;
x_auto  = tnx_auto + tvx_auto;
y_auto  = tnx_auto * tvx_auto;

auto_tune_I_log_parameters(1,1) = ti_auto ;  auto_tune_I_log_parameters(1,2) = ti_auto;
auto_tune_I_log_parameters(2,1) = tnx_auto;  auto_tune_I_log_parameters(2,2) = x_auto ;
auto_tune_I_log_parameters(3,1) = tvx_auto;  auto_tune_I_log_parameters(3,2) = y_auto ;

% storing ti_auto from the I controller for the FUZZY PI autotuning
% *************************************************************************
ti_auto_fuzzyInput = ti_auto;

disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
disp('++++++++++++++++++   Closed Loop Control Information ++++++++++++++')
disp('++++++++++++++++++           I Controller            ++++++++++++++')
fprintf('ti_auto : %1.5f - tnx_auto: %1.5f - tvx_auto: %1.5f\n',ti_auto,tnx_auto,tvx_auto)
fprintf('ti_auto : %1.5f - x_auto  : %1.5f - y_auto  : %1.5f\n',ti_auto,x_auto,y_auto)
fprintf('I_auto_tuned_ovs: %1.5f%%\n',ovsFinal_Fcl_auto_I)
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')

if (figureIndex == 0)
   figureIndex = 1 ;
end
% tnx_auto : tuning
% -------------------------------------------------------------------------
% Initial Conditions
% -------------------------------------------------------------------------
stepNo_tnx = 0           ;
tnx_auto = tsx_auto / 2  ;  tvx_auto = 0        ;
plant.tsx = tsx_auto     ;  plant.tnx = tnx_auto;  plant.tvx = tvx_auto;

[ovrst_r2 Fcl_auto_PI ti_2_auto tnx_auto tvx_auto] = auto_tune_param_calculateOvs(plant);
error2 = referenceOvrst - ovrst_r2;
relaxCounter = 1;
while ((ovrst_r2 < lowerlimitOvrst) || (ovrst_r2 > upperlimitOvrst))
    stepNo_tnx = stepNo_tnx + 1;
    if (error2 > 0)
        % if overshoot < 4.47 (reference value)
        % increase tnx
        % ---------------------------------------------------------
        tnx_auto = abs(tnx_auto + (tnx_auto*switchParam));                
        plant.tnx = tnx_auto;
        [ovrst_r2 Fcl_auto_PI ti_2_auto tnx_auto tvx_auto] = auto_tune_param_calculateOvs(plant);
        k = mod(stepNo_tnx,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_auto_PI,'k')
           title('Automatic Tuning [based on overshoot] - PI Controller')
           hold on
        end
        error2 = referenceOvrst - ovrst_r2;
    elseif (error2 < 0)
        % if overshoot > 4.47 (reference value)
        % decrease tnx
        % ---------------------------------------------------------
        tnx_auto = abs(tnx_auto - (tnx_auto*switchParam));                
        plant.tnx = tnx_auto;
        [ovrst_r2 Fcl_auto_PI ti_2_auto tnx_auto tvx_auto] = auto_tune_param_calculateOvs(plant);
        k = mod(stepNo_tnx,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_auto_PI,'k')
           title('Automatic Tuning [based on overshoot] - PI Controller')
           hold on
        end
        error2 = referenceOvrst - ovrst_r2;
            % relaxing the band
            % ----------------------------------------------------------------------------------------------
            if mod(stepNo_tnx,relaxCounter*100) == 0
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
    % Logging the ovs values
    % -------------------------------------------------------------
    ovs_PI_auto_tune(stepNo_tnx) = ovrst_r2;

    n = mod(stepNo_tnx,1);
    if (n == 0)
      fprintf('step_tsx: %d - ti_auto: %1.5f - tnx_auto: %1.5f - tvx_auto: %1.5f\n',stepNo_tnx,ti_2_auto,tnx_auto,tvx_auto)
    end
end
    if stepNo_tnx == 0
       stepNo_tnx = stepNo_tnx + 1;
       % Logging the ovs values
       % -------------------------------------------------------------
       ovs_PI_auto_tune(stepNo_tnx) = ovrst_r2;
    end

% Logging Convergence Steps for PI Controller - Automatic Tuning Real Zeros
% *************************************************************************
step_Conv_autoTuning_realZeros_PI_Controller = stepNo_tnx;

% adding figure properties
% *************************************************************************
title('Automatic Tuning [based on overshoot] - PI Controller')      
grid on

% *************************************************************************
% *************************************************************************
Fcl_auto_PI_struct = stepinfo(Fcl_auto_PI)         ;
ovsFinal_Fcl_auto_PI = Fcl_auto_PI_struct.Overshoot;
ti_auto = 2*kp*(tsx_auto - tnx_auto - tvx_auto)    ;
x_auto  = tnx_auto + tvx_auto;
y_auto  = tnx_auto * tvx_auto;

auto_tune_PI_log_parameters = zeros(3,2);

auto_tune_PI_log_parameters(1,1) = ti_auto ;  auto_tune_PI_log_parameters(1,2) = ti_auto;
auto_tune_PI_log_parameters(2,1) = tnx_auto;  auto_tune_PI_log_parameters(2,2) = x_auto ;
auto_tune_PI_log_parameters(3,1) = tvx_auto;  auto_tune_PI_log_parameters(3,2) = y_auto ;


disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
disp('++++++++++++++++++   Closed Loop Control Information ++++++++++++++')
disp('++++++++++++++++++           PI Controller           ++++++++++++++')
fprintf('ti_auto : %1.5f - tnx_auto: %1.5f - tvx_auto: %1.5f\n',ti_auto,tnx_auto,tvx_auto)
fprintf('ti_auto : %1.5f - x_auto  : %1.5f - y_auto  : %1.5f\n',ti_auto,x_auto,y_auto)
fprintf('PI_auto_tuned_ovs: %1.5f%%\n',ovsFinal_Fcl_auto_PI)
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')

% tvx_auto : tuning
% -------------------------------------------------------------------------
stepNo_tvx = 0                      ;
tvx_auto = (tsx_auto - tnx_auto) / 2;
% tvx_auto = 2*(tsx_auto - tnx_auto) ;
plant.tsx = tsx_auto                ;  plant.tnx = tnx_auto;  plant.tvx = tvx_auto;

[ovrst_r3 Fcl_auto_PID ti_3_auto tnx_auto tvx_auto] = auto_tune_param_calculateOvs(plant);

error3 = referenceOvrst - ovrst_r3  ;
relaxCounter = 1;
auto_tune_PID_log_parameters = zeros(3,2);

while ((ovrst_r3 < lowerlimitOvrst) || (ovrst_r3 > upperlimitOvrst))
    stepNo_tvx = stepNo_tvx + 1;
    if (error3 > 0)
        % if overshoot < 4.47 (reference value)
        % increase tvx
        % ---------------------------------------------------------
        tvx_auto = abs(tvx_auto + (tvx_auto*switchParam));              
        plant.tvx = tvx_auto;
        [ovrst_r3 Fcl_auto_PID ti_3_auto tnx_auto tvx_auto] = auto_tune_param_calculateOvs(plant);
        k = mod(stepNo_tvx,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_auto_PID,'k')
           title('Automatic Tuning [based on overshoot] - PID Controller')
           hold on
        end
        error3 = referenceOvrst - ovrst_r3;
    elseif (error3 < 0)
        % if overshoot > 4.47 (reference value)
        % decrease tvx
        % ---------------------------------------------------------
        tvx_auto = abs(tvx_auto - (tvx_auto*switchParam));           
        plant.tvx = tvx_auto;
        [ovrst_r3 Fcl_auto_PID ti_3_auto tnx_auto tvx_auto] = auto_tune_param_calculateOvs(plant);
        k = mod(stepNo_tvx,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_auto_PID,'k')
           title('Automatic Tuning [based on overshoot] - PID Controller')
           hold on
        end
        error3 = referenceOvrst - ovrst_r3;
            % relaxing the band
            % ----------------------------------------------------------------------------------------------
            if mod(stepNo_tvx,relaxCounter*100) == 0
               relaxCounter = relaxCounter + 1;
                % widen the band
                % ----------------------------------------------------------------------------------------------
                upperlimitOvrst = upperlimitOvrst + 0.01*upperlimitOvrst; 
                lowerlimitOvrst = lowerlimitOvrst - 0.01*lowerlimitOvrst;
                disp('----------------------------------------------------------------------');
                fprintf('Relaxing the Reference Band so that the PID controller tuning....\n')
                fprintf('converges faster...\n')
                fprintf('New Upper Limit: %2.5f - New Lower Limit: %2.5f \n',upperlimitOvrst,lowerlimitOvrst)
                disp('----------------------------------------------------------------------');
            end
    end
    % Logging the ovs values
    % -------------------------------------------------------------
    ovs_PID_auto_tune(stepNo_tvx) = ovrst_r3;

    n = mod(stepNo_tvx,1);
    if (n == 0)
      fprintf('step_tsx: %d - ti_auto: %1.5f - tnx_auto: %1.5f - tvx_auto: %1.5f\n',stepNo_tvx,ti_3_auto,tnx_auto,tvx_auto)
    end
end
    if stepNo_tvx == 0
       stepNo_tvx = stepNo_tvx + 1;
       ovs_PID_auto_tune(stepNo_tvx) = ovrst_r3;
    end

% Logging the overshoot while tuning
% *************************************************************************
ovs_Final_auto_tune = [ovs_I_auto_tune ovs_PI_auto_tune ovs_PID_auto_tune];
    
% Logging Convergence Steps for PID Controller - Automatic Tuning Real Zeros
% *************************************************************************
step_Conv_autoTuning_realZeros_PID_Controller = stepNo_tvx;

% adding figure properties
% *************************************************************************
title('Automatic Tuning [based on overshoot] - PID Controller')      
grid on

% *************************************************************************
% *************************************************************************
Fcl_auto_PID_struct = stepinfo(Fcl_auto_PID)         ;
ovsFinal_Fcl_auto_PID = Fcl_auto_PID_struct.Overshoot;
ti_auto = 2*kp*(tsx_auto - tnx_auto - tvx_auto)      ;
x_auto  = tnx_auto + tvx_auto                        ;
y_auto  = tnx_auto * tvx_auto                        ;

auto_tune_PID_log_parameters(1,1) = ti_auto ;  auto_tune_PID_log_parameters(1,2) = ti_auto;
auto_tune_PID_log_parameters(2,1) = tnx_auto;  auto_tune_PID_log_parameters(2,2) = x_auto ;
auto_tune_PID_log_parameters(3,1) = tvx_auto;  auto_tune_PID_log_parameters(3,2) = y_auto ;

disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
disp('++++++++++++++++++   Closed Loop Control Information ++++++++++++++')
disp('++++++++++++++++++          PID Controller           ++++++++++++++')
fprintf('ti_auto : %1.5f - tnx_auto: %1.5f - tvx_auto: %1.5f\n',ti_auto,tnx_auto,tvx_auto)
fprintf('ti_auto : %1.5f - x_auto  : %1.5f - y_auto  : %1.5f\n',ti_auto,x_auto,y_auto)
fprintf('PID_auto_tuned_ovs: %1.5f%%\n',ovsFinal_Fcl_auto_PID)
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
return
% Calling the optimal Controller
% =====================================================================================
% [Gc_pi Gp stepInformation_pi ti_MO_pi x_MO_pi y_MO_pi] = auto_tune_param_main_pi(plant)      ;
% [Gc_pid Gp stepInformation_pid ti_MO_pid x_MO_pid y_MO_pid] = auto_tune_param_main_pid(plant);
% Gc_pi_optimal = Gc_pi                                ;
% 
% t = 0:0.001:100;
% figure(1)
% step(Fcl_auto_PI)
% figure(2)
% step(Fcl_auto_PID)
% return
% -------------------------------------------------------------------------
% PID Automatic Tuning based on 'y' estimation
% tvx_fuzzy : tuning
% -------------------------------------------------------------------------
% Initial Conditions
% -------------------------------------------------------------------------
stepNo_tvx = 0                       ;
% x_auto_marg = tnx_auto + (tnx_auto/2);
x_auto_marg = tnx_auto/3;
plant.tsx = tsx_auto                 ;  plant.x = x_auto_marg;  plant.x2 = tnx_auto;

figureIndex = figureIndex + 1        ;
[ovrst_r3 Fcl_PID_auto_marg ti_3_auto_marg x_auto_marg y_auto_marg] = auto_tune_param_calculateOvs_x_y_margaris(plant);

error3 = referenceOvrst - ovrst_r3   ;
relaxCounter = 1;
while ((ovrst_r3 < lowerlimitOvrst) || (ovrst_r3 > upperlimitOvrst))
    stepNo_tvx = stepNo_tvx + 1;
    if (error3 > 0)                                                  % Me vasi ta diagrammata tis selidas 13, otan overshoot < 4.47(onomastiki timi)
        x_auto_marg = abs(x_auto_marg + (x_auto_marg*switchParam));  % isxiei tvx < tv_onomastiko. Ara to tvx tha prepei na afksithei
        plant.x = x_auto_marg;
        [ovrst_r3 Fcl_PID_auto_marg ti_3_auto_marg x_auto_marg y_auto_marg] = auto_tune_param_calculateOvs_x_y_margaris(plant);
        k = mod(stepNo_tvx,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_PID_auto_marg,'k')
           title('Automatic Tuning - PID Controller - "y" Estimation (Margaris)')
           hold on
        end
        error3 = referenceOvrst - ovrst_r3   ;
    elseif (error3 < 0)                                               %Me vasi ta diagrammata tis selidas 13, otan overshoot > 4.47(onomastiki timi)
        x_auto_marg = abs(x_auto_marg - (x_auto_marg*switchParam));               %isxiei tvx > tv_onomastiko. Ara to tvx tha prepei na meiothei
        plant.x = x_auto_marg;
        [ovrst_r3 Fcl_PID_auto_marg ti_3_auto_marg x_auto_marg y_auto_marg] = auto_tune_param_calculateOvs_x_y_margaris(plant);
        k = mod(stepNo_tvx,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_PID_auto_marg,'k')
           title('Automatic Tuning - PID Controller - "y" Estimation (Margaris)')
           hold on
        end
        error3 = referenceOvrst - ovrst_r3   ;
            % relaxing the band
            % ----------------------------------------------------------------------------------------------
            if mod(stepNo_tvx,relaxCounter*100) == 0
               relaxCounter = relaxCounter + 1;
                % widen the band
                % ----------------------------------------------------------------------------------------------
                upperlimitOvrst = upperlimitOvrst + 0.01*upperlimitOvrst; 
                lowerlimitOvrst = lowerlimitOvrst - 0.01*lowerlimitOvrst;
                disp('----------------------------------------------------------------------');
                fprintf('Relaxing the Reference Band so that the PID controller tuning....\n')
                fprintf('converges faster...\n')
                fprintf('New Upper Limit: %2.5f - New Lower Limit: %2.5f \n',upperlimitOvrst,lowerlimitOvrst)
                disp('----------------------------------------------------------------------');
            end
    end
   
    % Logging the ovs values
    % -------------------------------------------------------------
    ovs_PID_auto_tune_ypred(stepNo_tvx) = ovrst_r3;

    n = mod(stepNo_tvx,1);
    if (n == 0)
      fprintf('step_x3: %d - ti_auto_marg: %1.5f - x_auto_marg: %1.5f - y_auto_marg: %1.5f \n',stepNo_tvx,ti_3_auto_marg,x_auto_marg,y_auto_marg)
    end
end
    if stepNo_tvx == 0
       stepNo_tvx = stepNo_tvx + 1;
       ovs_PID_auto_tune_ypred(stepNo_tvx) = ovrst_r3;
    end 
% Logging the overshoot while tuning
% *************************************************************************
ovs_Final_auto_tune_margaris = [ovs_I_auto_tune ovs_PI_auto_tune ovs_PID_auto_tune_ypred];

% Logging Convergence Steps for PID Controller - Automatic Tuning Complex Zeros
% *************************************************************************
step_Conv_autoTuning_complexZeros_PID_Controller = stepNo_tvx;

% adding figure properties
% *************************************************************************
title('Automatic Tuning - PID Controller - Y Estimation (Margaris)')      
grid on

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
Fcl_PID_auto_marg_struct = stepinfo(Fcl_PID_auto_marg)          ;
ovsFinal_Fcl_auto_marg_PID = Fcl_PID_auto_marg_struct.Overshoot ;

disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
disp('++++++++++++++++++   Closed Loop Control Information ++++++++++++++')
disp('++++++++++++++++++          PID Controller           ++++++++++++++')
disp('++++++++++++++++++  Automatic Tuning: "Y" Estimation ++++++++++++++')
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
fprintf('ti_auto_marg : %1.5f - x_auto_marg: %1.5f - y_auto_marg: %1.5f\n',ti_3_auto_marg,x_auto_marg,y_auto_marg)
fprintf('PID_auto_tuned_marg_ovs: %1.5f%%\n',ovsFinal_Fcl_auto_marg_PID)
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')

log_ovs_while_auto_tuning(1) = lowerlimitOvrst;
log_ovs_while_auto_tuning(2) = upperlimitOvrst;
log_ovs_while_auto_tuning(3) = referenceOvrst ;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% ######################        FUZZY TUNING   ############################
% ######################         PI Control    ############################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
tr_choice = input('Train fuzzy PI?: y/n:','s');

if  strcmp(tr_choice,'y') == 1
        
    % Classifier Training
    % *********************************************************************
    % Classifier Output: [1] if PI control exists
    % Classifier Output: [0] if PI control DOES NOT exist --> switc to
    % I-Lag controller
    % *********************************************************************
    % Train the classifier for identifying whether or not an optimal PI can
    % be developed
    % *********************************************************************
    
    % Classifier training
    % *********************************************************************
    load ov_pi_zeros.txt
    
    ov_pi_control = ov_pi_zeros;

    %Get the number of examples
    % *********************************************************************
    numExamples = size(ov_pi_zeros);
    
    %Formulate the training data set
    index = 1;
    k = 1    ;
    Xtrain_pi_control = ov_pi_control;
    
    for i = 1:numExamples
        if ((Xtrain_pi_control(i,3) < 3) || (Xtrain_pi_control(i,3) > 8))
            Xtrain_pi_control(i,3) = 0;
        end
    end
    
    for i=1:numExamples
        if Xtrain_pi_control(i,3) ~= 0
           Xtrain_pi_control(i,3) = 1 ;
        end
    end

    % Use the Substractive Clustering method to generate a new Sugeno
    % type fuzzy system
    % *********************************************************************
    
    disp('Training the PI control classifier ....')
 
    % We use the default values of Jiang for the method, fismat contains
    % the initial fuzzy system
    % *********************************************************************
    fismat = genfis2(Xtrain_pi_control(:,1:2),Xtrain_pi_control(:,3),0.2);

    %Train the system using the anfis editor
    % *********************************************************************
    [fismat_control_pi,error] = anfis(Xtrain_pi_control,fismat,[100],[0 0 0 0]);

    %Calculate the ouput of ANFIS for the training data set
    % *********************************************************************
    trn_out_fismat_control_pi = evalfis(Xtrain_pi_control(:,1:2),fismat_control_pi);


    trn_out_fismat_control_pi = round(trn_out_fismat_control_pi);
    correct = 0;
    for i = 1:numExamples(1)
        if (trn_out_fismat_control_pi(i) == Xtrain_pi_control(i,3))
            correct=correct+1;
        end
    end

    Percentage = correct/numExamples(1)*100;
    % print out the result of the training
    % *********************************************************************
    fprintf('PI Classifier Percentage = %3.5f\n',Percentage)
    
    % **********************   PREDICTION PROBLEM  ************************
    % Predict the desired PI overshoot when 
    % the plant can be controlled by a PI controller
    % *********************************************************************

    ov_pi_pred = ov_pi_zeros;
        
    % Formulate the training data set
    % *********************************************************************
    index = 1;
    k = 1    ;

    for i = 1:numExamples
        if (Xtrain_pi_control(i,3) == 1)
            Xtrain_pred_pi(index,:) = ov_pi_pred(i,:);
            index = index + 1                        ;
        end 
    end

    disp('Training fuzzy PI overshoot estimator....')
    % Grid partition with Gaussian 3 MF and linear type consequent
    % *********************************************************************
    fismat = genfis1(Xtrain_pred_pi(:,1:3),[3 3],'gaussmf','linear');

    % fismat=genfis2(Xtrain_pred_pi(:,1:2),Xtrain_pred_pi(:,3),0.2); 
    % We use the default values of Jiang for the method, fismat contains
    % the initial fuzzy system

    % Train the system using the anfis editor
    % *********************************************************************
    [fismat_pred_pi,error] = anfis(Xtrain_pred_pi,fismat,[100],[0 0 0 0]);    

else

    % Load the already trained fuzzy PI overshoot estimator
    % *********************************************************************
    load fismat_pi
    load ov_pi_zeros.txt

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                   Fuzzy PI Overshoot Estimation
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

[Gc_pi Gp stepInformation_pi ti_MO_pi x_MO_pi y_MO_pi] = auto_tune_param_main_pi(plant);
Gc_pi_optimal = Gc_pi                                ;
Ffp_MO_PI_optimal = Gc_pi_optimal*Gp                 ;

Fcl_MO_PI_optimal = feedback(Ffp_MO_PI_optimal,kh)   ;
So_MO_PI_optimal  = 1 - Fcl_MO_PI_optimal            ; 
Si_MO_PI_optimal  = series(So_MO_PI_optimal,Gp)      ;
ovsPI_optimalstruct = stepinfo(Fcl_MO_PI_optimal)    ;

ovPIoptimal = ovsPI_optimalstruct.Overshoot          ;


max_pi(1) = max(ov_pi_zeros(:,1));
min_pi(1) = min(ov_pi_zeros(:,1));
max_pi(2) = max(ov_pi_zeros(:,2));
min_pi(2) = min(ov_pi_zeros(:,2));

if (ti_auto_fuzzyInput/kp < min_pi(1))
    fuzzyPI(1) = min_pi(1);
elseif (ti_auto_fuzzyInput/kp > max_pi(1))
    fuzzyPI(1) = max_pi(1);
else
    fuzzyPI(1) = ti_auto_fuzzyInput/kp;
end



if (plant_tss < min_pi(2))
    fuzzyPI(2) = min_pi(2);
elseif (plant_tss > max_pi(2))
    fuzzyPI(2) = max_pi(2);
else
    fuzzyPI(2) = plant_tss;
end


LagFlagPI = evalfis([fuzzyPI(1) fuzzyPI(2)],fismat_control_pi);
% LagFlagPI = round(LagFlagPI)                                  ;

if (LagFlagPI >= 0.5)
    
    % return the desired ovs for the PI Controller and tune the PI controller
    % *********************************************************************
    ovPIdesired = evalfis([fuzzyPI(1) fuzzyPI(2)],fismat_pred_pi); 

    disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
    disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    disp('++++++++++++++++++   Fuzzy PI Overshoot Estimation  +++++++++++++++')
    disp('++++++++++++++++++          PI  Controller          +++++++++++++++')
    disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    fprintf('PI optimal OVS: %1.5f\n',ovPIoptimal)
    fprintf('PI desired OVS: %1.5f\n',ovPIdesired)
    disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')

    % tnx_fuzzy : tuning
    % -------------------------------------------------------------------------
    % Initial Conditions
    % -------------------------------------------------------------------------
    stepNo_tnx = 0            ;
    tnx_fuzzy = tsx_auto / 2  ;  tvx_fuzzy = 0        ;  tsx_fuzzy = tsx_auto ;
    plant.tsx = tsx_auto      ;  plant.tnx = tnx_fuzzy;  plant.tvx = tvx_fuzzy;

    [ovrst_r2 Fcl_PI_fuzzy ti_2_fuzzy tnx_fuzzy tvx_fuzzy] = auto_tune_param_calculateOvs(plant);

    error2 = ovPIdesired - ovrst_r2      ;
    lowerlimitOvrst = ovPIdesired * 0.99 ;
    upperlimitOvrst = ovPIdesired * 1.01 ;

    figureIndex = figureIndex + 1        ; %increasing the Figure Index by 1 [figureIndex = 3]
    relaxCounter = 1;
    while ((ovrst_r2 < lowerlimitOvrst) || (ovrst_r2 > upperlimitOvrst))
        stepNo_tnx = stepNo_tnx + 1;

        if (error2 > 0)                                                        %Me vasi ta diagrammata tis selidas 13, otan overshoot < 4.47(onomastiki timi)
            tnx_fuzzy = abs(tnx_fuzzy + (tnx_fuzzy*switchParam));              %isxiei tnx_fuzzy < tn_onomastiko. Ara to tnx_fuzzy tha prepei na afksithei

            plant.tnx = tnx_fuzzy;
            [ovrst_r2 Fcl_PI_fuzzy ti_2_fuzzy tnx_fuzzy tvx_fuzzy] = auto_tune_param_calculateOvs(plant);
            k = mod(stepNo_tnx,plotNo);
            if (k == 0)
               figure(figureIndex)
               step(Fcl_PI_fuzzy,'k')
               title('Fuzzy Tuning - PI Controller')
               hold on
            end
            error2 = ovPIdesired - ovrst_r2;
        elseif (error2 < 0)                                                    %Me vasi ta diagrammata tis selidas 13, otan overshoot > 4.47(onomastiki timi)
            tnx_fuzzy = abs(tnx_fuzzy - (tnx_fuzzy*switchParam));              %isxiei tnx_fuzzy > tn_onomastiko. Ara to tnx_fuzzy tha prepei na meiothei

            plant.tnx = tnx_fuzzy;
            [ovrst_r2 Fcl_PI_fuzzy ti_2_fuzzy tnx_fuzzy tvx_fuzzy] = auto_tune_param_calculateOvs(plant);
            k = mod(stepNo_tnx,plotNo);
            if (k == 0)
               figure(figureIndex)
               step(Fcl_PI_fuzzy,'k')
               title('Fuzzy Tuning - PI Controller')
               hold on
            end
            error2 = ovPIdesired - ovrst_r2;
                % relaxing the band
                % ----------------------------------------------------------------------------------------------
                if mod(stepNo_tvx,relaxCounter*100) == 0
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

        % Logging the ovs values
        % -------------------------------------------------------------
        ovs_PI_fuzzy_tune(stepNo_tnx) = ovrst_r2;

        n = mod(stepNo_tnx,1);
        if (n == 0)
          fprintf('step_tnx: %d - ti_fuzzy: %1.5f - tnx_fuzzy: %1.5f - tvx_fuzzy: %1.5f\n',stepNo_tnx,ti_2_fuzzy,tnx_fuzzy,tvx_fuzzy)
        end
    end
        if stepNo_tnx == 0
           stepNo_tnx = stepNo_tnx + 1;
           ovs_PI_fuzzy_tune(stepNo_tnx) = ovrst_r2;
        end

    % Logging Convergence Steps for PI Controller - Fuzzy Tuning (OVS Estimation)
    % *************************************************************************
    step_Conv_fuzzyTuning_complexZeros_PI_Controller = stepNo_tnx;

    % adding figure properties
    % *************************************************************************
    title('Fuzzy Automatic Tuning - PI Controller')      

else
    disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
    disp('No PI control exists....')
    disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
end

FlagOptimal_PI = 0;

% Calling the Optimal Controllers and Returning the Optimal TransFerFunctions   
% *************************************************************************
% ----------------------------------------------------
[Gc_pi Gp stepInformation_pi ti_MO_pi x_MO_pi y_MO_pi] = auto_tune_param_main_pi(plant);
if isnan(stepInformation_pi.Overshoot) == 1
    FlagOptimal_PI = 1                                   ;
else
    Gc_pi_optimal = Gc_pi                                ;
    Ffp_MO_PI_optimal = Gc_pi_optimal*Gp                 ;

    Fcl_MO_PI_optimal = feedback(Ffp_MO_PI_optimal,kh)   ;
    So_MO_PI_optimal  = 1 - Fcl_MO_PI_optimal            ; 
    Si_MO_PI_optimal  = series(So_MO_PI_optimal,Gp)      ;
end
    
% *************************************************************************
% *************************************************************************
fuzzy_tune_PI_log_parameters = zeros(3,2);
if LagFlagPI >= 0.5
    Fcl_fuzzy_PI_struct = stepinfo(Fcl_PI_fuzzy)          ;
    ovsFinal_Fcl_fuzzy_PI = Fcl_fuzzy_PI_struct.Overshoot ;
    ti_fuzzy = 2*kp*(tsx_fuzzy - tnx_fuzzy - tvx_fuzzy)   ;
    x_fuzzy  = tnx_fuzzy + tvx_fuzzy                      ;
    y_fuzzy  = tnx_fuzzy * tvx_fuzzy                      ;

    fuzzy_tune_PI_log_parameters(1,1) = ti_fuzzy ;  fuzzy_tune_PI_log_parameters(1,2) = ti_fuzzy;
    fuzzy_tune_PI_log_parameters(2,1) = tnx_fuzzy;  fuzzy_tune_PI_log_parameters(2,2) = x_fuzzy ;
    fuzzy_tune_PI_log_parameters(3,1) = tvx_fuzzy;  fuzzy_tune_PI_log_parameters(3,2) = y_fuzzy ;
end

disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
disp('++++++++++++++++++   Closed Loop Control Information ++++++++++++++')
disp('++++++++++++++++++          PI Controller            ++++++++++++++')
disp('++++++++++++++++++          Fuzzy Tuning             ++++++++++++++')
if LagFlagPI >= 0.5
    fprintf('ti_fuzzy : %1.5f - tnx_fuzzy: %1.5f - tvx_fuzzy: %1.5f\n',ti_fuzzy,tnx_fuzzy,tvx_fuzzy)
    fprintf('ti_fuzzy : %1.5f - x_fuzzy  : %1.5f - y_fuzzy  : %1.5f\n',ti_fuzzy,x_fuzzy,y_fuzzy)
else
    fprintf('No PI Control Exists (based on the Fuzzy model)\n')
end
if FlagOptimal_PI == 0
    fprintf('ti_optim : %1.5f - x_optim  : %1.5f - y_optim  : %1.5f\n',ti_MO_pi,x_MO_pi,y_MO_pi)
else
    fprintf('No PI Optimal Control Exists\n')
end
if LagFlagPI >= 0.5
    fprintf('PI_fuzzy_tuned_ovs: %1.5f%%\n',ovsFinal_Fcl_fuzzy_PI)
    disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
end

% Initialize Log Matrix
% ------------------------------------------------------
log_ovs_while_fuzzy_pi_tuning = zeros(5,1);

if LagFlagPI >= 0.5
    log_ovs_while_fuzzy_pi_tuning(1) = lowerlimitOvrst;
    log_ovs_while_fuzzy_pi_tuning(2) = upperlimitOvrst;
    log_ovs_while_fuzzy_pi_tuning(3) = ovsFinal_Fcl_fuzzy_PI ; % ...where it finally tuned
    log_ovs_while_fuzzy_pi_tuning(4) = ovPIdesired           ; % ...desired (fromt the fuzzy estimator)
end
if FlagOptimal_PI == 0
    log_ovs_while_fuzzy_pi_tuning(5) = ovPIoptimal           ; % ...optimal
end

disp('press any key to continue...')
pause


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% ######################        FUZZY TUNING   ############################
% ######################        PID Control    ############################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
tr_choice = input('Train fuzzy PID?: y/n:','s');

if  strcmp(tr_choice,'y') == 1
      
    % Load the final with the training inputs
    % (ti_I ov_PI ti_PI Settling Time op-loop) - desiredPID ovs
    % ---------------------------------------------------------------------
    load ov_pid_zeros.txt

    % **********************  CLASSICICATION PROBLEM  *********************
    % *********************************************************************
    % Identify whether or not the system can be
    % controlled by a PI controller
    % *********************************************************************
    ov_pid_control = ov_pid_zeros;

    % Get the number of examples
    % *********************************************************************
    numExamples = size(ov_pid_zeros);

    % Formulate the training data set
    % *********************************************************************
    index = 1;
    k = 1    ;

    Xtrain_pid_control = ov_pid_control;

    for i = 1:numExamples
        if ((Xtrain_pid_control(i,5) < 3)||(Xtrain_pid_control(i,5) > 8.5))
            Xtrain_pid_control(i,5) = 0;
        else
            Xtrain_pid_control(i,5) = 1;
        end
    end

    % Use the Substractive Clustering method to generate a new Sugeno 
    % type fuzzy system
    % *********************************************************************
    disp('Training the PID control classifier ....')

    % Grid partition with Gaussian 3 MF and linear type consequent
    % fismat=genfis1(Xtrain_control(:,1:5),[3 3 3 3],'gaussmf','linear');

    % We use the default values of Jiang for the method, fismat contains
    % the initial fuzzy system
    % *********************************************************************
    fismat = genfis2(Xtrain_pid_control(:,1:4),Xtrain_pid_control(:,5),0.2); 

    % Train the system using the anfis editor
    % *********************************************************************
    [fismat_control_pid,error] = anfis(Xtrain_pid_control,fismat,[100],[0 0 0 0]);

    % Calculate the ouput of ANFIS for the training data set
    % *********************************************************************
    trn_out_fismat_control_pid = evalfis(Xtrain_pid_control(:,1:4),fismat_control_pid);


    trn_out_fismat_control_pid = round(trn_out_fismat_control_pid);
    correct = 0;
    for i = 1:numExamples(1)
        if (trn_out_fismat_control_pid(i) == Xtrain_pid_control(i,5))
           correct = correct + 1;
        end
    end

    Percentage = correct/numExamples(1)*100    ;
    fprintf('PID Classification Percentage = %3.f\n',Percentage)
   
    
    % **********************   PREDICTION PROBLEM *************************
    % *********************************************************************
    % Predict the desired PID overshoot when the plant
    % can be controlled by a PID controller
    % *********************************************************************
    ov_pid_pred = ov_pid_zeros;

    % Get the number of examples
    % *********************************************************************
    numExamples = size(ov_pid_zeros);

    % Formulate the training data set
    % *********************************************************************
    index = 1;
    k = 1    ;

    for i=1:numExamples
        if (Xtrain_pid_control(i,5) == 1)
            Xtrain_pred_pid(index,:) = ov_pid_pred(i,:);
            index = index + 1;
        end 
    end

    % Use the Substractive Clustering method to generate a new Sugeno type 
    % fuzzy system
    % *********************************************************************
    disp('Training fuzzy PID overshoot estimator....')

    % Grid partition with Gaussian 3 MF and linear type consequent
    % fismat=genfis1(Xtrain_pred_pid(:,1:3),[3 3],'gaussmf','linear');

    % We use the default values of Jiang for the method, fismat contains 
    % the initial fuzzy system
    % *********************************************************************
    fismat = genfis2(Xtrain_pred_pid(:,1:4),Xtrain_pred_pid(:,5),0.2);

    % Train the system using the anfis editor
    % *********************************************************************
    [fismat_pred_pid,error1] = anfis(Xtrain_pred_pid,fismat,[100],[0 0 0 0]);

else
    % Load the already trained fuzzy PI overshoot estimator
    % ---------------------------------------------------------------------
    load ov_pid_zeros.txt
    load fismat_PID
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                   Fuzzy PID Overshoot Estimation
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% SOS!!!!!: you have to measure the ACTUAL overshoot after the PI Fuzzy 
% tuning and NOT the desired overshoot generated by the ANFIS. In that, let
% the PI controller be fuzzy tuned, and then measure the OVS from the
% closed loop transfer function. There might be a small deviation beteween
% the actual ovs and the desired ovs.

Fcl_local_PI_fuzzyStruct = stepinfo(Fcl_PI_fuzzy);
actualOvs_Fuzzy_PI = Fcl_local_PI_fuzzyStruct.Overshoot;

% Check the inputs for extreem values as inputs to fuzzy for PID!!!!!!
% ti from I control
% ---------------------------------------------------------------------

[min_pid(1),min_ind(1)] = min(ov_pid_zeros(:,1));
[min_pid(2),min_ind(2)] = min(ov_pid_zeros(:,2));
[min_pid(3),min_ind(3)] = min(ov_pid_zeros(:,3));
[min_pid(4),min_ind(4)] = min(ov_pid_zeros(:,4));
[max_pid(1),max_ind(1)] = max(ov_pid_zeros(:,1));
[max_pid(2),max_ind(2)] = max(ov_pid_zeros(:,2));
[max_pid(3),max_ind(3)] = max(ov_pid_zeros(:,3));
[max_pid(4),max_ind(4)] = max(ov_pid_zeros(:,4));

if (ti_auto_fuzzyInput/kp < min_pid(1))
    fuzzyPID(1) = min_pid(1);
elseif (ti_auto_fuzzyInput/kp > max_pid(1))
    fuzzyPID(1) = max_pid(1);
else
    fuzzyPID(1) = ti_auto_fuzzyInput/kp;
end
   
% ovPI from PI control
% ---------------------------------------------------------------------

if (actualOvs_Fuzzy_PI < min_pid(2))
    fuzzyPID(2) = min_pid(2);
elseif (actualOvs_Fuzzy_PI > max_pid(2))
    fuzzyPID(2) = max_pid(2);
else
    fuzzyPID(2) = actualOvs_Fuzzy_PI;
end
   
% ti from PI control
% ---------------------------------------------------------------------
if ((2*kp*(tsx_fuzzy - tnx_fuzzy))/kp < min_pid(3))
     fuzzyPID(3) = min_pid(3);
elseif ((2*kp*(tsx_fuzzy - tnx_fuzzy))/kp > max_pid(3))
     fuzzyPID(3) = max_pid(3);
else
    fuzzyPID(3) = (2*kp*(tsx_fuzzy - tnx_fuzzy))/kp;
end

% Open Loop Settling Time
% ---------------------------------------------------------------------
if (plant_tss < min_pid(4))
    fuzzyPID(4) = min_pid(4);
elseif (plant_tss > max_pid(4))
    fuzzyPID(4) = max_pid(4);
else
    fuzzyPID(4) = plant_tss ;
end


LagFlagPID = evalfis([fuzzyPID(1) fuzzyPID(2) fuzzyPID(3) fuzzyPID(4)],fismat_control_pid);
% LagFlagPID = round(LagFlagPID);

if (LagFlagPID >= 0.5)
    % return the desired ovs for the PID Controller
    % *********************************************************************
    ovPIDdesired = evalfis([fuzzyPID(1) fuzzyPID(2) fuzzyPID(3) fuzzyPID(4)],fismat_pred_pid); 
else
    disp('No PID control...')
    disp('Switching to PI-LAG is not ready yet...')
    ovPIDdesired = 'No PID control required...';
end

% ovPIDdesired = evalfis([fuzzyPID(1) fuzzyPID(2) fuzzyPID(3) fuzzyPID(4)],fismat_PID); % return the desired ovs for the PI Controller
% ovPIDdesired = evalfis([actualOvs_Fuzzy_PI 2*kp*(tsx_fuzzy - tnx_fuzzy) plant_tss],fismat_PID); % return the desired ovs for the PI Controller

% % Calling the Optimal Controllers and Returning the Optimal TransFerFunctions   
% % ***************************************************************************
FlagOptimal_PID = 0;

[Gc_pid Gp stepInformation_pid ti_MO_pid x_MO_pid y_MO_pid] = auto_tune_param_main_pid(plant);
if isnan(stepInformation_pid.Overshoot) == 1
    FlagOptimal_PID = 1;
else
    Gc_pid_optimal = Gc_pid                              ;
    Ffp_MO_PID_optimal = Gc_pid_optimal*Gp               ;
    
    Fcl_MO_PID_optimal = feedback(Ffp_MO_PID_optimal,kh) ;
    structoptimal = stepinfo(Fcl_MO_PID_optimal);
    ovsoptimal = structoptimal.Overshoot                 ;
    ovPIDoptimal = ovsoptimal                            ; 
end

% -------------------------------------------------------------------------
% Y - Prediction using ANFIS based on x_PI and x_PID
% -------------------------------------------------------------------------
tr_choice = input('Train y fuzzy estimation for PID?: y/n:','s');

if  strcmp(tr_choice,'y') == 1
    
    %Inputs x_PI x_PID - output y_estimation
    % ---------------------------------------------------------------------
    load y_pred.txt

    % Get the number of examples
    % ---------------------------------------------------------------------
    numExamples = size(y_pred);

    %Formulate the training data set
    % ---------------------------------------------------------------------
    index = 1;
    k = 1    ;

    % First we use the Tsc=0.01,0.11,0.12...0.11 values and a=0.1,0.12,0.14,...0.9 for
    % each Tsc to formulate the training data set
    % ---------------------------------------------------------------------

    for i = 1:1:91*9;
        Xtrain_ypred(k,:) = y_pred(index,:);
        k = k + 1;
        index = index + 1;
    end

    % Then we use the Tsc = 0.21,0.31,...0.81 values and a=0.1,0.12,0.14,...0.9 for
    % each Tsc to formulate the training data set
    % ---------------------------------------------------------------------
    while (index < numExamples(1))
        for i = 1:1:91
            Xtrain_ypred(k,:) = y_pred(index,:);
            k = k + 1           ;
            index = index + 1   ;
        end
            index = index + 91*9;
    end


    [min_y(1),min_ind(1)] = min(y_pred(:,1));
    [min_y(2),min_ind(2)] = min(y_pred(:,2));
    [max_y(1),max_ind(1)] = max(y_pred(:,1));
    [max_y(2),max_ind(2)] = max(y_pred(:,2));
    Xtrain_ypred(k,:) = y_pred(min_ind(1),:);

    k = k + 1;
    Xtrain_ypred(k,:) = y_pred(min_ind(2),:);

    k = k + 1;
    Xtrain_ypred(k,:) = y_pred(max_ind(1),:);

    k = k + 1;
    Xtrain_ypred(k,:) = y_pred(max_ind(2),:);


    % Use the Substractive Clustering method to generate a new Sugeno 
    % type fuzzy system
    % ---------------------------------------------------------------------
    disp('Training fuzzy y for PID estimator....')

    %Grid partition with Gaussian 3 MF and linear type consequent
    % ---------------------------------------------------------------------
    fismat = genfis1(Xtrain_ypred(:,1:3),[3 3],'gaussmf','linear');

    % fismat = genfis2(Xtrain_ypred(:,1:2),Xtrain_ypred(:,3),0.5); 
    % We use the default values of Jiang for the method, fismat contains 
    % the initial fuzzy system

    % Train the system using the anfis editor
    % ---------------------------------------------------------------------
    [fismat_ypred,error1] = anfis(Xtrain_ypred,fismat,[100],[0 0 0 0]);

    % Calculate the ouput of ANFIS for the training data set
    % ---------------------------------------------------------------------
    trn_out_fismat_ypred = evalfis(Xtrain_ypred(:,1:2),fismat_ypred);

    %Plot ANFIS training output versus desired output
    % ---------------------------------------------------------------------
    figureIndex = figureIndex + 1;
    figure(figureIndex)
    plot(Xtrain_ypred(:,3),'o')
    hold on;
    plot(trn_out_fismat_ypred,'r*')
    title(['y estimation Training Examples - SSE = ',num2str(error1(end))])
    ylabel('y')
    xlabel('Example')
    grid on
    save fismat_ypred fismat_ypred y_pred
end

disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('+++++++++++++++++   Fuzzy PID Overshoot Estimation  +++++++++++++++')
disp('+++++++++++++++++          PID  Controller          +++++++++++++++')
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
if FlagOptimal_PID == 1
    fprintf('No OPTIMAL PID Control exists...\n')
else
    fprintf('PID optimal OVS: %1.5f\n',ovPIDoptimal)
end

if LagFlagPID >= 0.5
    fprintf('PID desired OVS: %1.5f\n',ovPIDdesired)
    disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
        % tvx_fuzzy : tuning
        % -------------------------------------------------------------------------
        % Initial Conditions
        % -------------------------------------------------------------------------
        stepNo_tvx = 0                       ;
        x_fuzzy = tnx_fuzzy + (tnx_fuzzy/2)  ;
        % x_fuzzy = 0.05*(tnx_fuzzy + (tnx_fuzzy/2))  ;
        plant.tsx = tsx_fuzzy                ;  plant.x = x_fuzzy;  plant.x2 = tnx_fuzzy;
        figureIndex = figureIndex + 1        ;
        [ovrst_r3 Fcl_PID_fuzzy ti_3_fuzzy x_fuzzy y_fuzzy] = auto_tune_param_calculateOvs_x_y(plant);

        error3 = ovPIDdesired - ovrst_r3     ;
        lowerlimitOvrst = ovPIDdesired * 0.99;
        upperlimitOvrst = ovPIDdesired * 1.01;
        relaxCounter = 1;
        while ((ovrst_r3 < lowerlimitOvrst) || (ovrst_r3 > upperlimitOvrst))
            stepNo_tvx = stepNo_tvx + 1;

            if (error3 > 0)                                                  %Me vasi ta diagrammata tis selidas 13, otan overshoot < 4.47(onomastiki timi)
                x_fuzzy = abs(x_fuzzy + (x_fuzzy*switchParam));              %isxiei tvx < tv_onomastiko. Ara to tvx tha prepei na afksithei
                plant.x = x_fuzzy;
                [ovrst_r3 Fcl_PID_fuzzy ti_3_fuzzy x_fuzzy y_fuzzy] = auto_tune_param_calculateOvs_x_y(plant);
                k = mod(stepNo_tvx,plotNo);
                if (k == 0)
                   figure(figureIndex)
                   step(Fcl_PID_fuzzy,'k')
                   title('Fuzzy Tuning - PID Controller')
                   hold on
                end
                error3 = ovPIDdesired - ovrst_r3;
            elseif (error3 < 0)                                               %Me vasi ta diagrammata tis selidas 13, otan overshoot > 4.47(onomastiki timi)
                x_fuzzy = abs(x_fuzzy - (x_fuzzy*switchParam));               %isxiei tvx > tv_onomastiko. Ara to tvx tha prepei na meiothei
                plant.x = x_fuzzy;
                [ovrst_r3 Fcl_PID_fuzzy ti_3_fuzzy x_fuzzy y_fuzzy] = auto_tune_param_calculateOvs_x_y(plant);
                k = mod(stepNo_tvx,plotNo);
                if (k == 0)
                   figure(figureIndex)
                   step(Fcl_PID_fuzzy,'k')
                   title('Fuzzy Tuning - PID Controller')
                   hold on
                end
                error3 = ovPIDdesired - ovrst_r3;
                    % relaxing the band
                    % ----------------------------------------------------------------------------------------------
                    if mod(stepNo_tvx,relaxCounter*100) == 0
                       relaxCounter = relaxCounter + 1;
                        % widen the band
                        % ----------------------------------------------------------------------------------------------
                        upperlimitOvrst = upperlimitOvrst + 0.01*upperlimitOvrst; 
                        lowerlimitOvrst = lowerlimitOvrst - 0.01*lowerlimitOvrst;
                        disp('----------------------------------------------------------------------');
                        fprintf('Relaxing the Reference Band so that the PID controller tuning....\n')
                        fprintf('converges faster...\n')
                        fprintf('New Upper Limit: %2.5f - New Lower Limit: %2.5f \n',upperlimitOvrst,lowerlimitOvrst)
                        disp('----------------------------------------------------------------------');
                    end
            end

            % Logging the ovs values
            % -------------------------------------------------------------
            ovs_PID_fuzzy_tune(stepNo_tvx) = ovrst_r3;

            n = mod(stepNo_tvx,1);
            if (n == 0)
              fprintf('step_tvx: %d - ti_fuzzy: %1.5f - x_fuzzy: %1.5f - y_fuzzy: %1.5f \n',stepNo_tvx,ti_3_fuzzy,x_fuzzy,y_fuzzy)
            end
        end
            if stepNo_tvx == 0
               stepNo_tvx = stepNo_tvx + 1;
               ovs_PID_fuzzy_tune(stepNo_tvx) = ovrst_r3;
            end

        % Logging the overshoot while tuning
        % *************************************************************************
        ovs_Final_fuzzy_tune = [ovs_I_auto_tune ovs_PI_fuzzy_tune ovs_PID_fuzzy_tune];

        % Logging Convergence Steps for PID Controller - Fuzzy Tuning (OVS Estimation)
        % *************************************************************************
        step_Conv_fuzzyTuning_complexZeros_PID_Controller = stepNo_tvx;

        % adding figure properties
        % *************************************************************************
        title('Fuzzy Automatic Tuning - PID Controller')      
        grid on

        % *************************************************************************
        % *************************************************************************
        Fcl_fuzzy_PID_struct = stepinfo(Fcl_PID_fuzzy)         ;
        ovsFinal_Fcl_fuzzy_PID = Fcl_fuzzy_PID_struct.Overshoot;
        ti_fuzzy = 2*kp*(tsx_fuzzy - x_fuzzy)                  ;

        fuzzy_tune_PID_log_parameters(1,1) = ti_fuzzy;
        fuzzy_tune_PID_log_parameters(2,1) = x_fuzzy ;
        fuzzy_tune_PID_log_parameters(3,1) = y_fuzzy ;
        % *************************************************************************
        % *************************************************************************
else
    fprintf('PID desired OVS does not exist...\n')
    disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
    
end     % OF LagFlagPID >= 0.5

% Calling the optimal Function for the PID Controller
% *************************************************************************
[Gc_pid Gp stepInformation_pid ti_MO_pid x_MO_pid y_MO_pid] = auto_tune_param_main_pid(plant);
Gc_pid_optimal = Gc_pid                              ;
Ffp_MO_PID_optimal = Gc_pid_optimal*Gp               ;
Fcl_MO_PID_optimal = feedback(Ffp_MO_PID_optimal,kh) ;

ti_optim = ti_MO_pid;    % storing the values
x_optim  = x_MO_pid ;
y_optim  = y_MO_pid ; 

% *************************************************************************
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
disp('++++++++++++++++++   Closed Loop Control Information ++++++++++++++')
disp('++++++++++++++++++         PID Controller            ++++++++++++++')
disp('++++++++++++++++++          Fuzzy Tuning             ++++++++++++++')
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
if LagFlagPID >= 0.5
    fprintf('ti_fuzzy     : %1.5f - x_fuzzy      : %1.5f - y_fuzzy      : %1.5f\n',ti_fuzzy,x_fuzzy,y_fuzzy)
    fprintf('PID_fuzzy_tuned_ovs: %1.5f%%\n',ovsFinal_Fcl_fuzzy_PID)
    fprintf('PID_fuzzy [desired ovs]: %1.5f%%\n',ovPIDdesired)
    disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
else
    fprintf('NO Fuzzy PID Control Exists....\n')
end
if FlagOptimal_PID == 0
    fprintf('ti_optim     : %1.5f - x_optim      : %1.5f - y_optim      : %1.5f\n',ti_optim,x_optim,y_optim)
else
    fprintf('NO OPTIMAL PID Control Exists....\n')
end

fprintf('tp1 : %1.5f \t tz1 : %1.5f \n',tp1,tz1)
fprintf('tp2 : %1.5f \t tz2 : %1.5f \n',tp2,tz2)
fprintf('tp3 : %1.5f \t tz3 : %1.5f \n',tp3,tz3)
fprintf('tp4 : %1.5f \t tz4 : %1.5f \n',tp4,tz4)
fprintf('tp5 : %1.5f \t tp6 [Controller Unmodelled dynamic]: %1.5f \n',tp5,tp6)
fprintf('kp  : %1.5f \t kh : %1.5f \n',kp,kh)

% Initialize Log Matrix
% ------------------------------------------------------
log_ovs_while_fuzzy_pid_tuning = zeros(5,1);
if LagFlagPID >= 0.5
   log_ovs_while_fuzzy_pid_tuning(1) = lowerlimitOvrst;
   log_ovs_while_fuzzy_pid_tuning(2) = upperlimitOvrst;
   log_ovs_while_fuzzy_pid_tuning(3) = ovsFinal_Fcl_fuzzy_PID ; % ...where it finally tuned
   log_ovs_while_fuzzy_pid_tuning(4) = ovPIDdesired           ; % ...desired (fromt the fuzzy estimator)
end
if FlagOptimal_PID == 0
   log_ovs_while_fuzzy_pid_tuning(5) = ovPIDoptimal           ; % ...optimal
end
 


checkExistence = exist('commandLogInformation.txt','file');
if checkExistence == 2
    %delete old Log File
    delete 'commandLogInformation.txt'
end

diary('commandLogInformation.txt')

% ##############################################################################################
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    Convergence Results  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% ----------------------------------------------------------------------------------------------
% ##############################################################################################
% ----------------------------------  AUTOMATIC TUNING  ----------------------------------------
% ----------------------------------------------------------------------------------------------

fprintf('##################################################################################\n')
fprintf('###########################  CONVERGENCE RESULTS #################################\n')
fprintf('##################################################################################\n')

fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n')
fprintf('Automatic Tuning [Real Poles Controller] - I   Controller tuned after : %2.1f \n',step_Conv_autoTuning_realZeros_I_Controller)
fprintf('OVS: Fcl Response : %1.5f \n',ovsFinal_Fcl_auto_I)

fprintf('Automatic Tuning [Real Poles Controller] - PI  Controller tuned after : %2.1f \n',step_Conv_autoTuning_realZeros_PI_Controller)
fprintf('OVS: Fcl Response : %1.5f \n',ovsFinal_Fcl_auto_PI)

fprintf('Automatic Tuning [Real Poles Controller] - PID Controller tuned after : %2.1f \n',step_Conv_autoTuning_realZeros_PID_Controller)
fprintf('OVS: Fcl Response : %1.5f \n',ovsFinal_Fcl_auto_PID)
fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n')


fprintf('Automatic Tuning [Complex Poles Controller] - I   Controller tuned after : %2.1f \n',step_Conv_autoTuning_realZeros_I_Controller)
fprintf('OVS: Fcl Response : %1.5f \n',ovsFinal_Fcl_auto_I)

fprintf('Automatic Tuning [Complex Poles Controller] - PI  Controller tuned after : %2.1f \n',step_Conv_autoTuning_realZeros_PI_Controller)
fprintf('OVS: Fcl Response : %1.5f \n',ovsFinal_Fcl_auto_PI)

fprintf('Automatic Tuning [Complex Poles Controller] - PID Controller tuned after : %2.1f \n',step_Conv_autoTuning_complexZeros_PID_Controller)
fprintf('OVS: Fcl Response : %1.5f \n',ovsFinal_Fcl_auto_marg_PID)
fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n')


fprintf('Fuzzy Tuning - I   Controller tuned after : %1.5f \n',step_Conv_autoTuning_realZeros_I_Controller)
fprintf('OVS: Fcl Response : %1.5f \n',ovsFinal_Fcl_auto_I)

if LagFlagPI >= 0.5
    fprintf('Fuzzy Tuning - PI  Controller tuned after : %1.5f \n',step_Conv_fuzzyTuning_complexZeros_PI_Controller)
    fprintf('OVS: Fcl Response : %1.5f \n',ovsFinal_Fcl_auto_PI)
else
    fprintf('No PI Control exists [based on fuzzy estimation]... \n')
end


if LagFlagPID >= 0.5
    fprintf('Fuzzy Tuning - PID Controller tuned after : %1.5f \n',step_Conv_fuzzyTuning_complexZeros_PID_Controller)
    fprintf('OVS: Fcl Response : %1.5f \n',ovsFinal_Fcl_auto_marg_PID)
    fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n')
else
    fprintf('No PID Control exists [based on fuzzy estimation]... \n')
end


fprintf('##################################################################################\n')
fprintf('##################################################################################\n')

fprintf('\n')
fprintf('##################################################################################\n')
fprintf('###########################  PARAMETER TUNING RESULTS ############################\n')
fprintf('##################################################################################\n')
fprintf('1. Automatic Tuning based on ovs:4.47 %%\n')
fprintf('------------------------------------------------------------------------------------\n')
fprintf('\n')


fprintf('++++++++++++++++++   Closed Loop Control Information ++++++++++++++\n')
fprintf('++++++++++++++++++           I Controller            ++++++++++++++\n')
fprintf('ti_auto : %1.5f - tnx_auto: %1.5f - tvx_auto: %1.5f\n',auto_tune_I_log_parameters(1,1),auto_tune_I_log_parameters(2,1),auto_tune_I_log_parameters(3,1))
fprintf('ti_auto : %1.5f - x_auto  : %1.5f - y_auto  : %1.5f\n',auto_tune_I_log_parameters(1,2),auto_tune_I_log_parameters(2,2),auto_tune_I_log_parameters(3,2))
fprintf('I_auto_tuned_ovs: %1.5f%%\n',ovsFinal_Fcl_auto_I)
fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n')


fprintf('++++++++++++++++++   Closed Loop Control Information ++++++++++++++\n')
fprintf('++++++++++++++++++           PI Controller           ++++++++++++++\n')
fprintf('ti_auto : %1.5f - tnx_auto: %1.5f - tvx_auto: %1.5f\n',auto_tune_PI_log_parameters(1,1),auto_tune_PI_log_parameters(2,1),auto_tune_PI_log_parameters(3,1))
fprintf('ti_auto : %1.5f - x_auto  : %1.5f - y_auto  : %1.5f\n',auto_tune_PI_log_parameters(1,2),auto_tune_PI_log_parameters(2,2),auto_tune_PI_log_parameters(3,2))
fprintf('PI_auto_tuned_ovs: %1.5f%%\n',ovsFinal_Fcl_auto_PI)
fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n')


fprintf('++++++++++++++++++   Closed Loop Control Information ++++++++++++++\n')
fprintf('++++++++++++++++++          PID Controller           ++++++++++++++\n')
fprintf('ti_auto : %1.5f - tnx_auto: %1.5f - tvx_auto: %1.5f\n',auto_tune_PID_log_parameters(1,1),auto_tune_PID_log_parameters(2,1),auto_tune_PID_log_parameters(3,1))
fprintf('ti_auto : %1.5f - x_auto  : %1.5f - y_auto  : %1.5f\n',auto_tune_PID_log_parameters(1,2),auto_tune_PID_log_parameters(2,2),auto_tune_PID_log_parameters(3,2))
fprintf('PID_auto_tuned_ovs: %1.5f%%\n',ovsFinal_Fcl_auto_PID)
fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n')


% ###############################################################################################
fprintf('\n')
fprintf('2. Automatic Tuning based on ovs:4.47 \n')
fprintf('[Complex Zeros Controller with Y estimation]:\n')
fprintf('------------------------------------------------------------------------------------\n')
fprintf('\n')

fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')
fprintf('++++++++++++++++++   Closed Loop Control Information ++++++++++++++\n')
fprintf('++++++++++++++++++          PID Controller           ++++++++++++++\n')
fprintf('++++++++++++++++++  Automatic Tuning: "Y" Estimation ++++++++++++++\n')
fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')
fprintf('ti_auto_marg : %1.5f - x_auto_marg: %1.5f - y_auto_marg: %1.5f\n',ti_3_auto_marg,x_auto_marg,y_auto_marg)
fprintf('PI_auto_tuned_marg_ovs: %1.5f%%\n',ovsFinal_Fcl_auto_marg_PID)
fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')

% ###############################################################################################
fprintf('\n')
fprintf('3. Fuzzy Automatic Tuning based\n')
fprintf('------------------------------------------------------------------------------------\n')
fprintf('\n')

fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('++++++++++++++++++   Fuzzy PI Overshoot Estimation  +++++++++++++++\n')
fprintf('++++++++++++++++++          PI  Controller          +++++++++++++++\n')
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('PI optimal OVS: %1.5f\n',ovPIoptimal)
fprintf('PI desired OVS: %1.5f\n',ovPIdesired)
fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n')

fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('+++++++++++++++++   Fuzzy PID Overshoot Estimation  +++++++++++++++\n')
fprintf('+++++++++++++++++          PID  Controller          +++++++++++++++\n')
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('PID optimal OVS: %1.5f\n',ovPIDoptimal)
if (LagFlagPID >= 0.5)
    fprintf('PID desired OVS: %1.5f\n',ovPIDdesired)
else
    fprintf('PID desired OVS: %s\n',ovPIDdesired)
end
fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n')


fprintf('++++++++++++++++++   Closed Loop Control Information ++++++++++++++\n')
fprintf('++++++++++++++++++          PI Controller            ++++++++++++++\n')
fprintf('++++++++++++++++++          Fuzzy Tuning             ++++++++++++++\n')
if LagFlagPI >= 0.5
    fprintf('ti_fuzzy : %1.5f - tnx_fuzzy: %1.5f - tvx_fuzzy: %1.5f\n',fuzzy_tune_PI_log_parameters(1,1),fuzzy_tune_PI_log_parameters(2,1),fuzzy_tune_PI_log_parameters(3,1))
    fprintf('ti_fuzzy : %1.5f - x_fuzzy  : %1.5f - y_fuzzy  : %1.5f\n',fuzzy_tune_PI_log_parameters(1,2),fuzzy_tune_PI_log_parameters(2,2),fuzzy_tune_PI_log_parameters(3,2))
    fprintf('PI_fuzzy_tuned_ovs: %1.5f%%\n',ovsFinal_Fcl_fuzzy_PI)
    fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n')
else
    fprintf('NO PI Control exists [based on fuzzy estimation]...\n')
end

if FlagOptimal_PI == 0
    fprintf('ti_optim : %1.5f - x_optim  : %1.5f - y_optim  : %1.5f\n',ti_MO_pi,x_MO_pi,y_MO_pi)
else
    fprintf('NO Optimal PI Control exists...\n')
end

fprintf('++++++++++++++++++   Closed Loop Control Information ++++++++++++++\n')
fprintf('++++++++++++++++++         PID Controller            ++++++++++++++\n')
fprintf('++++++++++++++++++          Fuzzy Tuning             ++++++++++++++\n')
if LagFlagPID >= 0.5
    fprintf('ti_fuzzy     : %1.5f - x_fuzzy      : %1.5f - y_fuzzy      : %1.5f\n',fuzzy_tune_PID_log_parameters(1,1),fuzzy_tune_PID_log_parameters(2,1),fuzzy_tune_PID_log_parameters(3,1))
    fprintf('ti_optim     : %1.5f - x_optim      : %1.5f - y_optim      : %1.5f\n',ti_optim,x_optim,y_optim)
    fprintf('PID_fuzzy_tuned_ovs: %1.5f%%\n',ovsFinal_Fcl_fuzzy_PID)
else
    fprintf('NO PID Control exists [based on fuzzy estimation]...\n')
end


% ###############################################################################################
fprintf('\n')
fprintf('3. Optimal Tuning\n')
fprintf('------------------------------------------------------------------------------------')
fprintf('\n')

fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n')
fprintf('Optimal Tuning - I   Controller \n')
fprintf('OVS: Fcl Response : %1.5f \n',Fcl_MO_I_optimal_ovs)
if FlagOptimal_PI == 0
    fprintf('Optimal Tuning - PI   Controller \n')
    fprintf('OVS: Fcl Response : %1.5f \n',ovPIoptimal)
else
    fprintf('NO Optimal PI Control exists...\n')
end

if FlagOptimal_PID == 0
    fprintf('Optimal Tuning - PID   Controller \n')
    fprintf('OVS: Fcl Response : %1.5f \n',ovPIDoptimal)
    fprintf('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n')
else
    fprintf('NO Optimal PID Control exists...\n')
end

% ###############################################################################################

diary off

% *************************************************************************
% Restoring the Autotuned Transfer Functions for the I, PI, PID Controller   
% *************************************************************************
Fcl_autoTuned_I = Fcl_auto_I                         ;
So_MO_I_autotuned  = 1 - Fcl_autoTuned_I             ; 
Si_MO_I_autotuned  = series(So_MO_I_autotuned,Gp)    ;
 
Fcl_autoTuned_PI = Fcl_auto_PI                       ;
So_MO_PI_autotuned  = 1 - Fcl_autoTuned_PI           ; 
Si_MO_PI_autotuned  = series(So_MO_PI_autotuned,Gp)  ;

Fcl_autoTuned_PID = Fcl_auto_PID                     ; 
So_MO_PID_autotuned  = 1 - Fcl_autoTuned_PID         ; 
Si_MO_PID_autotuned  = series(So_MO_PID_autotuned,Gp);

Fcl_PID_auto_marg = Fcl_PID_auto_marg                ;
So_MO_PID_auto_marg = 1 - Fcl_PID_auto_marg          ;
Si_MO_PID_auto_marg = series(So_MO_PID_auto_marg,Gp) ;

% Calling the Optimal Controllers and Returning the Optimal TransFerFunctions   
% *************************************************************************
[Gc_pid Gp stepInformation_pid ti_MO_pid x_MO_pid y_MO_pid] = auto_tune_param_main_pid(plant);
Gc_pid_optimal = Gc_pid                              ;
Ffp_MO_PID_optimal = Gc_pid_optimal*Gp               ;

Fcl_MO_PID_optimal = feedback(Ffp_MO_PID_optimal,kh) ;
So_MO_PID_optimal  = 1 - Fcl_MO_PID_optimal          ; 
Si_MO_PID_optimal  = series(So_MO_PID_optimal,Gp)    ;
% ----------------------------------------------------
[Gc_pi Gp stepInformation_pi ti_MO_pi x_MO_pi] = auto_tune_param_main_pi(plant);
Gc_pi_optimal = Gc_pi                                ;
Ffp_MO_PI_optimal = Gc_pi_optimal*Gp                 ;

Fcl_MO_PI_optimal = feedback(Ffp_MO_PI_optimal,kh)   ;
So_MO_PI_optimal  = 1 - Fcl_MO_PI_optimal            ; 
Si_MO_PI_optimal  = series(So_MO_PI_optimal,Gp)      ;
% ----------------------------------------------------
[Gc_i Gp stepInformation_i ti_MO_i] = auto_tune_param_main_i(plant)    ;
Gc_i_optimal = Gc_i                                  ;
Ffp_MO_I_optimal = Gc_i_optimal*Gp                   ;

Fcl_MO_I_optimal = feedback(Ffp_MO_I_optimal,kh)     ;
So_MO_I_optimal  = 1 - Fcl_MO_I_optimal              ; 
Si_MO_I_optimal  = series(So_MO_I_optimal,Gp)        ;

% ----------------------------------------------------
% Restoring the Fuzzy Autotuned Transfer Functions for the I, PI, PID Controller   
% *************************************************************************
Fcl_I_fuzzy = Fcl_auto_I;
Fcl_fuzzy_autoTuned_I = Fcl_I_fuzzy                              ; % Note that the I autotuned Fcl is the same
So_MO_I_fuzzy_autotuned  = 1 - Fcl_fuzzy_autoTuned_I             ; % with the I Fuzzy autotuned Fcl
Si_MO_I_fuzzy_autotuned  = series(So_MO_I_fuzzy_autotuned,Gp)    ;


Fcl_fuzzy_autoTuned_PI = Fcl_PI_fuzzy                            ;
So_MO_PI_fuzzy_autotuned  = 1 - Fcl_fuzzy_autoTuned_PI           ;  
Si_MO_PI_fuzzy_autotuned  = series(So_MO_PI_fuzzy_autotuned,Gp)  ;

Fcl_fuzzy_autoTuned_PID = Fcl_PID_fuzzy                          ;
So_MO_PID_fuzzy_autotuned  = 1 - Fcl_fuzzy_autoTuned_PID         ; 
Si_MO_PID_fuzzy_autotuned  = series(So_MO_PID_fuzzy_autotuned,Gp);

% Calculating the "settling time" for each one of the tuning methods
% I, PI, PID Controller
% *************************************************************************
optimal_I_struct = stepinfo(Fcl_MO_I_optimal)                            ;
optimal_PI_struct = stepinfo(Fcl_MO_PI_optimal)                          ;
optimal_PID_struct = stepinfo(Fcl_MO_PID_optimal)                        ;

autotune_I_struct = stepinfo(Fcl_autoTuned_I)                            ;
autotune_PI_struct = stepinfo(Fcl_autoTuned_PI)                          ;
autotune_PID_struct = stepinfo(Fcl_autoTuned_PID)                        ;
autotune_PID_struct_marg = stepinfo(Fcl_PID_auto_marg);

fuzzy_autotune_I_struct = stepinfo(Fcl_fuzzy_autoTuned_I)                ;
fuzzy_autotune_PI_struct = stepinfo(Fcl_fuzzy_autoTuned_PI)              ;
fuzzy_autotune_PID_struct = stepinfo(Fcl_fuzzy_autoTuned_PI)             ;


settlingtime_I_optimal = optimal_I_struct.SettlingTime                   ;
settlingtime_PI_optimal = optimal_PI_struct.SettlingTime                 ; 
settlingtime_PID_optimal = optimal_PID_struct.SettlingTime               ;

settlingtime_I_autotuned = autotune_I_struct.SettlingTime                ;
settlingtime_PI_autotuned = autotune_PI_struct.SettlingTime              ;
settlingtime_PID_autotuned = autotune_PID_struct.SettlingTime            ;
settlingtime_PID_autotuned_marg = autotune_PID_struct_marg.SettlingTime  ;

settlingtime_I_fuzzy_autotuned = fuzzy_autotune_I_struct.SettlingTime    ;
settlingtime_PI_fuzzy_autotuned = fuzzy_autotune_PI_struct.SettlingTime  ;
settlingtime_PID_fuzzy_autotuned = fuzzy_autotune_PID_struct.SettlingTime;

close all

% % -------------------------------------------------------------------------
% %*************************    I Controller   ******************************
% %*************************       Optimal     ******************************
% % -------------------------------------------------------------------------
NoofSamples = 3000;


FinalTime_I_optimal = 1.2*settlingtime_I_optimal                ;
FinalTime_I_autotuned = 1.2*settlingtime_I_autotuned            ;
FinalTime_I_fuzzy_autotuned = 1.2*settlingtime_I_fuzzy_autotuned;

sortMatr_I = [FinalTime_I_optimal FinalTime_I_autotuned FinalTime_I_fuzzy_autotuned];
FinalTime_I = max(sortMatr_I);



% % -------------------------------------------------------------------------
% %*************************    PI Controller  ******************************
% %*************************       Optimal     ******************************
% % -------------------------------------------------------------------------
FinalTime_PI_optimal = 1.2*settlingtime_PI_optimal                ;
FinalTime_PI_autotuned = 1.2*settlingtime_PI_autotuned            ;
FinalTime_PI_fuzzy_autotuned = 1.2*settlingtime_PI_fuzzy_autotuned;
sortMatr_PI = [FinalTime_PI_optimal FinalTime_PI_autotuned FinalTime_PI_fuzzy_autotuned];
FinalTime_PI = max(sortMatr_PI);



% % -------------------------------------------------------------------------
% %*************************    PID Controller ******************************
% %*************************       Optimal     ******************************
% % -------------------------------------------------------------------------
FinalTime_PID_optimal = 1.2*settlingtime_PID_optimal                    ;
FinalTime_PID_autotuned = 1.2*settlingtime_PID_autotuned                ;
FinalTime_PID_fuzzy_autotuned_marg = 1.2*settlingtime_PID_autotuned_marg; 
FinalTime_PID_fuzzy_autotuned = 1.2*settlingtime_PID_fuzzy_autotuned    ;

sortMatr_PID = [FinalTime_PID_optimal FinalTime_PID_autotuned FinalTime_PID_fuzzy_autotuned FinalTime_PID_fuzzy_autotuned_marg];
FinalTime_PID = max(sortMatr_PID);



ltiview(Fcl_MO_I_optimal,'--b',Fcl_auto_I,'r',Fcl_I_fuzzy,'g')
ltiview(Fcl_MO_PI_optimal,'--b',Fcl_auto_PI,'r',Fcl_PI_fuzzy,'g')
ltiview(Fcl_MO_PID_optimal,'--b',Fcl_auto_PID,'r',Fcl_PID_fuzzy,'g',Fcl_PID_auto_marg,'k')

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% ************************  SAME TIME SCALE *******************************
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
matrFinalTime = [FinalTime_I FinalTime_PI FinalTime_PID];
FinalTime = max(matrFinalTime)                          ;  
t1 = 0:(3*FinalTime/NoofSamples):3*FinalTime            ;
t2 = 0:(3*FinalTime/NoofSamples):2*FinalTime            ;
t3 = 0:(3*FinalTime/NoofSamples):FinalTime              ;

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% ************************  PID Controller  *******************************
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% PID Optimal Controller
% -------------------------------------------------------------------------
y_optimal_Fcl_PID = step(Fcl_MO_PID_optimal,t1)         ;
y_optimal_So_PID = step(So_MO_PID_optimal,t2)           ;
y_optimal_Si_PID = step(Si_MO_PID_optimal,t3)           ;

y_optimal_PID = zeros(NoofSamples + 1,1);
for i = 1:1:NoofSamples + 1
    y_optimal_PID(i) = y_optimal_Fcl_PID(i);
    if i > NoofSamples/3
        y_optimal_PID(i) = y_optimal_PID(i) + y_optimal_So_PID(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_optimal_PID(i) = y_optimal_PID(i) + y_optimal_Si_PID(i - (2/3)*NoofSamples);
        end
    end
end

% PID Autotuned Controller
% -------------------------------------------------------------------------
y_autotuned_Fcl_PID = step(Fcl_autoTuned_PID,t1)             ;
y_autotuned_So_PID = step(So_MO_PID_autotuned,t2)            ;
y_autotuned_Si_PID = step(Si_MO_PID_autotuned,t3)            ;

y_autotuned_PID = zeros(NoofSamples + 1,1);

for i = 1:1:NoofSamples + 1
    y_autotuned_PID(i) = y_autotuned_Fcl_PID(i);
    if i > NoofSamples/3
        y_autotuned_PID(i) = y_autotuned_PID(i) + y_autotuned_So_PID(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_autotuned_PID(i) = y_autotuned_PID(i) + y_autotuned_Si_PID(i - (2/3)*NoofSamples);
        end
    end
end

% PID Auto Tuned Controller [Margaris Estimation]
% -------------------------------------------------------------------------
y_auto_autotuned_marg_Fcl_PID = step(Fcl_PID_auto_marg,t1) ;
y_auto_autotuned_marg_So_PID = step(So_MO_PID_auto_marg,t2);
y_auto_autotuned_marg_Si_PID = step(Si_MO_PID_auto_marg,t3);

y_auto_marg_PID = zeros(NoofSamples + 1,1);

for i = 1:1:NoofSamples + 1
    y_auto_marg_PID(i) = y_auto_autotuned_marg_Fcl_PID(i);
    if i > NoofSamples/3
        y_auto_marg_PID(i) = y_auto_autotuned_marg_Fcl_PID(i) + y_auto_autotuned_marg_So_PID(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
           y_auto_marg_PID(i) = y_auto_autotuned_marg_Fcl_PID(i) + y_auto_autotuned_marg_Si_PID(i - (2/3)*NoofSamples);
        end
    end
end

% PID Fuzzy Auto Tuned Controller
% -------------------------------------------------------------------------
if LagFlagPID >= 0.5

    y_fuzzy_autotuned_Fcl_PID = step(Fcl_fuzzy_autoTuned_PID,t1) ;
    y_fuzzy_autotuned_So_PID = step(So_MO_PID_fuzzy_autotuned,t2);
    y_fuzzy_autotuned_Si_PID = step(Si_MO_PID_fuzzy_autotuned,t3);

    y_fuzzy_autotuned_PID = zeros(NoofSamples + 1,1);

    for i = 1:1:NoofSamples + 1
        y_fuzzy_autotuned_PID(i) = y_fuzzy_autotuned_Fcl_PID(i);
        if i > NoofSamples/3
            y_fuzzy_autotuned_PID(i) = y_fuzzy_autotuned_PID(i) + y_fuzzy_autotuned_So_PID(i - NoofSamples/3);
            if i > (NoofSamples/3)*2
                y_fuzzy_autotuned_PID(i) = y_fuzzy_autotuned_PID(i) + y_fuzzy_autotuned_Si_PID(i - (2/3)*NoofSamples);
            end
        end
    end
end


% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% ************************  PI Controller  *******************************
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% PI Optimal Controller
% -------------------------------------------------------------------------
y_optimal_Fcl_PI = step(Fcl_MO_PI_optimal,t1)              ;
y_optimal_So_PI = step(So_MO_PI_optimal,t2)                ;
y_optimal_Si_PI = step(Si_MO_PI_optimal,t3)                ;

y_optimal_PI = zeros(NoofSamples + 1,1);
for i = 1:1:NoofSamples + 1
    y_optimal_PI(i) = y_optimal_Fcl_PI(i);
    if i > NoofSamples/3
        y_optimal_PI(i) = y_optimal_PI(i) + y_optimal_So_PI(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_optimal_PI(i) = y_optimal_PI(i) + y_optimal_Si_PI(i - (2/3)*NoofSamples);
        end
    end
end

% PI Autotuned Controller
% -------------------------------------------------------------------------
y_autotuned_Fcl_PI = step(Fcl_autoTuned_PI,t1)             ;
y_autotuned_So_PI = step(So_MO_PI_autotuned,t2)            ;
y_autotuned_Si_PI = step(Si_MO_PI_autotuned,t3)            ;

y_autotuned_PI = zeros(NoofSamples + 1,1);
for i = 1:1:NoofSamples + 1
    y_autotuned_PI(i) = y_autotuned_Fcl_PI(i);
    if i > NoofSamples/3
        y_autotuned_PI(i) = y_autotuned_PI(i) + y_autotuned_So_PI(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_autotuned_PI(i) = y_autotuned_PI(i) + y_autotuned_Si_PI(i - (2/3)*NoofSamples);
        end
    end
end


% PI Fuzzy Auto Tuned Controller
% -------------------------------------------------------------------------
if LagFlagPI >= 0.5
    y_fuzzy_autotuned_Fcl_PI = step(Fcl_fuzzy_autoTuned_PI,t1) ;
    y_fuzzy_autotuned_So_PI = step(So_MO_PI_fuzzy_autotuned,t2);
    y_fuzzy_autotuned_Si_PI = step(Si_MO_PI_fuzzy_autotuned,t3);

    y_fuzzy_autotuned_PI = zeros(NoofSamples + 1,1);
    for i = 1:1:NoofSamples + 1
        y_fuzzy_autotuned_PI(i) = y_fuzzy_autotuned_Fcl_PI(i);
        if i > NoofSamples/3
            y_fuzzy_autotuned_PI(i) = y_fuzzy_autotuned_PI(i) + y_fuzzy_autotuned_So_PI(i - NoofSamples/3);
            if i > (NoofSamples/3)*2
                y_fuzzy_autotuned_PI(i) = y_fuzzy_autotuned_PI(i) + y_fuzzy_autotuned_Si_PI(i - (2/3)*NoofSamples);
            end
        end
    end

end
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% ************************  I Controller  *******************************
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% I Optimal Controller
% -------------------------------------------------------------------------
y_optimal_Fcl_I = step(Fcl_MO_I_optimal,t1)              ;
y_optimal_So_I = step(So_MO_I_optimal,t2)                ;
y_optimal_Si_I = step(Si_MO_I_optimal,t3)                ;

y_optimal_I = zeros(NoofSamples + 1,1);
for i = 1:1:NoofSamples + 1
    y_optimal_I(i) = y_optimal_Fcl_I(i);
    if i > NoofSamples/3
        y_optimal_I(i) = y_optimal_I(i) + y_optimal_So_I(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_optimal_I(i) = y_optimal_I(i) + y_optimal_Si_I(i - (2/3)*NoofSamples);
        end
    end
end

% I Autotuned Controller
% -------------------------------------------------------------------------
y_autotuned_Fcl_I = step(Fcl_autoTuned_I,t1)             ;
y_autotuned_So_I = step(So_MO_I_autotuned,t2)            ;
y_autotuned_Si_I = step(Si_MO_I_autotuned,t3)            ;

y_autotuned_I = zeros(NoofSamples + 1,1);
for i = 1:1:NoofSamples + 1
    y_autotuned_I(i) = y_autotuned_Fcl_I(i);
    if i > NoofSamples/3
        y_autotuned_I(i) = y_autotuned_I(i) + y_autotuned_So_I(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_autotuned_I(i) = y_autotuned_I(i) + y_autotuned_Si_I(i - (2/3)*NoofSamples);
        end
    end
end

% I Fuzzy Auto Tuned Controller
% -------------------------------------------------------------------------
y_fuzzy_autotuned_Fcl_I = step(Fcl_fuzzy_autoTuned_I,t1) ;
y_fuzzy_autotuned_So_I = step(So_MO_I_fuzzy_autotuned,t2);
y_fuzzy_autotuned_Si_I = step(Si_MO_I_fuzzy_autotuned,t3);

y_fuzzy_autotuned_I = zeros(NoofSamples + 1,1);
for i = 1:1:NoofSamples + 1
    y_fuzzy_autotuned_I(i) = y_fuzzy_autotuned_Fcl_I(i);
    if i > NoofSamples/3
        y_fuzzy_autotuned_I(i) = y_fuzzy_autotuned_I(i) + y_fuzzy_autotuned_So_I(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_fuzzy_autotuned_I(i) = y_fuzzy_autotuned_I(i) + y_fuzzy_autotuned_Si_I(i - (2/3)*NoofSamples);
        end
    end
end
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% ************************  STEP RESPONSE PLOTS   *************************
    plotResults()
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% &&&&&&&&&&&&&&&&&&&&&&&&&&& BODE PLOTS &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Automatic Tuning
% -------------------------------------------------------------------------
figureIndex = figureIndex + 1                                                             ;
auto_tune_param_plot_ovs_fluctStruct.figureIndex = figureIndex                            ;
auto_tune_param_plot_ovs_fluctStruct.log_ovs_while_auto_tuning = log_ovs_while_auto_tuning;
auto_tune_param_plot_ovs_fluctStruct.step_Conv_autoTuning_realZeros_I_Controller = step_Conv_autoTuning_realZeros_I_Controller          ;
auto_tune_param_plot_ovs_fluctStruct.step_Conv_autoTuning_realZeros_PI_Controller = step_Conv_autoTuning_realZeros_PI_Controller        ;
auto_tune_param_plot_ovs_fluctStruct.step_Conv_autoTuning_realZeros_PID_Controller = step_Conv_autoTuning_realZeros_PID_Controller      ;
auto_tune_param_plot_ovs_fluctStruct.step_Conv_autoTuning_complexZeros_PID_Controller = step_Conv_autoTuning_complexZeros_PID_Controller;



auto_tune_param_plot_ovs_fluctStruct.ovs_I_auto_tune =   ovs_I_auto_tune              ;
auto_tune_param_plot_ovs_fluctStruct.ovs_PI_auto_tune =  ovs_PI_auto_tune             ;
auto_tune_param_plot_ovs_fluctStruct.ovs_PID_auto_tune = ovs_PID_auto_tune            ;
auto_tune_param_plot_ovs_fluctStruct.ovs_PID_auto_tune_ypred = ovs_PID_auto_tune_ypred;

% calling the ovs plot function
% -------------------------------------------------------------------------
auto_tune_param_plot_ovs_fluct(auto_tune_param_plot_ovs_fluctStruct)

            
% Fuzzy Automatic Tuning
% -------------------------------------------------------------------------
figureIndex = figureIndex + 1                                                         ;

auto_tune_param_fuzzy_plot_ovs_fluctStruct.figureIndex = figureIndex                            ;
auto_tune_param_fuzzy_plot_ovs_fluctStruct.log_ovs_while_fuzzy_pi_tuning = log_ovs_while_fuzzy_pi_tuning;
auto_tune_param_fuzzy_plot_ovs_fluctStruct.log_ovs_while_fuzzy_pid_tuning = log_ovs_while_fuzzy_pid_tuning;

auto_tune_param_fuzzy_plot_ovs_fluctStruct.step_Conv_fuzzyTuning_complexZeros_PI_Controller = step_Conv_fuzzyTuning_complexZeros_PI_Controller        ;
auto_tune_param_fuzzy_plot_ovs_fluctStruct.step_Conv_fuzzyTuning_complexZeros_PID_Controller = step_Conv_fuzzyTuning_complexZeros_PID_Controller      ;
auto_tune_param_fuzzy_plot_ovs_fluctStruct.ovs_PI_fuzzy_tune =  ovs_PI_fuzzy_tune                 ;
auto_tune_param_fuzzy_plot_ovs_fluctStruct.ovs_PID_fuzzy_tune = ovs_PID_fuzzy_tune                ;


% calling the ovs plot function
% -------------------------------------------------------------------------
auto_tune_param_fuzzy_plot_ovs_fluct(auto_tune_param_fuzzy_plot_ovs_fluctStruct)

% calling the ovs plot function for the whole tuning
% -------------------------------------------------------------------------
figureIndex = figureIndex + 1                                                         ;

auto_tune_param_all_plot_ovs_fluctStruct.ovs_Final_auto_tune = ovs_Final_auto_tune ;
auto_tune_param_all_plot_ovs_fluctStruct.ovs_Final_auto_tune_margaris = ovs_Final_auto_tune_margaris;
auto_tune_param_all_plot_ovs_fluctStruct.ovs_Final_fuzzy_tune = ovs_Final_fuzzy_tune;
auto_tune_param_all_plot_ovs_fluctStruct.figureIndex = figureIndex;
auto_tune_param_plot_ovs_fluct_all(auto_tune_param_all_plot_ovs_fluctStruct)

return

% Checking the ANFIS' for PI and PID using a large number of systems
% *************************************************************************

% Train the classifier for identifying whether or not an optimal PI can
% be developed
% *************************************************************************

% Classifier training
% *************************************************************************
load ov_pi_zeros_check.txt
load ov_pi_zeros.txt
load fismat_pi

ov_pi_control = ov_pi_zeros_check;

% Get the number of examples
% *************************************************************************
numExamples = size(ov_pi_zeros_check);

% Formulate the training data set
% *************************************************************************
index = 1;
k = 1    ;
Xtrain_pi_control = ov_pi_control;

for i = 1:numExamples
    if ((Xtrain_pi_control(i,3) < 3) || (Xtrain_pi_control(i,3) > 8))
        Xtrain_pi_control(i,3) = 0;
    end
end

for i = 1:numExamples
    if Xtrain_pi_control(i,3) ~= 0
       Xtrain_pi_control(i,3) = 1;
    end
end


max_pi(1) = max(ov_pi_zeros(:,1));
min_pi(1) = min(ov_pi_zeros(:,1));
max_pi(2) = max(ov_pi_zeros(:,2));
min_pi(2) = min(ov_pi_zeros(:,2));


for i = 1:numExamples
    for j = 1:2
        if (Xtrain_pi_control(i,j) < min_pi(j))
            fuzzyPI(j) = min_pi(j);
        elseif (Xtrain_pi_control(i,j) > max_pi(j))
            fuzzyPI(j) = max_pi(j);
        else
            fuzzyPI(j) = Xtrain_pi_control(i,j);
        end
    end

    %Calculate the ouput of ANFIS for the training data set
    % *********************************************************************
    trn_out_fismat_control_pi(i) = evalfis(fuzzyPI,fismat_control_pi);
    
end


% Calculate the ouput of ANFIS for the training data set
% trn_out_fismat_control_pi = evalfis(Xtrain_pi_control(:,1:2),fismat_control_pi);
% *************************************************************************
trn_out_fismat_control_pi = round(trn_out_fismat_control_pi);
correct = 0;

for i = 1:numExamples(1)
    if (trn_out_fismat_control_pi(i) == Xtrain_pi_control(i,3))
        correct = correct + 1;
    end
end


Percentage = correct/numExamples(1)*100;

% **************************  PREDICTION PROBLEM  *************************
% *************************************************************************
% Predict the desired PI overshoot when the plant
% can be controlled by a PI controller

%%%ONLY FOR THE WARNINGS TO GO OFF!!!
% *************************************************************************
ov_pi_control_warn = ov_pi_zeros;

% Get the number of examples
% *************************************************************************
numExamples_warn=size(ov_pi_zeros);
    
%Formulate the training data set
% *************************************************************************
index = 1;
k = 1    ;
Xtrain_pi_control_warn = ov_pi_control_warn;
    
    for i = 1:numExamples_warn
        if ((Xtrain_pi_control_warn(i,3) < 3) || (Xtrain_pi_control_warn(i,3) > 8))
            Xtrain_pi_control_warn(i,3) = 0 ;
        end
    end
    
    for i = 1:numExamples_warn
        if Xtrain_pi_control_warn(i,3) ~= 0
           Xtrain_pi_control_warn(i,3) = 1 ;
        end
    end

    ov_pi_pred_warn=ov_pi_zeros;
    
    
%Formulate the training data set
% *************************************************************************
index = 1;
k = 1    ;

for i = 1:numExamples_warn
   if (Xtrain_pi_control_warn(i,3) == 1)
       Xtrain_pred_pi_warn(index,:) = ov_pi_pred_warn(i,:);
       index = index + 1;
   end
end

max_pi(1) = max(Xtrain_pred_pi_warn(:,1));
min_pi(1) = min(Xtrain_pred_pi_warn(:,1));
max_pi(2) = max(Xtrain_pred_pi_warn(:,2));
min_pi(2) = min(Xtrain_pred_pi_warn(:,2));    
    
ov_pi_pred = ov_pi_zeros_check;


% Formulate the training data set
% *************************************************************************
index = 1;
k = 1    ;

for i = 1:numExamples
    if (Xtrain_pi_control(i,3) == 1)
        Xtrain_pred_pi(index,:) = ov_pi_pred(i,:);
        index = index + 1;
    end 
end
numExamples = size(Xtrain_pred_pi);


for i = 1:numExamples
    for j = 1:2
        if (Xtrain_pred_pi(i,j) < min_pi(j))
            fuzzyPI(j) = min_pi(j);
        elseif (Xtrain_pred_pi(i,j) > max_pi(j))
            fuzzyPI(j) = max_pi(j);
        else
            fuzzyPI(j) = Xtrain_pred_pi(i,j);
        end
    end

    % Calculate the ouput of ANFIS for the training data set
    % *********************************************************************
    trn_out_fismat_pred_pi(i) = evalfis(fuzzyPI,fismat_pred_pi);
end

% Calculate the ouput of ANFIS for the training data set
% *************************************************************************
% trn_out_fismat_pred_pi=evalfis(Xtrain_pred_pi(:,1:2),fismat_pred_pi);
% *************************************************************************

e = trn_out_fismat_pred_pi' - Xtrain_pred_pi(:,3);
error = sqrt(mse(e))

figure(1)
hist(abs(e),100)
title('Checking errors distribution')
xlabel('Error in the Prediction of the Desired Oveshooot')

% *************************************************************************
load ov_pid_zeros_check.txt
load ov_pid_zeros.txt
load fismat_pid

% Load the final with the training inputs
% (ti_I ov_PI ti_PI Settling Time op-loop) - desiredPID ovs
% ---------------------------------------------------------------------

% *************************   CLASSIFICATION PROBLEM   ********************
% *************************************************************************
% Recognize whether or not the system can be controlled by a PI controller
% *************************************************************************
ov_pid_control = ov_pid_zeros_check;

% Get the number of examples
% *************************************************************************
numExamples = size(ov_pid_zeros_check);

% Formulate the training data set
% *************************************************************************
index = 1;
k = 1    ;

Xtrain_pid_control = ov_pid_control;

for i = 1:numExamples
    if ((Xtrain_pid_control(i,5) < 3)||( Xtrain_pid_control(i,5) > 8.5 ))
        Xtrain_pid_control(i,5) = 0 ;
    else
        Xtrain_pid_control(i,5) = 1 ;
    end

end


[min_pid(1),min_ind(1)] = min(ov_pid_zeros(:,1));
[min_pid(2),min_ind(2)] = min(ov_pid_zeros(:,2));
[min_pid(3),min_ind(3)] = min(ov_pid_zeros(:,3));
[min_pid(4),min_ind(4)] = min(ov_pid_zeros(:,4));
[max_pid(1),max_ind(1)] = max(ov_pid_zeros(:,1));
[max_pid(2),max_ind(2)] = max(ov_pid_zeros(:,2));
[max_pid(3),max_ind(3)] = max(ov_pid_zeros(:,3));
[max_pid(4),max_ind(4)] = max(ov_pid_zeros(:,4));


for i = 1:numExamples
    for j = 1:4
        if (Xtrain_pid_control(i,j) < min_pid(j))
            fuzzyPID(j) = min_pid(j);
        elseif (Xtrain_pid_control(i,j) > max_pid(j))
            fuzzyPID(j) = max_pid(j);
        else
            fuzzyPID(j) = Xtrain_pid_control(i,j);
        end
    end
    % Calculate the ouput of ANFIS for the training data set
    % *********************************************************************
    trn_out_fismat_control_pid(i) = evalfis(fuzzyPID,fismat_control_pid);
    
end

% %Calculate the ouput of ANFIS for the training data set
% *************************************************************************
% trn_out_fismat_control_pid=evalfis(Xtrain_pid_control(:,1:4),fismat_control_pid);
% *************************************************************************
trn_out_fismat_control_pid = round(trn_out_fismat_control_pid);
correct = 0;

for i = 1:numExamples(1)
    if (trn_out_fismat_control_pid(i) == Xtrain_pid_control(i,5))
        correct = correct + 1;
    end
end

Percentage = correct/numExamples(1)*100  ;

% ****************************  PREDICTION PROBLEM ************************
% *************************************************************************
% Predict the desired PI overshoot when the plant can be controlled by a PI
% controller
% *************************************************************************

% ONLY FOR THE WARNINGS TO GO OFF!!!
% *************************************************************************
ov_pid_control_warn = ov_pid_zeros;

% Get the number of examples
% *************************************************************************
numExamples_warn=size(ov_pid_zeros);
    
% Formulate the training data set
% *************************************************************************
index = 1 ;
k = 1     ;
Xtrain_pid_control_warn = ov_pid_control_warn;

for i=1:numExamples_warn
   if ((Xtrain_pid_control_warn(i,5) < 3) || (Xtrain_pid_control_warn(i,5) > 8.5))
      Xtrain_pid_control_warn(i,5) = 0;
   end
end


for i = 1:numExamples_warn
   if Xtrain_pid_control_warn(i,5) ~= 0
      Xtrain_pid_control_warn(i,5) = 1 ;
   end
end

ov_pid_pred_warn = ov_pid_zeros;
    
    
% Formulate the training data set
% *************************************************************************
index = 1;
k = 1    ;

for i = 1:numExamples_warn
   if (Xtrain_pid_control_warn(i,5) == 1)
       Xtrain_pred_pid_warn(index,:) = ov_pid_pred_warn(i,:);
       index = index + 1;
   end 
end

[min_pid(1),min_ind(1)] = min(Xtrain_pred_pid_warn(:,1));
[min_pid(2),min_ind(2)] = min(Xtrain_pred_pid_warn(:,2));
[min_pid(3),min_ind(3)] = min(Xtrain_pred_pid_warn(:,3));
[min_pid(4),min_ind(4)] = min(Xtrain_pred_pid_warn(:,4));
[max_pid(1),max_ind(1)] = max(Xtrain_pred_pid_warn(:,1));
[max_pid(2),max_ind(2)] = max(Xtrain_pred_pid_warn(:,2));
[max_pid(3),max_ind(3)] = max(Xtrain_pred_pid_warn(:,3));
[max_pid(4),max_ind(4)] = max(Xtrain_pred_pid_warn(:,4));

ov_pid_pred = ov_pid_zeros_check;
% *************************************************************************

%Get the number of examples
% *************************************************************************
numExamples = size(ov_pid_zeros_check);

% Formulate the training data set
% *************************************************************************
index = 1 ;
k = 1     ; 

for i = 1:numExamples
    if (Xtrain_pid_control(i,5) == 1)
        Xtrain_pred_pid(index,:) = ov_pid_pred(i,:);
        index = index + 1;
    end 
end

numExamples = size(Xtrain_pred_pid);

for i = 1:numExamples
   for j = 1:4
        if (Xtrain_pred_pid(i,j) < min_pid(j))
            fuzzyPID(j) = min_pid(j);
        elseif (Xtrain_pred_pid(i,j) > max_pid(j))
            fuzzyPID(j) = max_pid(j);
        else
            fuzzyPID(j) = Xtrain_pred_pid(i,j);
        end
   end

    %Calculate the ouput of ANFIS for the training data set
    % *********************************************************************
    trn_out_fismat_pred_pid(i) = evalfis(fuzzyPID,fismat_pred_pid);
    
end

% %Calculate the ouput of ANFIS for the training data set
% *************************************************************************
% trn_out_fismat_pred_pid=evalfis(Xtrain_pred_pid(:,1:4),fismat_pred_pid);
% *************************************************************************
e_pid = trn_out_fismat_pred_pid' - Xtrain_pred_pid(:,5);

error_pid = sqrt(mse(e_pid));

% *************************************************************************
figure(2)
hist(abs(e_pid),100)
title('Checking errors distribution')
xlabel('Error in the Prediction of the Desired Oveshooot')

% -------------------------------------------------------------------------
%EOF: auto_tune_param_main.m
