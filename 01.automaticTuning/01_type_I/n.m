% *********************************************************************************************************** 
%               Copyright (C) 2008
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_2ndOrder_main.m																			   																		  	
%  Project:     Automatic tuning of the parameters for PI,PID controllers
%  
%  Purpose:     main script for the automatic tuning based on the three methods																		   																		
%  Author :     kostas g. papadopoulos																	   																		
% 																										   																		
%  History:     Date: 07.07.2008  date last modified
% 																										  																		
%  Contact:     leonidas droukas   ,       kostas g. papadopoulos    
%               leon_drouk@yahoo.gr,       kpapadop@eng.auth.gr      
% 																										  																		
%  Place:	    Aristotle University of Thessaloniki, Thessaloniki, Greece							   																		
% 
% ************************************************************************************************************
figureIndex = 0;
disp('********************************************************************')
disp('enter "y" for random plant or "n" for manual plant')
disp('********************************************************************')
randPlant = input('y/n:','s');
switch randPlant
  case {'n'}
      disp('----------------------------------------------------------------------');
      disp('Enter T')
      disp('----------------------------------------------------------------------');
      T = input('T = ');    
      disp('----------------------------------------------------------------------');
      fprintf('Enter unmodelled dynamic Ts [Ts is the sum unmodelled dynamics of \n')
      fprintf('the controller and the plant]\n')
      disp('----------------------------------------------------------------------');
      Ts = input('Ts = ');
      disp('----------------------------------------------------------------------');
      disp('Enter kp')
      disp('----------------------------------------------------------------------');
      kp = input('kp = ');
      disp('----------------------------------------------------------------------');
      disp('Enter "zeta"')
      disp('----------------------------------------------------------------------');
      zeta = input('zeta = ');
      while ((zeta <= 0) || (zeta >= 1))
          disp('----------------------------------------------------------------------');
          disp('re-enter "zeta" parameter: [0 < zeta < 1]')
          disp('----------------------------------------------------------------------');
          zeta = input('zeta = ');
      end
  case {'y'}
      % generate random plant parameters
      % --------------------------------------------------------------------------------------------
      kp = rand    ;
      T = rand     ;
      Ts = 0.1*rand;
      zeta = rand;
      while zeta > 0.9
            zeta = rand ;
      end
          
      while ((zeta <= 0) || (zeta >= 1))
          zeta = rand;
      end
end

% -------------------------------------------------------------------------------
% kp = 0.6791;
% T  = 0.3955;
% Ts = 0.0367;
% zeta = 0.0377;
% -------------------------------------------------------------------------------

kh = 1       ;  
plant.kh = kh;

kp_init = kp ;   T_init = T      ;
Ts_init = Ts ;   zeta_init = zeta; 

plant.kp_arx = kp_init;    plant.T_arx = T_init      ;
plant.Ts_arx = Ts_init;    plant.zeta_arx = zeta_init;   

% Based on the open loop experiment of the plant we estimate the following charecteristics
% --------------------------------------------------------------------------------------------------
% kp  : plant gain
% zeta: zeta parameter of the second order term of the plant
% T   : plant time constant
% --------------------------------------------------------------------------------------------------
[kp_est zeta_est T_est kp_est_good zeta_est_good T_est_good] = auto_tune_2ndOrder_kp_zeta_T_estimation(plant);

plant.kp_est = kp_est;               plant.T_est = T_est                ;        
                                     plant.zeta_est = zeta_est          ; 
plant.kp_est_good = kp_est_good;     plant.T_est_good = T_est_good      ;        
                                     plant.zeta_est_good = zeta_est_good;

fprintf('zeta-init: %5.4f, T-init: %5.4f \n',zeta_init,T_init)                                     
fprintf('zeta-estg: %5.4f, T-estg: %5.4f \n',zeta_est_good,T_est_good)                                     
fprintf('zeta-est : %5.4f, T-est : %5.4f \n',zeta_est,T_est)

% =================================================================
% kp_est_good  : the error is the result of equation approximation;
% T_est_good   : the error is the result of equation approximation;
% zeta_est_good: the error is the result of equation approximation; 
% =================================================================
% Forcing linear error (artificial error)
% =================================================================
% kp_est  : (1 + flagError_kp)*plant_loc.kp_init    ;
% T_est   : (1 + flagError_zeta)*plant_loc.zeta_init;
% zeta_est: (1 + flagError_T)*T_init                ; 
% =================================================================
                           
% check plant open loop response
% --------------------------------------------------------------------------------------------------
checkOpenLoopPlant = input('check the plant open loop step response - [y/n]:','s');
figureIndex = figureIndex + 1;
if strcmp(checkOpenLoopPlant ,'y') == 1
    num_Gp_est = kp_est                      ;  
    den_Gp_est = [T_est^2 2*zeta_est*T_est 1];
    Gp_est = tf(num_Gp_est,den_Gp_est)       ;

    num_Gp_est_good = kp_est_good                      ;  
    den_Gp_est_good = [T_est_good^2 2*zeta_est_good*T_est_good 1];
    Gp_est_good = tf(num_Gp_est_good,den_Gp_est_good)       ;

    figure(figureIndex)
    tspan = 0:0.01:50;
    y_Gp_est_good = step(Gp_est_good,tspan);
    y_Gp_est = step(Gp_est,tspan)          ;
    plot(tspan,y_Gp_est_good,'k',tspan,y_Gp_est,'r','LineWidth',3)
    grid
    figureIndex = figureIndex + 1          ;

    % Frequency Response
    % ================================
    figure(figureIndex)
    w = logspace(-2,2,6000);
    bodemag(Gp_est_good,Gp_est,w)     ;
    grid on
    hold on
    
    % Get the handle to the line object
    % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    HFig = gcf;
    hC = get(HFig,'children')           ;
    hg = get(hC(3),'children')          ;     % Grab a handle to the axes
                                         % Check the Type property to make sure
    hgC_Fcl_good = get(hg(1),'children');     % Grab a handle to the hg group
    hgC_Fcl  = get(hg(2),'children');

    % set the property you wish to modify
    % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    set(hgC_Fcl,'Color','k')              % changes the color
    set(hgC_Fcl,'LineStyle','-')          % changes the linestyle
    set(hgC_Fcl,'LineWidth',3)            % changes the Linewidth
        
    set(hgC_Fcl_good,'Color',[0.39 0.47 0.59])  % changes the color
    set(hgC_Fcl_good,'LineStyle','-')           % changes the linestyle
    set(hgC_Fcl_good,'LineWidth',3)             % changes the Linewidth
end

disp('----------------------------------------------------------------------');
plotNo = input('plottimes = ');
upperlimitOvrst = 4.33;
lowerlimitOvrst = 4.31;
fprintf('Upper Limit Overshoot: %5.2f \n',upperlimitOvrst)
fprintf('Lower Limit Overshoot: %5.2f \n',lowerlimitOvrst)
disp('----------------------------------------------------------------------');


% Tx: initialization
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Tx_init = T_est;
Tx = Tx_init   ;

% Tuning Ti_1 (Tsx) until overshoot of the closed loop step response reaches 4.32% ovs
% -------------------------------------------------------------------------------------------------
stepNo_Tsx = 0  ;

% Initializing the Controller Unmodelled dynamic with Ts_init = Tsc + Tsp
% -------------------------------------------------------------------------------------------------
Tsx = Ts_init      ;
Tsx_good = Ts_init ;
plant.Tsx = Tsx    ;
plant.Tsx_good = Tsx_good  ;

% checking if I control is ok (if not switch to I-LAG)
% -------------------------------------------------------------------------------------------------
[ovrst_loc_1 Fcl_loc_1 Ti_1] = auto_tune_2ndOrder_calculate_overshoot_1(plant)                       ;
[ovrst_loc_1_good Fcl_loc_1_good Ti_1_good] = auto_tune_2ndOrder_calculate_overshoot_1_goodEst(plant);

% i proti timi pou dokimazetai gia to Tsx einai isi me Ts. Afto iparxei
% periptosi na odigisei se astathes sistima kleistou vrogxou me apotelesma 
% na min mporei na sinexisei i diadikasia aftomatis rithmisis tou elegkti. Gia afto to
% logo elegxoume an ginetai kati tetoio kai se periptosi pou ginetai,
% allazoume tin proti timi tou Tsx.
% --------------------------------------------------------------------------------------------------
counter = 0     ;
counter_good = 0;
relaxCounter = 1;
while ( isnan(ovrst_loc_1) == 1)
        Tsx = Tsx + Tx       ;
        plant.Tsx = Tsx      ;
        counter = counter + 1;
        [ovrst_loc_1 Fcl_loc_1 Ti_1] = auto_tune_2ndOrder_calculate_overshoot_1(plant);
end
while ( isnan(ovrst_loc_1_good) == 1)
        Tsx_good = Tsx_good + Tx   ;
        plant.Tsx_good = Tsx_good  ;
        counter_good = counter_good + 1;
        [ovrst_loc_1_good Fcl_loc_1_good Ti_1_good] = auto_tune_2ndOrder_calculate_overshoot_1_goodEst(plant);
end
% ---------------------------------------------------------------------------------------------------------
% Automatic tuning with the wrong estimation of kp, zeta, T parameters after the open loop plant experiment
% ---------------------------------------------------------------------------------------------------------
fprintf('Starting the automatic tuning with the bad estimation of the plant paremeters...\n')
disp('----------------------------------------------------------------------')
% ---------------------------------------------------------------------------------------------------------
figureIndex = figureIndex + 1;
error_1 = upperlimitOvrst - ovrst_loc_1;
while ((ovrst_loc_1 < lowerlimitOvrst) || (ovrst_loc_1 > upperlimitOvrst))
    stepNo_Tsx = stepNo_Tsx + 1;
    if (error_1 < 0) 
        Tsx = Tsx + (0.025*Tsx);                     
        plant.Tsx = Tsx;
        [ovrst_loc_1 Fcl_loc_1 Ti_1] = auto_tune_2ndOrder_calculate_overshoot_1(plant);
        k = mod(stepNo_Tsx,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_loc_1,'k')
           title('T_{i1} - tuning - I Controller')      
           hold on
        end
        error_1 = upperlimitOvrst - ovrst_loc_1;  
    elseif (error_1 > 0)               
        Tsx = Tsx - (0.025*Tsx);                      
        plant.Tsx = Tsx;                                
        [ovrst_loc_1 Fcl_loc_1 Ti_1] = auto_tune_2ndOrder_calculate_overshoot_1(plant);
        k = mod(stepNo_Tsx,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_loc_1,'k')
           title('T_{i1} - tuning - I Controller')      
           hold on
        end
        error_1 = upperlimitOvrst - ovrst_loc_1;
    end 
    Fcltr{stepNo_Tsx} = Fcl_loc_1;
%     if stepNo_Tsx == 100
%         return
%     end
    n = mod(stepNo_Tsx,1);
    if (n == 0)
      fprintf('step_Tsx: %d - Ti_1: %2.5f - overshoot: %2.5f \n',stepNo_Tsx,Ti_1,ovrst_loc_1)
    end
    % relaxing the band
    % ----------------------------------------------------------------------------------------------
    if mod(stepNo_Tsx,relaxCounter*100) == 0
        relaxCounter = relaxCounter + 1;

    % widen the band
    % ----------------------------------------------------------------------------------------------
       upperlimitOvrst = upperlimitOvrst + 0.01*upperlimitOvrst; 
       lowerlimitOvrst = lowerlimitOvrst - 0.01*lowerlimitOvrst;
       disp('----------------------------------------------------------------------');
       fprintf('Relaxing the Reference Band so that the tuning converges faster...\n')
       fprintf('New Upper Limit: %2.5f - New Lower Limit: %2.5f \n',upperlimitOvrst,lowerlimitOvrst)
       disp('----------------------------------------------------------------------');
    end
end
overshootTunedStruct = stepinfo(Fcl_loc_1)     ;
overshootTuned = overshootTunedStruct.Overshoot;
Ti_1_rythm = Ti_1                              ;
hold off
% ---------------------------------------------------------------------------------------------------------
disp('----------------------------------------------------------------------')
fprintf('Starting the automatic tuning with the good estimation of the plant paremeters...\n')
disp('----------------------------------------------------------------------')

% ---------------------------------------------------------------------------------------------------------
% Automatic tuning with the good estimation of kp, zeta, T parameters after the open loop plant experiment
% ---------------------------------------------------------------------------------------------------------
error_1 = upperlimitOvrst - ovrst_loc_1;
figureIndex = figureIndex + 1;
while ((ovrst_loc_1_good < lowerlimitOvrst) || (ovrst_loc_1_good > upperlimitOvrst))
    stepNo_Tsx = stepNo_Tsx + 1;
    if (error_1 < 0) 
        Tsx = Tsx + (0.025*Tsx);                     
        plant.Tsx = Tsx;
         [ovrst_loc_1_good Fcl_loc_1_good Ti_1_good] = auto_tune_2ndOrder_calculate_overshoot_1_goodEst(plant);
        k = mod(stepNo_Tsx,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_loc_1_good,'k')
           title('T_{i1_good} - tuning - I Controller')      
           hold on
        end
        error_1 = upperlimitOvrst - ovrst_loc_1_good;  
    elseif (error_1 > 0)               
        Tsx = Tsx - (0.025*Tsx);                      
        plant.Tsx = Tsx;                                
         [ovrst_loc_1_good Fcl_loc_1_good Ti_1_good] = auto_tune_2ndOrder_calculate_overshoot_1_goodEst(plant);
        k = mod(stepNo_Tsx,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_loc_1_good,'k')
           title('T_{i1_good} - tuning - I Controller')      
           hold on
        end
        error_1 = upperlimitOvrst - ovrst_loc_1_good;
    end 
    Fcltr{stepNo_Tsx} = Fcl_loc_1_good;
%     if stepNo_Tsx == 100
%         return
%     end
    n = mod(stepNo_Tsx,1);
    if (n == 0)
      fprintf('step_Tsx: %d - Ti_1: %2.5f - overshoot: %2.5f \n',stepNo_Tsx,Ti_1_good,ovrst_loc_1_good)
    end
    % relaxing the band
    % ----------------------------------------------------------------------------------------------
    if mod(stepNo_Tsx,relaxCounter*100) == 0
        relaxCounter = relaxCounter + 1;

    % widen the band
    % ----------------------------------------------------------------------------------------------
       upperlimitOvrst = upperlimitOvrst + 0.01*upperlimitOvrst; 
       lowerlimitOvrst = lowerlimitOvrst - 0.01*lowerlimitOvrst;
       disp('----------------------------------------------------------------------');
       fprintf('Relaxing the Reference Band so that the tuning converges faster...\n')
       fprintf('New Upper Limit: %2.5f - New Lower Limit: %2.5f \n',upperlimitOvrst,lowerlimitOvrst)
       disp('----------------------------------------------------------------------');
    end
end
overshootTunedStruct = stepinfo(Fcl_loc_1_good)     ;
overshootTuned_good = overshootTunedStruct.Overshoot;
Ti_1_rythm_good = Ti_1_good                         ;

% Print out the results of the tuning
% ----------------------------------------------------------------------------------------------
disp('----------------------------------------------------------------------');
fprintf('Overshoot Fcl(s) after Tuning: %2.5f - Ti_1: %2.5f\n',overshootTuned,Ti_1_rythm)       
disp('----------------------------------------------------------------------');
fprintf('Overshoot Good Fcl(s) after Tuning: %2.5f - Ti_1_good: %2.5f\n',overshootTuned_good,Ti_1_rythm_good)       
disp('----------------------------------------------------------------------');
ltiview(Fcl_loc_1_good,'k',Fcl_loc_1,'r')  


% calculate Final Tx
% --------------------------------------------------------------------------------------------------
% Tsx : is the final value in the while Loop
% me vasi tis proseggiseis ton kp , zeta kai T pou eginan pio pano kathos kai
% tin rithmismeni timi tis statheras Ti_1 , ipologizetai mia
% proseggistiki timi gia to Ts tou sistimatos
% --------------------------------------------------------------------------------------------------
Ts_est = (Ti_1_rythm/(2*kp_est*kh)) - counter*Tx - 2*zeta_est*T_est;
Ts_est_good = (Ti_1_rythm_good/(2*kp_est_good*kh)) - counter_good*Tx - 2*zeta_est_good*T_est_good;

plant.Ts_est = Ts_est;   
plant.Ts_est_good = Ts_est_good; 

% Ti_2 tuning
% --------------------------------------------------------------------------------------------------
stepNo_Ti_2 = 0          ;
Ti_2 = 2*kp_est*kh*Ts_est;
plant.Ti_2 = Ti_2        ; 
[ovrst_loc_2 Fcl_loc_2 Gp_loc] = auto_tune_2ndOrder_calculate_overshoot_2(plant);

error_2 = upperlimitOvrst - ovrst_loc_2 ;
figureIndex = figureIndex + 1;
relaxCounter = 1             ;
while ((ovrst_loc_2 < lowerlimitOvrst) || (ovrst_loc_2 > upperlimitOvrst))
    stepNo_Ti_2= stepNo_Ti_2 + 1;
    if (error_2 < 0) 
        Ti_2 = Ti_2 + (0.025*Ti_2);                     
        plant.Ti_2 = Ti_2;
        [ovrst_loc_2 Fcl_loc_2 Gp_loc] = auto_tune_2ndOrder_calculate_overshoot_2(plant);
        k = mod(stepNo_Ti_2,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_loc_2,'k')
           title('T_{i2} - tuning - PID Controller')      
           hold on
        end
        error_2 = upperlimitOvrst - ovrst_loc_2;
    elseif (error_2 > 0)  
        Ti_2 = Ti_2 - (0.025*Ti_2);                      
        plant.Ti_2 = Ti_2;                                
        [ovrst_loc_2 Fcl_loc_2 Gp_loc] = auto_tune_2ndOrder_calculate_overshoot_2(plant);
        k = mod(stepNo_Ti_2,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_loc_2,'k')
           title('T_{i2} - tuning - PID Controller')      
           hold on
        end
        error_2 = upperlimitOvrst - ovrst_loc_2;
    end 
    n = mod(stepNo_Ti_2,1);
    if (n == 0)
      fprintf('step_Tsx: %d - Ti_2: %2.5f - overshoot: %2.5f \n',stepNo_Ti_2,Ti_2,ovrst_loc_2)
    end

    % relaxing the band
    % ----------------------------------------------------------------------------------------------
    if mod(stepNo_Ti_2,relaxCounter*100) == 0
        relaxCounter = relaxCounter + 1;
    % widen the band
    % ----------------------------------------------------------------------------------------------
       upperlimitOvrst = upperlimitOvrst + 0.01*upperlimitOvrst; 
       lowerlimitOvrst = lowerlimitOvrst - 0.01*lowerlimitOvrst;
       disp('----------------------------------------------------------------------');
       fprintf('Relaxing the Reference Band so that the tuning converges faster...\n')
       fprintf('New Upper Limit: %2.5f - New Lower Limit: %2.5f \n',upperlimitOvrst,lowerlimitOvrst)
       disp('----------------------------------------------------------------------');
    end
end
overshootTunedStruct2 = stepinfo(Fcl_loc_2)      ;
overshootTuned2 = overshootTunedStruct2.Overshoot;


% Ti_2 tuning [good estimation]
% --------------------------------------------------------------------------------------------------
stepNo_Ti_2 = 0                         ;
Ti_2_good = 2*kp_est_good*kh*Ts_est_good;
plant.Ti_2_good = Ti_2_good             ; 
[ovrst_loc_2_good Fcl_loc_2_good Gp_loc_good] = auto_tune_2ndOrder_calculate_overshoot_2_good(plant);

error_2 = upperlimitOvrst - ovrst_loc_2 ;
figureIndex = figureIndex + 1;
relaxCounter = 1             ;
while ((ovrst_loc_2 < lowerlimitOvrst) || (ovrst_loc_2 > upperlimitOvrst))
    stepNo_Ti_2= stepNo_Ti_2 + 1;
    if (error_2 < 0) 
        Ti_2_good = Ti_2_good + (0.025*Ti_2_good);                     
        plant.Ti_2_good = Ti_2_good;
        [ovrst_loc_2_good Fcl_loc_2_good Gp_loc_good] = auto_tune_2ndOrder_calculate_overshoot_2_good(plant);
        k = mod(stepNo_Ti_2,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_loc_2_good,'k')
           title('T_{i2} - tuning - PID Controller')      
           hold on
        end
        error_2 = upperlimitOvrst - ovrst_loc_2;
    elseif (error_2 > 0)  
        Ti_2_good = Ti_2_good - (0.025*Ti_2_good);                      
        plant.Ti_2 = Ti_2;                                
        [ovrst_loc_2_good Fcl_loc_2_good Gp_loc_good] = auto_tune_2ndOrder_calculate_overshoot_2_good(plant);
        k = mod(stepNo_Ti_2,plotNo);
        if (k == 0)
           figure(figureIndex)
           step(Fcl_loc_2_good,'k')
           title('T_{i2} - tuning - PID Controller')      
           hold on
        end
        error_2 = upperlimitOvrst - ovrst_loc_2;
    end 
    n = mod(stepNo_Ti_2,1);
    if (n == 0)
      fprintf('step_Tsx: %d - Ti_2: %2.5f - overshoot: %2.5f \n',stepNo_Ti_2,Ti_2_good,ovrst_loc_2_good)
    end

    % relaxing the band
    % ----------------------------------------------------------------------------------------------
    if mod(stepNo_Ti_2,relaxCounter*100) == 0
        relaxCounter = relaxCounter + 1;
    % widen the band
    % ----------------------------------------------------------------------------------------------
       upperlimitOvrst = upperlimitOvrst + 0.01*upperlimitOvrst; 
       lowerlimitOvrst = lowerlimitOvrst - 0.01*lowerlimitOvrst;
       disp('----------------------------------------------------------------------');
       fprintf('Relaxing the Reference Band so that the tuning converges faster...\n')
       fprintf('New Upper Limit: %2.5f - New Lower Limit: %2.5f \n',upperlimitOvrst,lowerlimitOvrst)
       disp('----------------------------------------------------------------------');
    end
end
overshootTunedStruct2_good = stepinfo(Fcl_loc_2_good)      ;
overshootTuned2_good = overshootTunedStruct2_good.Overshoot;


ltiview(Fcl_loc_1_good,'k',Fcl_loc_2_good,'--k',Fcl_loc_1,'r',Fcl_loc_2,'--r')
ltiview(Fcl_loc_1_good,'k',Fcl_loc_2_good,'--k')

% Print out the results of the tuning
% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
disp('----------------------------------------------------------------------');
fprintf('Overshoot Fcl(s) after Tuning: %2.5f - Ti_2: %2.5f\n',overshootTuned2,Ti_2)       
disp('----------------------------------------------------------------------');
Gp = Gp_loc;
So_MO_loc_1  = 1 - Fcl_loc_1         ; 
Si_MO_loc_1  = series(So_MO_loc_1,Gp);

So_MO_loc_2  = 1 - Fcl_loc_2         ; 
Si_MO_loc_2  = series(So_MO_loc_2,Gp);

So_MO_loc_1_good  = 1 - Fcl_loc_1_good         ; 
Si_MO_loc_1_good  = series(So_MO_loc_1_good,Gp);

So_MO_loc_2_good  = 1 - Fcl_loc_2_good         ; 
Si_MO_loc_2_good  = series(So_MO_loc_2_good,Gp);

tspan = 0:0.001:100                ;
y_Fcl_loc_1 = step(Fcl_loc_1,tspan);
y_Fcl_loc_2 = step(Fcl_loc_2,tspan);
y_Fcl_loc_1_good = step(Fcl_loc_1_good,tspan);
y_Fcl_loc_2_good = step(Fcl_loc_2_good,tspan);

return
figureIndex = figureIndex + 1;
figure(figureIndex)
plot(tspan,y_Fcl_loc_1_good,'k',tspan,y_Fcl_loc_2_good,'--k',...
     tspan,y_Fcl_loc_1,'r',tspan,y_Fcl_loc_2,'--r','LineWidth',3)
grid on

figureIndex = figureIndex + 1;
figure(figureIndex)
plot(tspan,y_Fcl_loc_1_good,'k',tspan,y_Fcl_loc_1,'r','LineWidth',3)
grid on

figureIndex = figureIndex + 1;
figure(figureIndex)
plot(tspan,y_Fcl_loc_2_good,'--k',...
     tspan,y_Fcl_loc_2,'--r','LineWidth',3)
grid on
return
% Plot the results
% ----------------------------------------------------------------------------------------------
Fcl_1_struct = stepinfo(Fcl_loc_1)                             ;
Fcl_2_struct = stepinfo(Fcl_loc_2)                             ;

Fcl_1_good_struct = stepinfo(Fcl_loc_1_good)                   ;
Fcl_2_good_struct = stepinfo(Fcl_loc_2_good)                   ;

settlingtime_Fcl_1_struct = Fcl_1_struct.SettlingTime          ;
settlingtime_Fcl_2_struct = Fcl_2_struct.SettlingTime          ;

settlingtime_Fcl_1_good_struct = Fcl_1_good_struct.SettlingTime;
settlingtime_Fcl_2_good_struct = Fcl_2_good_struct.SettlingTime;


FinalTime_Fcl_1 = 2*settlingtime_Fcl_1_struct      ;
FinalTime_Fcl_2 = 2*settlingtime_Fcl_2_struct      ;

FinalTime_Fcl_1_good = 2*settlingtime_Fcl_1_good_struct      ;
FinalTime_Fcl_2_good = 2*settlingtime_Fcl_2_good_struct      ;

sortMatr_Fcl = [FinalTime_Fcl_1 FinalTime_Fcl_2...
                FinalTime_Fcl_1_good FinalTime_Fcl_2_good]   ;
FinalTime = max(sortMatr_Fcl)                      ;

NoofSamples = 3000                                 ;
t1 = 0:(3*FinalTime/NoofSamples):3*FinalTime       ;
t2 = 0:(3*FinalTime/NoofSamples):2*FinalTime       ;
t3 = 0:(3*FinalTime/NoofSamples):FinalTime         ;
% ----------------------------------------------------------------------------------------------
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% ----------------------------------------------------------------------------------------------
y_ti_1_Fcl_good = step(Fcl_loc_1_good,t1)   ;
y_ti_1_So_good = step(So_MO_loc_1_good,t2)  ;
y_ti_1_Si_good = step(Si_MO_loc_1_good,t3)  ;

y_ti_1_good = zeros(NoofSamples + 1,1);
for i = 1:1:NoofSamples + 1
    y_ti_1_good(i) = y_ti_1_Fcl_good(i);
    if i > NoofSamples/3
        y_ti_1_good(i) = y_ti_1_good(i) + y_ti_1_So_good(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_ti_1_good(i) = y_ti_1_good(i) + y_ti_1_Si_good(i - (2/3)*NoofSamples);
        end
    end
end
figureIndex = figureIndex + 1;
figure(figureIndex)
plot(t1,y_ti_1_good,'k','LineWidth',3)
hold on
% ----------------------------------------------------------------------------------------------
y_ti_2_Fcl_good = step(Fcl_loc_2_good,t1)   ;
y_ti_2_So_good = step(So_MO_loc_2_good,t2)  ;
y_ti_2_Si_good = step(Si_MO_loc_2_good,t3)  ;

y_ti_2_good = zeros(NoofSamples + 1,1) ;
for i = 1:1:NoofSamples + 1
    y_ti_2_good(i) = y_ti_2_Fcl_good(i);
    if i > NoofSamples/3
        y_ti_2_good(i) = y_ti_2_good(i) + y_ti_2_So_good(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_ti_2_good(i) = y_ti_2_good(i) + y_ti_2_Si_good(i - (2/3)*NoofSamples);
        end
    end
end
plot(t1,y_ti_2_good,'r','LineWidth',3)
% ----------------------------------------------------------------------------------------------
x_min_fig3 = 0; 
x_max_fig3 = 3*FinalTime;
y_min_fig3 = 0;

sortMAX = [y_ti_1_Si_good y_ti_2_Si_good];
chooseMaxtemp = max(max(sortMAX));
chooseMax = max(2,1 + chooseMaxtemp);
if chooseMax == 2
    y_max_fig3 = 2.1;
else
  y_max_fig3 = 1.05*chooseMax;
end
axis([x_min_fig3 x_max_fig3 y_min_fig3 y_max_fig3])
grid
legend('F_{cl1} - T_{i1} tuning','F_{cl2} - T_{i2} tuning',0)
title('Automatic Tuning Second Order System - Good Open Loop Estimation Parameters')

% ----------------------------------------------------------------------------------------------
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% ----------------------------------------------------------------------------------------------
y_ti_1_Fcl = step(Fcl_loc_1,t1)   ;
y_ti_1_So = step(So_MO_loc_1,t2)  ;
y_ti_1_Si = step(Si_MO_loc_1,t3)  ;

y_ti_1 = zeros(NoofSamples + 1,1);
for i = 1:1:NoofSamples + 1
    y_ti_1(i) = y_ti_1_Fcl(i);
    if i > NoofSamples/3
        y_ti_1(i) = y_ti_1(i) + y_ti_1_So(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_ti_1(i) = y_ti_1(i) + y_ti_1_Si(i - (2/3)*NoofSamples);
        end
    end
end
figureIndex = figureIndex + 1;
figure(figureIndex)
plot(t1,y_ti_1,'k','LineWidth',3)
hold on
% ----------------------------------------------------------------------------------------------
y_ti_2_Fcl = step(Fcl_loc_2,t1)   ;
y_ti_2_So = step(So_MO_loc_2,t2)  ;
y_ti_2_Si = step(Si_MO_loc_2,t3)  ;

y_ti_2 = zeros(NoofSamples + 1,1) ;
for i = 1:1:NoofSamples + 1
    y_ti_2(i) = y_ti_2_Fcl(i)     ;
    if i > NoofSamples/3
        y_ti_2(i) = y_ti_2(i) + y_ti_2_So(i - NoofSamples/3);
        if i > (NoofSamples/3)*2
            y_ti_2(i) = y_ti_2(i) + y_ti_2_Si(i - (2/3)*NoofSamples);
        end
    end
end
plot(t1,y_ti_2,'r','LineWidth',3)
% ----------------------------------------------------------------------------------------------
x_min_fig3 = 0; 
x_max_fig3 = 3*FinalTime;
y_min_fig3 = 0;

sortMAX = [y_ti_1_Si y_ti_2_Si];
chooseMaxtemp = max(max(sortMAX));
chooseMax = max(2,1 + chooseMaxtemp);
if chooseMax == 2
    y_max_fig3 = 2.1;
else
  y_max_fig3 = 1.05*chooseMax;
end
axis([x_min_fig3 x_max_fig3 y_min_fig3 y_max_fig3])
grid
legend('F_{cl1} - T_{i1} tuning','F_{cl2} - T_{i2} tuning',0)
title('Automatic Tuning Second Order System - Bad Open Loop Estimation Parameters')

% ----------------------------------------------------------------------------------------------
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% ----------------------------------------------------------------------------------------------
figureIndex = figureIndex + 1;
figure(figureIndex)
plot(t1,y_ti_1_good,'k',t1,y_ti_1,'r','LineWidth',3)
grid on

figureIndex = figureIndex + 1;
figure(figureIndex)
plot(t1,y_ti_2_good,'k',t1,y_ti_2,'r','LineWidth',3)
grid on

ltiview(Fcl_loc_1_good,'k',Fcl_loc_1,'r')
ltiview(Fcl_loc_2_good,'k',Fcl_loc_2,'r')
ltiview(Si_MO_loc_2_good,'k',Si_MO_loc_2,'r')
ltiview(Gp_est_good,'k',Gp_est,'r')

% ----------------------------------------------------------------------------------------------
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% ----------------------------------------------------------------------------------------------
figureIndex = figureIndex + 1;
figure(figureIndex)
[y_1_good tspan] = step(Fcl_loc_1_good,tspan);
[y_1 tspan] = step(Fcl_loc_1,tspan)          ;
plot(tspan,y_1_good,'k',tspan,y_1,'r','LineWidth',3)
grid on
% ltiview(Fcl_loc_1_good,'k',Fcl_loc_1,'r')

figureIndex = figureIndex + 1;
figure(figureIndex)
[y_2_good tspan] = step(Fcl_loc_2_good,tspan);
[y_2 tspan] = step(Fcl_loc_2,tspan)          ;
plot(tspan,y_2_good,'k',tspan,y_2,'r','LineWidth',3)
grid on
% ltiview(Fcl_loc_2_good,'k',Fcl_loc_2,'r')

figureIndex = figureIndex + 1;
figure(figureIndex)
[y_Si_1_good tspan] = step(Si_MO_loc_1_good,tspan);
[y_Si_1 tspan] = step(Si_MO_loc_1,tspan)          ;
plot(tspan,y_Si_1_good,'k',tspan,y_Si_1,'r','LineWidth',3)
grid on
% ltiview(Si_MO_loc_1_good,'k',Si_MO_loc_1,'r')

figureIndex = figureIndex + 1;
figure(figureIndex)
[y_Si_2_good tspan] = step(Si_MO_loc_2_good,tspan);
[y_Si_2 tspan] = step(Si_MO_loc_2,tspan)          ;
plot(tspan,y_Si_2_good,'k',tspan,y_Si_2,'r','LineWidth',3)
grid on
% ltiview(Si_MO_loc_2_good,'k',Si_MO_loc_2,'r')