% Diploma thesis: Section: 2.3.1.1.
% ----------------------------------------------------------------
clc
clear all
close all
% *************************************************************************
inclDelStr = input('inclDel - y/n:','s');
if strcmp(inclDelStr,'y') == 1
    inclDel = 1;
    scaleDel = input('scaleDel = ');
else
    inclDel = 0 ;
    scaleDel = 0;
end
inclZerStr = input('inclZer - y/n:','s');
if strcmp(inclZerStr,'y') == 1
    inclZer = 1               ;
    zeroValue = 0.45          ;
    b = inclZer*zeroValue     ;
else
    inclZer = 0          ;
    zeroValue = 0        ; 
    b = inclZer*zeroValue;
end
    stepA = 0.02                ;
    stepB = 0.02                ;
    VarStep= input('variable Step - y/n:','s');
if strcmp(VarStep,'y') == 1
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    a = [0.01 0.1 0.3 0.5 0.7 0.9];
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
else
    a = 0.01:stepA:0.9 + stepA ;
end

tn_init = 7.7            ;
plant.tn = tn_init       ;
scale = input('scale = ');
NoofSamples = 3000       ;
y1 = zeros(NoofSamples,1)  ; t1 = zeros(NoofSamples,1);
y2 = zeros(NoofSamples,1)  ; t2 = zeros(NoofSamples,1);
y3 = zeros(NoofSamples,1)  ; t3 = zeros(NoofSamples,1);
y4 = zeros(NoofSamples,1)  ; t4 = zeros(NoofSamples,1);
y5 = zeros(NoofSamples,1)  ; t5 = zeros(NoofSamples,1);
y6 = zeros(NoofSamples,1)  ; t6 = zeros(NoofSamples,1);
tFinal = 0:0.001:50;  

% Initialize the optimal matrix
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       for i = 1:length(a)
          if mod(i,1) == 0
             fprintf('stepNo:%d of %d - scale:%d\n',i,length(a),scale);
          end

                kp = 1 ;
                kh = 1 ; 

                Tz1 = inclZer*b     ;  Tp1 = 1         ;   
                Tz2 = 0*b^2         ;  Tp2 = a(i)      ;
                Tz3 = inclZer*b^3   ;  Tp3 = a(i)^2    ;
                Tz4 = 0*inclZer*b^4 ;  Tp4 = a(i)^3    ;
                                       Tp5 = a(i)^4    ;
                                       Tp6 = scale*Tp1 ;
                Td = scaleDel       ;
                % -------------------------------------------------------------------------
                % Normalizing zeros/poles with Tp1
                % -------------------------------------------------------------------------
                tz1 = Tz1 / Tp1 ;       tp1 = Tp1 / Tp1 ;
                tz2 = Tz2 / Tp1 ;       tp2 = Tp2 / Tp1 ;
                tz3 = Tz3 / Tp1 ;       tp3 = Tp3 / Tp1 ;
                tz4 = Tz4 / Tp1 ;       tp4 = Tp4 / Tp1 ;
                                        tp5 = Tp5 / Tp1 ;
                                        
                % -------------------------------------------------------------------------
                % Normalizing the controller unmodelled dynamic/ time delay with Tp1
                % -------------------------------------------------------------------------
                tp6 = Tp6 / Tp1;    td = Td / Tp1  ;
                plant.tp1 = tp1;    plant.tz1 = tz1;
                plant.tp2 = tp2;    plant.tz2 = tz2;
                plant.tp3 = tp3;    plant.tz3 = tz3;
                plant.tp4 = tp4;    plant.tz4 = tz4;
                plant.tp5 = tp5;
                plant.tp6 = tp6;
                plant.td = td  ;    plant.kp = kp  ; plant.kh = kh  ;
                [Fcl_MO Fcl_optimal Gp_loc stepinformation ti_MO x_MO y_MO] = autotune_param_typeII_calculateOvs(plant);
                Gp = Gp_loc;

                % Logging the Optimal Controller
                % *********************************************
                Gp_log{i} = Gp;
                Fcl_optimal_PI{i} = Fcl_optimal                                               ;
                So_MO_optimal_PI{i}  = 1 - Fcl_optimal_PI{i}                                  ;
                Si_MO_optimal_PI{i}  = series(So_MO_optimal_PI{i},Gp_log{i})                  ;
                Fcl_optimal_PI_struct{i} = stepinfo(Fcl_optimal_PI{i})                        ;
                settlingtime_Fcl_optimal_PI_struct(i) = Fcl_optimal_PI_struct{i}.SettlingTime ;
                FinalTime_Fcl_optimal_PI(i) = 2*settlingtime_Fcl_optimal_PI_struct(i)         ;
                
                if a(i) > 0 && a(i) <=  0.0100
                    [y1 t1] = step(Fcl_optimal_PI{i},tFinal);
                    [y1_di t1] = step(Si_MO_optimal_PI{i},tFinal);
                elseif a(i) > 0.0100 && a(i) <=  0.1000
                    [y2 t2] = step(Fcl_optimal_PI{i},tFinal);
                    [y2_di t1] = step(Si_MO_optimal_PI{i},tFinal);
                elseif a(i) >0.1000 && a(i) <= 0.3000 
                    [y3 t3] = step(Fcl_optimal_PI{i},tFinal);
                    [y3_di t1] = step(Si_MO_optimal_PI{i},tFinal);
                elseif a(i) >0.3000 && a(i) <= 0.5000 
                    [y4 t4] = step(Fcl_optimal_PI{i},tFinal);
                    [y4_di t1] = step(Si_MO_optimal_PI{i},tFinal);
                elseif a(i) >0.5000 && a(i) <= 0.7000
                    [y5 t5] = step(Fcl_optimal_PI{i},tFinal);
                    [y5_di t1] = step(Si_MO_optimal_PI{i},tFinal);
                else  %a >0.7000 && a <=  0.9000
                    [y6 t6] = step(Fcl_optimal_PI{i},tFinal);
                    [y6_di t1] = step(Si_MO_optimal_PI{i},tFinal);
                end
       end
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figure(1)     
plot(t1,y1,'k','LineWidth',1.5)
hold on

plot(t1,y2,'-.k','LineWidth',1.5)
hold on

plot(t1,y3,'r','LineWidth',1.5)
hold on

plot(t1,y4,'-.r','LineWidth',1.5)
hold on

plot(t1,y5,'b','LineWidth',1.5)
hold on

plot(t1,y6,'-.b','LineWidth',1.5)
hold on

grid on
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figure(2)     
plot(t1,y1_di,'k','LineWidth',1.5)
hold on

plot(t1,y2_di,'-.k','LineWidth',1.5)
hold on

plot(t1,y3_di,'r','LineWidth',1.5)
hold on

plot(t1,y4_di,'-.r','LineWidth',1.5)
hold on

plot(t1,y5_di,'b','LineWidth',1.5)
hold on

plot(t1,y6_di,'-.b','LineWidth',1.5)
hold on

grid on
           
    
    
    
    
    
