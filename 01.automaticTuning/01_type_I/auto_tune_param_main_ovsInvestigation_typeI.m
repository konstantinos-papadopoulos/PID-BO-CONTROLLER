% *********************************************************************************************************** 
%               Copyright (C) 2007
%               Aristotle University of Thessaloniki
%               Depaertment of Electrical & Computer Engineering
%               Division of Electronics & Computer Engineering
% 
% ************************************************************************************************************
%  Title:       auto_tune_param_main.m																			   																		  	
%  Project:     Automatic tuning of the paramters for PI,PID controllers
%  
%  Purpose:     plant implementation																		   																		
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
clc
clear all
disp('-------------------------------------------------------------------')
disp('enter "y" in order to create data...')
disp('-------------------------------------------------------------------')
createData = input('y/n:','s');
switch lower(createData)
    case('y')
    % *************************************************************************
    scaleInit = 0.01;
    scaleStep = 0.01;
    scaleEnd = 0.91 ;
    % *************************************************************************

    include = 0;

    step = 0.01              ;
    a = 0.01:step:0.9 + step  ;
    b = 0.01:step:0.9 + step  ;

    count = 0 ;
    k = 0     ;
    l = 0     ;
    m = 0     ;
    scale = scaleInit;
    while scale <= scaleEnd
        % matrix initialization 
        % *********************************************************************
        alpha = zeros(length(a),1);
        beta = zeros(length(b),1) ;
        ovs = zeros(length(a),length(b)) ;
        settlingTime = zeros(length(a),length(b));
        
        ovs_pid = zeros(length(a),length(a));
        ti_pid = zeros(length(a),length(a));
        x_pid = zeros(length(a),length(a));
        y_pid = zeros(length(a),length(a));
                           
        ovs_pi = zeros(length(a),length(a));
        ti_pi = zeros(length(a),length(a));
        x_pi = zeros(length(a),length(a));
        
        ovs_i = zeros(length(a),length(a));
        ti_i = zeros(length(a),length(a));

        % *********************************************************************

        for i = 1:length(a)

          if mod(i,1) == 0
                fprintf('stepNo:%d of %d - scale:%d\n',i,length(a),scale);
          end

            for j = 1:length(b)

                % PlantGeneration
                %-------------------------------------------
                kp = 1 ;
                kh = 1 ; 

                Tz1 = 0          ;  Tp1 = 1       ;  
                Tz2 = 0*b(j)     ;  Tp2 = a(i)    ;
                Tz3 = 0*b(j)^2   ;  Tp3 = a(i)^2  ;
                Tz4 = 0*b(j)^3   ;  Tp4 = a(i)^3  ;
                                    Tp5 = 0*rand  ;
                Td = include*b(j);  Tp6 = scale*Tp1 ;
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
                tp6 = Tp6 / Tp1 ; td = Td / Tp1 ;

                plant.tp1 = tp1;    plant.tz1 = tz1;
                plant.tp2 = tp2;    plant.tz2 = tz2;
                plant.tp3 = tp3;    plant.tz3 = tz3;
                plant.tp4 = tp4;    plant.tz4 = tz4;
                plant.tp5 = tp5;
                plant.tp6 = tp6;
                plant.td = td  ;    plant.kp = kp  ; plant.kh = kh  ;
%                            %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                        
%                            %PID Controller
%                            %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                        
%                            [Gc_pid Gp stepInformation_pid ti_MO_pid x_MO_pid y_MO_pid stepinformation_Gp] = auto_tune_param_main_pid(plant);
                             SettlingTimePlant = auto_tune_param_main_SettlingTimePlant(plant);
                             settlingtimeFinalPlant = SettlingTimePlant.SettlingTime;
%                            Gc = Gc_pid;                        
%                            ovs_pid(i,j) = stepInformation_pid.Overshoot;
%                            ti_pid(i,j) = ti_MO_pid ;
%                            x_pid(i,j) = x_MO_pid;
%                            y_pid(i,j) = y_MO_pid;
%                            if isnan(ovs_pid(i,j)) == 1
%                                ovs_pid(i,j) = -1;
%                            end
%                            %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                        
%                            %PI Controller
%                            %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                        
%                            [Gc_pi Gp stepInformation_pi ti_MO_pi x_MO_pi] = auto_tune_param_main_pi(plant)  ;
%                            ovs_pi(i,j) = stepInformation_pi.Overshoot;
%                            ti_pi(i,j) = ti_MO_pi;
%                            x_pi(i,j) = x_MO_pi;
%                            if isnan(ovs_pi(i,j)) == 1
%                                ovs_pi(i,j) = -1;
%                            end
%                            %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                        
%                            %I Controller
%                            %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                        
%                            [Gc_i Gp stepInformation_i ti_MO_i] = auto_tune_param_main_i(plant)    ;
%                            ovs_i(i,j) = stepInformation_i.Overshoot;
%                            ti_i(i,j) = ti_MO_i;
%                            if isnan(ovs(i,j)) == 1
%                                ovs_i(i,j) = -1;
%                            end
                           settlingtimeGp(i,j) = settlingtimeFinalPlant;
                           
            end

        end
                    FilenameNo = scale ;
                    FilenameNo2str = num2str(scale); 
                    % -------------------------------------------------------------------------


                    fid = fopen(filename,'wt');     
                        for k = 1:length(a)
                            fprintf(fid,'%1.5f\n',settlingtimeGpFinal(k));
                        end   
                    fclose(fid);  

        % storing the Final Values [ovs,ti]
        % -------------------------------------------------------------------------
        ovsFinal_pid = ovs_pid(:,1);
        tiFinal_pid = ti_pid(:,1)  ;
        xfinal_pid = x_pid(:,1)    ;
        yfinal_pid = x_pid(:,1)    ;
        
        ovsFinal_pi = ovs_pi(:,1);
        tiFinal_pi = ti_pi(:,1)  ;
        xfinal_pi = x_pi(:,1)    ;
        
        ovsFinal_i = ovs_i(:,1)  ;
        tiFinal_i = ti_i(:,1)    ;
        % Measuring the Plant Sampling Time
        % -------------------------------------------------------------------------
        settlingtimeGpFinal = settlingtimeGp(:,1)            ;
        
        FilenameNo = scale ;
        FilenameNo2str = num2str(scale); 
        
        % -------------------------------------------------------------------------
%         filename = strcat('data_',FilenameNo2str,'.txt');
        filename = strcat('settlingTime_',FilenameNo2str,'.txt');
        fid = fopen(filename,'wt');     
            for i = 1:length(a)
%                 fprintf(fid,'%1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f\n',ovsFinal_i(i),tiFinal_i(i),ovsFinal_pi(i),tiFinal_pi(i),xfinal_pi(i),ovsFinal_pid(i),tiFinal_pid(i),xfinal_pid(i),yfinal_pid(i));
                fprintf(fid,'%1.5f\n',settlingtimeGpFinal(k));
            end   
        fclose(fid);   

        scale = scale + scaleStep;
    end
    otherwise
       
        % sum up the whole data
        % -------------------------------------------------------------------------
        scaleInit = 0.01;
        scaleStep = 0.01;
        scaleEnd = 0.9 ;
        scaleMatrix = scaleInit:scaleStep:scaleEnd;
        stepA = 0.01               ;
        a = 0.01:stepA:0.9 + stepA  ;
        LoopEnd = length(scaleMatrix)-1;

            for i = 1:LoopEnd
                    FilenameNo2str = num2str(scaleMatrix(i)); 
                    filename = strcat('data_',FilenameNo2str,'.txt');
                    load(filename)
                    ovs_I = data_0(:,1) ;
                    ovs_PI = data_0(:,3);
                    ti_I = data_0(:,2)    ;
                    figure(1)
                    plot(a,ovs_I)
                    hold on
                    figure(2)
                    plot(a,ovs_PI)
                    hold on
            end
            figure(1)
            grid on
            xlabel('a')
            ylabel('ovs_I')
            figure(2)
            grid on
            xlabel('a')
            ylabel('ovs_PI')
            
            %creating the Final Data for Training
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            counter = 0;
            loopEnd = length(scaleMatrix);

            % Matrix Initialization
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          
            ovs_IFinal  = zeros((loopEnd-1)*length(ovs_I),1);
            ti_IFinal   = zeros((loopEnd-1)*length(ovs_I),1);
            ovs_PIFinal = zeros((loopEnd-1)*length(ovs_I),1);
            ti_PIFinal  = zeros((loopEnd-1)*length(ovs_I),1);
            x_PIFinal   = zeros((loopEnd-1)*length(ovs_I),1);
            
            ovs_PIDFinal = zeros((loopEnd-1)*length(ovs_I),1);
            ti_PIDFinal  = zeros((loopEnd-1)*length(ovs_I),1);
            x_PIDFinal   = zeros((loopEnd-1)*length(ovs_I),1);
            y_PIDFinal   = zeros((loopEnd-1)*length(ovs_I),1);

            for i = 1:loopEnd
                FilenameNo2str = num2str(scaleMatrix(i)); 
                filename = strcat('data_',FilenameNo2str,'.txt');
                [ovs_I ti_I ovs_PI ti_PI x_PI ovs_PID ti_PID x_PID y_PID] = textread(filename);
                  for j = 1:length(ovs_I)
                      ovs_IFinal(j + counter) = ovs_I(j)    ;
                      ti_IFinal(j + counter) = ti_I(j)      ;
                      
                      ovs_PIFinal(j + counter) = ovs_PI(j)  ;
                      ti_PIFinal(j + counter) = ti_PI(j)    ;
                      x_PIFinal(j + counter) = x_PI(j)      ;

                      ovs_PIDFinal(j + counter) = ovs_PID(j);
                      ti_PIDFinal(j + counter) = ti_PID(j)  ;
                      x_PIDFinal(j + counter) = x_PID(j)    ;
                      y_PIDFinal(j + counter) = y_PID(j)    ;
                  end
                 counter = counter + length(ovs_I);  
            end
            %writing to Final Data
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            finalData = 'FinalData.txt';
            fid = fopen(finalData,'wt');     
            for k = 1:length(ti_IFinal)
                fprintf(fid,'%1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f\n',ovs_IFinal(k),ti_IFinal(k),ovs_PIFinal(k),ti_PIFinal(k),x_PIFinal(k),ovs_PIDFinal(k),ti_PIDFinal(k),x_PIDFinal(k),y_PIDFinal(k));
            end   
            fclose(fid);  
end   
