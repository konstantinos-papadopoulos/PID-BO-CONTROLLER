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
%  Contact:     kostas g. papadopoulos,    nikos e. mitrakis     konstantina mermkikli
%               kpapadop@eng.auth.gr       nmitr@ee.auth.gr      kmermikl@auth.gr
% 																										  																		
%  Place:	    Aristotle University of Thessaloniki, Thessaloniki, Greece							   																		
% 
% ************************************************************************************************************
function [figureIndexLocal] = train_anfis_pid_classifier(figureIndexLocal)
figureIndexLocal = figureIndexLocal + 1 ;
hold off
train_class = input('classifier_pi2pid_training - [y/n]:','s');
if strcmp(train_class,'y') == 1
    
    % Collect training and checking data
    % ------------------------------------------------------
    load classifier_train_pid.txt
    load classifier_check_pid.txt
    
    train_pid_class = classifier_train_pid;
    check_pid_class = classifier_check_pid;
    Xtrain = classifier_train_pid;
    Xcheck = classifier_check_pid;
    
    index = 1;
    k = 1    ; 
    fismat = genfis2(Xtrain(:,1:5),Xtrain(:,6),0.075); % We use the default values 
    [fismat1,error1] = anfis(Xtrain,fismat,[100])    ;
    fismat_classifier_pi2pid = fismat1               ;
    saveClass = input('Do you want to save the new trained classifier-  [y/n]:','s');
    if strcmp(train_class,'y') == 1
        save fismat_classifier_pi2pid fismat_classifier_pi2pid
    end
    trn_out_fismat1 = evalfis(Xtrain(:,1:5),fismat1);
else
    % Collect training and checking data
    % ------------------------------------------------------
    load classifier_train_pid.txt
    load classifier_check_pid.txt
    
    train_pid_class = classifier_train_pid;
    check_pid_class = classifier_check_pid;
    
    Xtrain = classifier_train_pid;
    Xcheck = classifier_check_pid;
    
    load fismat_classifier_pi2pid
    fismat1 = fismat_classifier_pi2pid;
 
    trn_out_fismat1 = evalfis(Xtrain(:,1:5),fismat1);
end

for i = 1:length(train_pid_class)
    if (trn_out_fismat1(i) >= 0.5)
        trn_out_fismat1(i) = 1;
    else
        trn_out_fismat1(i) = 0;
    end
end

% Calculate the 
% -------------------------------------------------------------------------
Perc = abs(trn_out_fismat1 - Xtrain(:,6));
total = sum(Perc)                        ;
CorrPercTrain = 100*(length(train_pid_class) - sum(Perc))/length(train_pid_class);

% Plot ANFIS training output versus desired output
% -------------------------------------------------------------------------
figureIndexLocal = figureIndexLocal + 1 ;
plot(Xtrain(:,6),'o')     ;
hold on                   ;
plot(trn_out_fismat1,'r*');
ylabel('Overshoot')
xlabel('Example')

% Plot ANFIS training output versus desired output for specific Tsc, a
% -------------------------------------------------------------------------
disp('Test ALL train plants')
reply = 'Y' ;
k = 2       ;
while (reply == 'Y')||(reply == 'y')
    reply = input('Do you want to test ANFIS output? Y/N: ','s');
    if ( reply == 'n' )||( reply == 'N' )
        break ;
    end

    desiredTsc = input('Enter value of Tsc (0.01,0.03,...,0.91), Tsc = ');
    desiredA = input('Enter value of a (0.01, 0.135, 0.26, 0.61), a = ') ;
    indexTsc = int16(((desiredTsc-0.01)/0.02)*84)                        ;
    
    if (desiredA == 0.01)
        indexA = 0  ;
    elseif (desiredA == 0.135)
        indexA = 21 ;
    elseif (desiredA == 0.26)
        indexA = 42 ;    
    else
        indexA = 63 ;        
    end
    
    anfis_in(1:21,1:5) = Xtrain(indexTsc+indexA+1:indexTsc+indexA+21,1:5);
    test_out_fismat1 = evalfis(anfis_in,fismat1);
    
    for i=1:21
        if (test_out_fismat1(i) >= 0.5)
            test_out_fismat1(i) = 1;
        else
            test_out_fismat1(i) = 0;
        end
    end
                    
    CorrPerf = (21 - sum(abs(Xtrain(indexTsc+indexA+1:indexTsc+indexA+21,6)-test_out_fismat1)))/21*100;
    figureIndexLocal = figureIndexLocal + 1 ;
    figure(figureIndexLocal)
    plot(Xtrain(indexTsc+indexA+1:indexTsc+indexA+21,6),'o')
    hold on;
    plot(test_out_fismat1,'r*')
    title(['Testing Examples - Tsc=',num2str(desiredTsc),', Performance=',num2str(CorrPerf),'%'])
    ylabel('Classifier Existence')
    xlabel('Example')
    k = k + 1 ;
end

% Plot ANFIS training output versus desired output for specific Tsc, a
% -------------------------------------------------------------------------
disp('Test ALL check plants')

test_all_out_fismat1 = evalfis(Xcheck(:,1:5),fismat1);
for i = 1:length(check_pid_class)
    if (test_all_out_fismat1(i) >= 0.5)
        test_all_out_fismat1(i) = 1;
    else
        test_all_out_fismat1(i) = 0;
    end
end

Perc = abs(test_all_out_fismat1 - Xcheck(:,6));
total = sum(Perc)                           ;
CorrPercCheck = 100*(length(check_pid_class)-sum(Perc))/length(check_pid_class);


disp('Check the checking systems')

reply='Y';
k = 2    ;
while (reply == 'Y')||(reply == 'y')
    reply = input('Do you want to test ANFIS output? Y/N: ','s');
    if ( reply == 'n' )||( reply == 'N' )
        break ;
    end
    desiredTsc = input('Enter value of Tsc (0.01,0.03,...,0.91), Tsc=');
    desiredA = input('Enter value of a (0.035, 0.06, 0.085, 0.11, 0.16, 0.185, 0.21, 0.235, 0.285, 0.31, 0.41, 0.51, 0.71, 0.81, 0.91), a=');
    indexTsc = int16(((desiredTsc-0.01)/0.02)*84);
    if (desiredA==0.035)
        indexA = 0 ;
    elseif (desiredA == 0.06)
        indexA = 21;
    elseif (desiredA == 0.085)
        indexA = 42;    
    elseif (desiredA == 0.11)
        indexA = 63;    
    elseif (desiredA == 0.16)
        indexA = 84;    
    elseif (desiredA == 0.185)
        indexA = 105;    
    elseif (desiredA == 0.21)
        indexA = 126;    
    elseif (desiredA == 0.235)
        indexA = 147;    
    elseif (desiredA == 0.285)
        indexA = 168;            
    elseif (desiredA == 0.31)
        indexA = 189;
    elseif (desiredA == 0.41)
        indexA = 210;    
    elseif (desiredA == 0.51)
        indexA = 231;    
    elseif (desiredA == 0.71)
        indexA = 252;    
    elseif (desiredA == 0.81)
        indexA = 273;    
    else
        indexA = 294;    
    end
    
    anfis_in(1:21,1:5) = Xcheck(indexTsc + indexA + 1:indexTsc + indexA + 21,1:5);
    test_out_fismat1 = evalfis(anfis_in,fismat1);
    for i = 1:21
        if (test_out_fismat1(i) >= 0.5)
            test_out_fismat1(i) = 1;
        else
            test_out_fismat1(i) = 0;
        end
    end

    CorrPerf = (21-sum(abs(Xtrain(indexTsc+indexA+1:indexTsc+indexA+21,6)-test_out_fismat1)))/21*100;
    figureIndexLocal = figureIndexLocal + 1 ;
    figure(figureIndexLocal)
    plot(Xtrain(indexTsc + indexA + 1:indexTsc + indexA + 21,6),'o')
    hold on          ;
    plot(test_out_fismat1,'r*')
    title(['Testing Examples - Tsc=',num2str(desiredTsc),', Performance=',num2str(CorrPerf),'%'])
    ylabel('Classifier Existence')
    xlabel('Example')
    k = k + 1        ;

end

