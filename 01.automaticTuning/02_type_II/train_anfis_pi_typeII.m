% Desired PI Overshoot Prediction for OSK
% ------------------------------------------------------
close all;
clear all;
clc;

load ov_pi_check.txt
load ov_pi_train.txt

% Get the number of examples
% ------------------------------------------------------
Xtrain = ov_pi_train;
Xcheck = ov_pi_check;


% numExamples=size(ov_pi3);

% Formulate the training data set
% ------------------------------------------------------
index = 1;
k = 1;


%Use the Substractive Clustering method to generate a new Sugeno type fuzzy
%system
% ------------------------------------------------------
disp('Generating FIS initial architecture....')

%Grid partition with Gaussian 3 MF and linear type consequent
%fismat=genfis1(Xtrain(:,1:3),[3 3],'gaussmf','linear');

fismat = genfis2(Xtrain(:,1:4),Xtrain(:,5),0.5); % 
We use the default values 
%of Jiang for the method, fismat contains the initial fuzzy system

%Train the system using the anfis editor
% ------------------------------------------------------
[fismat1,error1] = anfis(Xtrain,fismat,[100]);

%Calculate the ouput of ANFIS for the training data set
% ------------------------------------------------------
trn_out_fismat1=evalfis(Xtrain(:,1:4),fismat1);

%Plot ANFIS training output versus desired output
% ------------------------------------------------------
plot(Xtrain(:,5),'o')
hold on;
plot(trn_out_fismat1,'r*')
title(['Training Examples - SSE=',num2str(error1(end))])
ylabel('Overshoot')
xlabel('Example')


%Choose the values of Tsc to plot the anfis output
% ------------------------------------------------------
reply ='Y';
k = 2;
while ( reply == 'Y' )||( reply == 'y' )
    reply = input('Do you want to test ANFIS output? Y/N: ','s');
    if ( reply =='n' )||( reply == 'N' )
        break;
    end
    desiredTsc = input('Enter value of Tsc (0.01,0.03,...,0.91), Tsc=');
    index=int16(((desiredTsc - 0.01)/0.02)*21);
    anfis_in(1:21,1:4) = Xcheck(index + 1:index + 21,1:4);
    test_out_fismat1 = evalfis(anfis_in,fismat1);

    % calculate Sum Square Error     
    % ------------------------------------------------------
    SSE = sum((Xcheck(index+1:index+21,5) - test_out_fismat1).^2);
    figure(k)
    plot(Xcheck(index+1:index+21,5),'o')
    hold on;
    plot(test_out_fismat1,'r*')
    title(['Testing Examples - Tsc=',num2str(desiredTsc),',SSE=',num2str(SSE)])
    ylabel('Overshoot')
    xlabel('Example')
    k = k + 1;
end
