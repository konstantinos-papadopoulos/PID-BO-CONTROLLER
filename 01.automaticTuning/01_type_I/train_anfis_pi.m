%Desired PI Overshoot Prediction for OSK
close all;
clear all;
clc;

load ov_pi.txt

%Get the number of examples
numExamples = size(ov_pi);

%Formulate the training data set
index = 1;
k = 1    ;

%We use the Tsc=0.01,0.11,0.21...0.81 values and a=0.1,0.12,0.14,...0.9 for
%each Tsc to formulate the training data set
while (index<numExamples(1))
    for i=1:1:82
        Xtrain(k,:) = ov_pi(index,:);
        k = k + 1;
        index = index + 1;
    end
    index = index + 82*9;
end

%Use the Substractive Clustering method to generate a new Sugeno type fuzzy
%system
disp('Generating FIS initial architecture....')
%Grid partition with Gaussian 3 MF and linear type consequent
fismat = genfis1(Xtrain(:,1:3),[3 3],'gaussmf','linear');

%%%fismat=genfis2(Xtrain(:,1:2),Xtrain(:,3),0.5); % We use the default values 
%of Jiang for the method, fismat contains the initial fuzzy system

%Train the system using the anfis editor
[fismat1,error1] = anfis(Xtrain,fismat,[100]);

%Calculate the ouput of ANFIS for the training data set
trn_out_fismat1 = evalfis(Xtrain(:,1:2),fismat1);

%Plot ANFIS training output versus desired output
plot(Xtrain(:,3),'o')
hold on;
plot(trn_out_fismat1,'r*')
title(['Training Examples - SSE=',num2str(error1(end))])
ylabel('Overshoot')
xlabel('Example')

%Choose the values of Tsc to plot the anfis output
reply ='Y';
k = 2;
while (reply == 'Y')||(reply == 'y')
    reply = input('Do you want to test ANFIS output? Y/N: ','s');
    if (reply=='n')||(reply=='N')
        break;
    end
    desiredTsc = input('Enter value of Tsc (0.01,0.02,...,0.9), Tsc=');
    index = (desiredTsc-0.01)*100*82;
    Xtest(1:82,1:2) = ov_pi(index+1:index+82,1:2);
    test_out_fismat1 = evalfis(Xtest,fismat1);
    SSE = sum((ov_pi(index+1:index+82,3)-test_out_fismat1).^2);
    figure(k)
    plot(ov_pi(index+1:index+82,3),'o')
    hold on;
    plot(test_out_fismat1,'r*')
    title(['Testing Examples - Tsc=',num2str(desiredTsc),',SSE=',num2str(SSE)])
    ylabel('Overshoot')
    xlabel('Example')
    k = k + 1;
end

disp('The control procedure now is to follow.....')



