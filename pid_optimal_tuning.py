# Copyright (c) 2017 Konstantinos G. Papdopoulos. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the Eclipse Public License v1.0 which accompanies this distribution,
# and is available at http://www.eclipse.org/legal/epl-v10.html

class ControllerOptimal:
    """
    Creates an Optimal Controller object
    """
    def __init__(self, args):
        """
        Initializes the appropriate test component objects according to the
        test_type and the test configuration json object, in order to prepare
        the test for running

        :param args:
        :param json_conf:
        :param test_type:
        :type args:
        :type json_conf:
        :type test_type: str
        """
        self.plant.create_plant

    def autoTuneParamMainTspEstimation(args, plant):

        pass

    def pidOptimalTuningTypeI(args, plant):
        """Creates the controlled process determined by the user input
           Plant creation supports poles (max:5), zeros (max:5), time delay

        :param poles
        :param zeros
        :param time_delay
        :param user_defined_plant
        :
        :type
        :type
        :type
        :type
        """
w1 = n1 + p1            ;
w2 = n2 + n1*p1 +p2     ;
w3 = n2*p1 + n1*p2 + p3 ;
w4 = n2*p2 + n1*p3 + p4 ;
w5 = n2*p3 + n1*p4 + p5 ;
w6 = n2*p4 + n1*p5 + p6 ;
w7 = n2*p5 + n1*p6      ;
w8 = n2*p6              ;



A1 = w1^2- w2 - w1*q1 + q2                              ;
B1 = q1 -  w1                                           ;
C1 = w1^3 - w1^2*q1 - 2*w1*w2 + w2*q1 + w3 + w1*q2 - q3 ;

A2 = w2^2 - w2*q2 + w4 + w3*q1 + w1*q3 - 2*w1*w3 - q4 ;
B2 = w3 - w2*q1 + w1*q2 - q3 ;
C2 = (w1 - q1)*(w2^2 - 2*w1*w3 + 2*w4) - (w1*q4 + w3*q2 + w5 - w4*q1 + w2*q3);

C = [C1;C2]       ;
D = [A1 B1;A2 B2] ;
E = inv(D)        ;
F = E*C           ; 

x_MO  = 0         ;
y_MO  = 0         ;  
ti_MO = 2*kp*(w1 - q1 - x_MO) ;
ti_MOInit = ti_MO ;

num0Gc = 1                       ;
num1Gc = x_MO                    ;
num2Gc = y_MO                    ;
numGc_MO = [num2Gc num1Gc num0Gc];
    
den0Gc = 0                       ;
den1Gc = ti_MO                   ;
den2Gc = ti_MO*tp6               ;
denGc_MO = [den2Gc den1Gc den0Gc];
Gc_MO = tf(numGc_MO,denGc_MO)    ; 
Gc = Gc_MO                       ;
        pass

    def pidOptimalTuningTypePI(args, plant):
        """Creates the controlled process determined by the user input
           Plant creation supports poles (max:5), zeros (max:5), time delay

        :param poles
        :param zeros
        :param time_delay
        :param user_defined_plant
        :
        :type
        :type
        :type
        :type
        """

        pass

    def pidOptimalTuningTypePID(poles,zeros,time_delay,user_defined_plant):
        """Creates the controlled process determined by the user input
           Plant creation supports poles (max:5), zeros (max:5), time delay

        :param poles
        :param zeros
        :param time_delay
        :param user_defined_plant
        :
        :type
        :type
        :type
        :type
        """

        pass
