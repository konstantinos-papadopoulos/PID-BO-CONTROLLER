
# Copyright (c) 2016 Konstantinos G. Papdopoulos. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the Eclipse Public License v1.0 which accompanies this distribution,
# and is available at http://www.eclipse.org/legal/epl-v10.html

"""
Controlled process
"""

import control
import logging
import matplotlib
import matplotlib.pyplot as plt
import random

def create_plant(args):
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
    kp = random.random();
    polesOrder = args.poles
    zerosOrder = args.zeros
    time_delay = args.time_delay
    poles = []; polesSorted = [];
    zeros = []; zerosSorted = [];
    for i in range(int(polesOrder)+1):
        poles.append(random.random())
    for j in range(int(zerosOrder)+1):
        zeros.append(random.random())
    
    # Sort poles and zeros to identify dominant pole and dominant zeros
    # --------------------------------------------------------------------------
    polesSorted = sorted(poles, reverse=True)
    zerosSorted = sorted(zeros, reverse=True)

    # plant zeros coefficients
    # TODO: create a dynamic dictionary to tz_i based on the order of zeros 
    # coming from the user
    # --------------------------------------------------------------------------
    tz1 = zerosSorted[0]; tz2 = 0; tz3 = 0; 
    tz4 = 0; 
    
    q0 = 1                                                         ;
    q1 = tz1 + tz2 + tz3 + tz4                                     ;
    q2 = tz1*tz2 + tz1*tz3 + tz1*tz4 + tz2*tz3 + tz2*tz4 + tz3*tz4 ;
    q3 = tz1*tz2*tz3 + tz1*tz2*tz4 + tz2*tz3*tz4 + tz1*tz3*tz4     ;
    q4 = tz1*tz2*tz3*tz4 

    # plant poles coefficients dynamics
    # TODO: create a dynamic dictionary to tp_i based on the order of zeros 
    # coming from the user    
    # --------------------------------------------------------------------------
    tp1 = polesSorted[0]; tp2 = polesSorted[1]; tp3 = polesSorted[2]; 
    tp4 = polesSorted[3]; tp5 = polesSorted[4];

    pGp0 = 1; 
    pGp1 = sum(polesSorted); # tp1 + tp2 + tp3 + tp4 + tp5
    pGp2 = tp1*tp2 + tp1*tp3 + tp1*tp4 + tp1*tp5 + tp2*tp3 + \
           tp2*tp4 + tp2*tp5 + tp3*tp4 + tp3*tp5 + tp4*tp5;

    pGp3 = tp1*tp2*tp3 + tp1*tp2*tp4 + tp1*tp2*tp5 + tp1*tp3*tp4 + \
           tp1*tp3*tp5 + tp1*tp4*tp5 + tp2*tp3*tp4 + tp2*tp3*tp5 + tp2*tp4*tp5 \
           + tp3*tp4*tp5;

    pGp4 = tp1*tp2*tp3*tp4 + tp1*tp2*tp3*tp5 + tp1*tp2*tp4*tp5 + \
           tp1*tp3*tp4*tp5 + tp2*tp3*tp4*tp5;
    pGp5 = tp1*tp2*tp3*tp4*tp5

    numGp = [q4, q3, q2, q1, q0];
    denGp  = [pGp5, pGp4, pGp3, pGp2, pGp1, pGp0];
    plantGp = control.tf(numGp,denGp);

    
    (T , yout) = control.step_response(plantGp)
    plotStepResponse(T, yout)
    
    return plantGp

def plantRandom(args):
    alpha = 0.1

    pass

def plantDetermined(args):

    pass

def plantAlpha(args):
    
    pass

def plotStepResponse(T, yout):
    
    fig, ax = plt.subplots()
    ax.plot(T, yout)
    ax.set(xlabel='time (s)', ylabel='$y_{out}$',
           title='Plant $G_p$ Open Loop step response')
    plt.grid(color='k', linestyle='-.', linewidth=0.1)
    plt.show()


def plant_estimation(poles,zeros,time_delay,user_defined_plant):
    """The method takes the real plant as an input and returns the settling
       time of the plant's step response, the estimated DC gain (plant's steady
       state gain, an estimation of the plant's unmodelled dynamics and the
       plant transfer function [first order mode] based on the step response
       measurements). All above measurements will be used as an input for
       initializing the automatic tuning algorith,

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

def plant_normalize():
    """The method takes the real plant as an input and returns the settling
       time of the plant's step response, the estimated DC gain (plant's steady
       state gain, an estimation of the plant's unmodelled dynamics and the
       plant transfer function [first order mode] based on the step response
       measurements). All above measurements will be used as an input for
       initializing the automatic tuning algorith,

       Plant creation supports poles (max:5), zeros (max:5), time delay
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
