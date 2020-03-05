
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
    polesOrder = args.poles
    zerosOrder = args.zeros
    time_delay = args.time_delay
    poles = []; polesSorted = [];
    zeros = []; zerosSorted = [];
    for i in range(int(polesOrder)+1):
        poles.append(random.random())
    for j in range(int(zerosOrder)+1):
        zeros.append(random.random())

    polesSorted = sorted(poles, reverse=True)
    zerosSorted = sorted(zeros, reverse=True)
    plant = control.TransferFunction(zerosSorted, polesSorted)
    (T , yout) = control.step_response(plant)
    plot_step_response(T, yout)
    
    return (T , yout)

def plot_step_response(T, yout):
    
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
