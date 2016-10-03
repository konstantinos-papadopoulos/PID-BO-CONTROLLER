
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


def create_plant(poles,zeros,time_delay,user_defined_plant):
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

