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
