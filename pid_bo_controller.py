# Copyright (c) 2016 Konstantinos G. Papdopoulos. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the Eclipse Public License v1.0 which accompanies this distribution,
# and is available at http://www.eclipse.org/legal/epl-v10.html




"""
Main configuration script for choosing the control loop design
-The control loop can be of type-I, type-II, type-III, type-V
-The control loop design can be either of analog/digital
-The plant to be controlled (process) can be determined manually by the user
 or generated manually.
"""

import argparse
import plant


def main():
    """This is the main function where the main test application starts.
    """

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--control-loop-type',
                        required=True,
                        type=str,
                        dest='control_loop_type',
                        action='store',
                        help="type_I   (for zero steady state position error)\n"
                             "type_II  (for zero steady state position/velocity error)\n"
                             "type_III (for zero steady state position/velocity/acceleration error)\n"
                             "type_IV  (for elimination higher order errors )"
                             )
    parser.add_argument('--plant',
                        required=True,
                        type=str,
                        dest='plant',
                        action='store',
                        help="random\n"
                             "determined"
                             )
    parser.add_argument('--poles',
                        required=True,
                        type=str,
                        dest='poles',
                        action='store',
                        help="number of poles"
                             )
    parser.add_argument('--zeros',
                        required=True,
                        type=str,
                        dest='zeros',
                        action='store',
                        help="number of zeros"
                             )
    parser.add_argument('--time-delay',
                        required=True,
                        type=str,
                        dest='time_delay',
                        action='store',
                        help="time delay of the controlled process"
                             )
    parser.add_argument('--design_type',
                        required=False,
                        type=str,
                        dest='plant',
                        action='store',
                        help="analog\n"
                             "digital"
                             )
    parser.add_argument('--log-file',
                        dest='log_file',
                        action='store',
                        help='log file name')
    parser.add_argument('--output-dir',
                        required=False,
                        type=str,
                        dest='output_dir',
                        action='store',
                        help='result files output directory ')
    parser.add_argument('--logging-level',
                        type=str,
                        dest='logging_level',
                        action='store',
                        default='DEBUG',
                        help="log level set."
                             "possible values are:\n"
                             "INFO\n"
                             "DEBUG (default)\n"
                             "ERROR")

    args = parser.parse_args()

    var = input("Please enter something: ")
    print(var)
    print(args.poles)
    plant.create_plant(poles, zeros, time_delay, user_defined_plant)


if __name__ == '__main__':
    main()
