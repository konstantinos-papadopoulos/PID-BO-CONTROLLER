#!/bin/bash

# Copyright (c) 2017 Konstantinos G. Papadopoulos. All rights reserved.
# -----------------------------------------------------------------------------
# This program and the accompanying materials are made available under the
# terms of the Eclipse Public License v1.0 which accompanies this distribution,
# and is available at http://www.eclipse.org/legal/epl-v10.html

BASE_DIR="/opt"
VENV_DIR_PID_BO_CONTROLLER="venv_pid_bo_controller"

# PROXY value is passed either from Vagrantfile/Dockerfile
#------------------------------------------------------------------------------
if [ ! -z "$1" ]; then
    echo "Creating PROXY variable: $1"
    PROXY=$1
else
    echo "Empty PROXY variable"
    PROXY=""
fi
# Python installation and other necessary libraries for pip
#------------------------------------------------------------------------------
apt-get update && apt-get install -y \
    python-pip \
    python3-pip \
    python-virtualenv \
    python3-matplotlib

# Configure pip optionspip
#------------------------------------------------------------------------------
pip_options=""

if [ ! -z "$PROXY" ]; then
    pip_options=" --proxy==$PROXY $pip_options"
fi

pip3 $pip_options install --upgrade pip

# NSTAT node provisioning actions
#------------------------------------------------------------------------------
mkdir $BASE_DIR/$VENV_DIR_PID_BO_CONTROLLER
virtualenv --system-site-packages $BASE_DIR/$VENV_DIR_PID_BO_CONTROLLER

wget https://raw.githubusercontent.com/konstantinos-papadopoulos/PID-BO-CONTROLLER/develop/deploy/requirements.txt -P $BASE_DIR
source $BASE_DIR/$VENV_DIR_PID_BO_CONTROLLER/bin/activate
pip3 $pip_options install -r $BASE_DIR/requirements.txt
rm -rf $BASE_DIR/requirements.txt
deactivate
exit
# This step is required to run jobs with any user
#------------------------------------------------------------------------------
chmod 777 -R $BASE_DIR

