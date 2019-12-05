#!/bin/bash
#===============================================================================================================
# 
# Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
# 
# Copyright (C) 2018 Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
#
# Based on "VVTC" (Vectorized Variable stepsize Trajectory Code) by David Schwenke (NASA Ames Research Center). 
# 
# This program is free software; you can redistribute it and/or modify it under the terms of the 
# Version 2.1 GNU Lesser General Public License as published by the Free Software Foundation. 
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU Lesser General Public License for more details. 
# 
# You should have received a copy of the GNU Lesser General Public License along with this library; 
# if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
# 
#---------------------------------------------------------------------------------------------------------------
#===============================================================================================================

# ${1}   =   COARSEAIR_OUTPUT_DIR     =   Path to the output directory
# ${2}   =   COARSEAIR_WORKING_DIR    =   Path to the working directory
# ${3}   =   COARSEAIR_BIN_OUTPUT_DIR =   Path to the bin output directory
# ${4}   =   NNode                 =   Total Nb of Nodes
# ${5}   =   iNode                 =   Current Node
# ${6}   =   NProc                 =   Nb of Processors per Node
# ${7}   =   iProc                 =   Current Processor


COARSEAIR_WORKING_DIR_TEMP=${3}
cd ${COARSEAIR_WORKING_DIR_TEMP}
#. ~/.bash_profile
#COARSEAIR_UPDATE
#COARSEAIR_release

COARSEAIR_OUTPUT_DIR=${1}
COARSEAIR_BIN_OUTPUT_DIR=${2}
COARSEAIR_WORKING_DIR=${COARSEAIR_WORKING_DIR_TEMP}
NNode=${4}
iNode=${5}
NProc=${6}
iProc=${7}

#export RUN_DIR=${COARSEAIR_WORKING_DIR}
#export COARSEAIR_INPUT_DIR=${RUN_DIR}/input
#export COARSEAIR_INPUT_FILE=${RUN_DIR}/input/CoarseAIR.inp
#export COARSEAIR_DTB_DIR=${RUN_DIR}/dtb


echo "      [ComputeLevelsAtProc.sh]: COARSEAIR_OUTPUT_DIR     = "${COARSEAIR_OUTPUT_DIR}
echo "      [ComputeLevelsAtProc.sh]: COARSEAIR_BIN_OUTPUT_DIR = "${COARSEAIR_BIN_OUTPUT_DIR}
echo "      [ComputeLevelsAtProc.sh]: NNode                    = "${NNode}
echo "      [ComputeLevelsAtProc.sh]: iNode                    = "${iNode}
echo "      [ComputeLevelsAtProc.sh]: NProc                    = "${NProc}
echo "      [ComputeLevelsAtProc.sh]: iProc                    = "${iProc}


ComputeLevelsCommand="coarseair-computelevels.x"

cd ${COARSEAIR_OUTPUT_DIR}

echo "      [ComputeLevelsAtProc.sh]: Generating Energy Levels for Node "${iNode}", Processor "${iProc}
eval ${ComputeLevelsCommand} ${COARSEAIR_BIN_OUTPUT_DIR} ${NProc} ${iProc} ${NNode} ${iNode}
echo "      [ComputeLevelsAtProc.sh]: Done Generating Energy Levels for Node "${iNode}", Processor "${iProc}

exit 0
