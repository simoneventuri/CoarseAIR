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

COARSEAIR_WORKING_DIR_TEMP=${3}
cd ${COARSEAIR_WORKING_DIR_TEMP}
#. ~/.bash_profile
#COARSEAIR_UPDATE
#COARSEAIR_release

System=${1}
COARSEAIR_OUTPUT_DIR=${2}
COARSEAIR_WORKING_DIR=${COARSEAIR_WORKING_DIR_TEMP}
NNode=${4}
iNode=${5}
NProc=${6}
iProc=${7}
TranFlg=${8}
Tran=${9}                                                                                             
Tint=${10}
Tint=${Tran}
iLevels1=${11}
iLevels2=${12}


#export RUN_DIR=${COARSEAIR_WORKING_DIR}
#export COARSEAIR_INPUT_DIR=${RUN_DIR}/input
#export COARSEAIR_INPUT_FILE=${RUN_DIR}/input/CoarseAIR.inp
#export COARSEAIR_DTB_DIR=${RUN_DIR}/dtb


NProcTot=${NProc}

echo "      [RunTrajectoriesAtProc.sh]: COARSEAIR_OUTPUT_DIR = "${COARSEAIR_OUTPUT_DIR}
echo "      [RunTrajectoriesAtProc.sh]: iNode                = "${iNode}
echo "      [RunTrajectoriesAtProc.sh]: NProcTot             = "${NProcTot}
echo "      [RunTrajectoriesAtProc.sh]: iProc                = "${iProc}
echo "      [RunTrajectoriesAtProc.sh]: TranFlg              = "${TranFlg}
echo "      [RunTrajectoriesAtProc.sh]: Tran                 = "${Tran}
echo "      [RunTrajectoriesAtProc.sh]: Tint                 = "${Tint}
echo "      [RunTrajectoriesAtProc.sh]: iLevels1             = "${iLevels1}
echo "      [RunTrajectoriesAtProc.sh]: iLevels2             = "${iLevels2}

RunTrajectoriesCommand="coarseair-runtrajectories.x"

if [ ${TranFlg} -eq 0 ]; then 
  COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevels1}"_"${iLevels2}
elif [ ${TranFlg} -eq 1 ]; then 
  COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevels1}"_"${iLevels2}
elif [ ${TranFlg} -eq 2 ]; then 
  COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"EMu/Bins_"${iLevels1}"_"${iLevels2}
else 
  echo "      [RunTrajectoriesAtProc.sh]: ERROR! Wrong Model for the Tranlational Energy! "
  stop
fi
echo "      [RunTrajectoriesAtProc.sh]: COARSEAIR_BIN_OUTPUT_DIR = "${COARSEAIR_BIN_OUTPUT_DIR}

mkdir -p ${COARSEAIR_BIN_OUTPUT_DIR}/Node_${iNode}/Proc_${iProc}
cd ${COARSEAIR_BIN_OUTPUT_DIR}/Node_${iNode}/Proc_${iProc}
scp ../levels* ./

echo "      [RunTrajectoriesAtProc.sh]: Running Trajectories for Node "${iNode}" of "${NNode}", Processor "${iProc}
eval ${RunTrajectoriesCommand} '${COARSEAIR_BIN_OUTPUT_DIR}' ${Tran} ${Tint} ${NProcTot} ${iProc} ${NNode} ${iNode} ${iLevels1} ${iLevels2}
echo "      [RunTrajectoriesAtProc.sh]: Done Running Trajectories for Node "${iNode}", Processor "${iProc}

exit 0