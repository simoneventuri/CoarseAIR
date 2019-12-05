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

# ---  PARAMETERS ------------------------------------------------------------------------------- #
NNode=1                                                                                           # Nb of Nodes
iNode=1                                                                                           # Current Node
ProcType='ivy'                                                                                    # Only for Clusters; 'san'/'ivy' (Sandy Bridge/Ivy Bridge)
NProc=20                                                                                          # Nb of Processors
SlncFlg=0                                                                                         # =1 -> Silencing Bash File Echoes
MergeAllFlg=0                                                                                     # =1 -> Merging All the ASCI Traj Files in 1 File
RmTrajFlg=1                                                                                       # =1 -> Removing Traj Files from Single Processors
BinaryTrajFlg=0                                                                                   # =1 -> Statistics Reads Binary Traj Files
#-------------------------------------------------------------------------------------------------#


# ---  PATHS ------------------------------------------------------------------------------------ #
export RUN_DIR=$(pwd)
export COARSEAIR_WORKING_DIR=$(pwd)
export COARSEAIR_INPUT_DIR=$(pwd)/input
export COARSEAIR_INPUT_FILE=$(pwd)/input/CoarseAIR.inp
export COARSEAIR_OUTPUT_DIR=$(pwd)/Test
export COARSEAIR_BIN_OUTPUT_DIR=$(pwd)/Test
export $RmTrajFlg
#-------------------------------------------------------------------------------------------------#



###################################### DEFINING FUNCTIONS #########################################
function SetData() {
  mkdir -p ${COARSEAIR_OUTPUT_DIR}
  mkdir -p ${COARSEAIR_OUTPUT_DIR}/"../output"
  [[ -z "${COARSEAIR_BUILD_CONFIG}" ]] && echo "The CoarseAIR modulefile has not been loaded => Stopping" && exit 1
  
  if [ "${SlncFlg}" == "1" ]; then
    SlncStr=' &>/dev/null'
  else
    SlncStr=' '
  fi
  
  COARSEAIR_SH_DIR=${COARSEAIR_SOURCE_DIR}"/scripts/executing"
}


function PrintParameters() {
  echo ' '
  echo '------------------------------------------------------------------------------------------'
  echo ' CoarseAIR: Coarse-Grained Quasi-Classical Trajectories '
  echo '------------------------------------------------------------------------------------------'
  echo ' '
  echo ' Nb of Nodes         = ' $NNode
  echo ' Nb of Processors    = ' $NProc
  echo ' '
  echo '-----------------------'
  echo ' Loading CoarseAIR modulefile'
  echo '-----------------------'
  echo '  COARSEAIR_LIBRARY_NAME   = ' ${COARSEAIR_LIBRARY_NAME}
  echo '  COARSEAIR_VERSION        = ' ${COARSEAIR_VERSION}
  echo '  COARSEAIR_BUILD_CONFIG   = ' ${COARSEAIR_BUILD_CONFIG}
  echo '  COARSEAIR_SOURCE_DIR     = ' ${COARSEAIR_SOURCE_DIR}
  echo '  COARSEAIR_BUILD_DIR      = ' ${COARSEAIR_BUILD_DIR}
  echo '  COARSEAIR_INSTALL_DIR    = ' ${COARSEAIR_INSTALL_DIR}
  echo '  COARSEAIR_EXECUTABLE_DIR = ' ${COARSEAIR_EXECUTABLE_DIR}
  echo '  COARSEAIR_CMAKE_DIR      = ' ${COARSEAIR_CMAKE_DIR}
  echo '  COARSEAIR_LIBRARY_DIR    = ' ${COARSEAIR_LIBRARY_DIR}
  echo '  COARSEAIR_MODULES_DIR    = ' ${COARSEAIR_MODULES_DIR}
  echo '  COARSEAIR_DTB_DIR        = ' ${COARSEAIR_DTB_DIR}
  echo '-----------------------'
  echo ' '
  echo '-----------------------'
  echo ' Paths:'
  echo '-----------------------'
  echo ' Working directory    = ' ${COARSEAIR_WORKING_DIR}
  echo ' Input directory      = ' ${COARSEAIR_INPUT_DIR}
  echo ' CoarseAIR Input file    = ' ${COARSEAIR_INPUT_FILE}
  echo ' CoarseAIR .sh directory = ' ${COARSEAIR_SH_DIR}
  echo ' Output directory     = ' ${COARSEAIR_OUTPUT_DIR}
  echo '-----------------------'
}


function LoadReadBashInput() {
  #################################################################################
  InitializeCommand="coarseair-initialize.x "
  cd ${COARSEAIR_OUTPUT_DIR}
  eval ${InitializeCommand}
  ################################################################################
  source ${COARSEAIR_SH_DIR}/ReadBashInput.sh
  ReadBashInput
}


function LoadComputeRates() {
  source ${COARSEAIR_SH_DIR}/ComputeRates.sh
  PostTrajectoriesFromPBS
}


###################################################################################################
###################################################################################################


#################################### RUNNINQ CoarseAIR CODE #############################################
startTot=`date +%s`

SetData
echo " "


PrintParameters
echo " "


echo "[CoarseAIR]: Calling LoadReadBashInput"
LoadReadBashInput
echo " "


for Tran in "${Tran_vec[@]}"; do :
  echo '[CoarseAIR] '
  if [ ${TranFlg} -eq 0 ]; then 
    echo "[CoarseAIR]: ===   Translational Energy           = " ${Tran} " ============================ "
    #for Tint in "${Tint_Vec[@]}"; do :
    Tint=0.0
  else
    echo "[CoarseAIR]: ===   Translational Temperature      = " ${Tran} " ============================ "
    #for Tint in "${Tint_Vec[@]}"; do :
    Tint=${Tran}
  fi
    echo "[CoarseAIR]: ===== Internal Temperature           = " ${Tint} " ======================== "
    echo '[CoarseAIR] '     
    

    #=============================================================================================#
    if [ ${Mthd1} == "State-Specific" ] || [ ${Mthd1} == "Vib-Specific" ] ; then
      typeset -i NLevels1=$(cat "${COARSEAIR_OUTPUT_DIR}/${System}/${Molecule1}/NLevels.inp")
      echo "[CoarseAIR]: Found "${NLevels1}" Levels for Molecule 1"
      if [ ${VarB1} -le 0 ]; then
        MaxLevel1=${NLevels1}
        echo "[CoarseAIR]: Max Levels for Molecule 1 set to "${MaxLevel1}
      fi
    fi
    if [ ${NMolecules} -gt 1 ]; then
      if [ ${Mthd2} == "State-Specific" ] || [ ${Mthd2} == "Vib-Specific" ] ; then
        typeset -i NLevels2=$(cat ${COARSEAIR_OUTPUT_DIR}/${System}/${Molecule2}/'NLevels.inp')
        echo "[CoarseAIR]: Found "${NLevels2}" Levels for Molecule 2"
        if [ ${VarB2} -le 0 ]; then
          MaxLevel2=${NLevels2}
          echo "[CoarseAIR]: Max Levels for Molecule 2 set to "${MaxLevel2}
        fi
      fi
    fi
    echo " "
    #=============================================================================================#


    # --- Computing Rates ----------------------------------------------------------------------- #
    if [ ${RunTrajFlg} -eq 1 ] || [ ${PostFlg} -eq 1 ]; then
      echo "[CoarseAIR]: Calling LoadComputeRates"
      LoadComputeRates
      echo " "
    fi
    #=============================================================================================#

  
    echo "[CoarseAIR]: ======================================================================== "
  echo "[CoarseAIR]: ============================================================================ "
done
wait


endTot=`date +%s`
Totruntime=$((endTot-startTot))
echo "[CoarseAIR]: Done with CoarseAIR. Total RunTime = "${Totruntime}"s"
echo " "
###################################################################################################
###################################################################################################

exit 0
