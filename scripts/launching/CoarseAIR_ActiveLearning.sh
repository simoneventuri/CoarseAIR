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
ParNodes=0                                                                                        # Nb of Nodes =0 -> Through .sh; >0 -> Through .pbs
ProcType='ivy'                                                                                    # Only for Clusters; 'san'/'ivy' (Sandy Bridge/Ivy Bridge)
Queue='long'                                                                                      # Only for Clusters; 'devel'/'debug'/'low'/'medium'/'long'
WallTime=120                                                                                      # Only for Clusters; WallTime in hours (e.g., 120)
NProc=${1}
SlncFlg=0                                                                                         # =1 -> Silencing Bash File Echoes
MergeAllFlg=0                                                                                     # =1 -> Merging All the ASCI Traj Files in 1 File
RmTrajFlg=0                                                                                       # =1 -> Removing Traj Files from Single Processors
BinaryTrajFlg=0                                                                                   # =1 -> Statistics Reads Binary Traj Files
#-------------------------------------------------------------------------------------------------#


# ---  PATHS ------------------------------------------------------------------------------------ #
export RUN_DIR=***RUN_DIR***
export COARSEAIR_WORKING_DIR=${RUN_DIR}
export COARSEAIR_INPUT_DIR=***COARSEAIR_INPUT_DIR***
export COARSEAIR_INPUT_FILE=${COARSEAIR_INPUT_DIR}/CoarseAIR.inp
export COARSEAIR_OUTPUT_DIR=***COARSEAIR_OUTPUT_DIR***
export COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}
#-------------------------------------------------------------------------------------------------#



###################################### DEFINING FUNCTIONS #########################################
function SetData() {
  mkdir -p ${COARSEAIR_OUTPUT_DIR}
  #mkdir -p ${COARSEAIR_OUTPUT_DIR}/"../output"
  [[ -z "${COARSEAIR_BUILD_CONFIG}" ]] && echo "The CoarseAIR modulefile has not been loaded => Stopping" && exit 1
  
  if [ "${SlncFlg}" == "1" ]; then
    SlncStr=' &>/dev/null'
  else
    SlncStr=' '
  fi
  
  COARSEAIR_SH_DIR=${COARSEAIR_SOURCE_DIR}"/scripts/executing"

  if [ ${ParNodes} -le 1 ]; then
    NNode=1
  else
    NNode=${ParNodes}
  fi
}


function PrintParameters() {
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


function LoadGenerateLevels() {
  source ${COARSEAIR_SH_DIR}/GenerateLevels.sh
  GenerateLevels
}


function LoadPreprocLevels() {
  source ${COARSEAIR_SH_DIR}/PreprocLevels.sh
  PreprocLevels
}


echo

function LoadComputeRates() {
  source ${COARSEAIR_SH_DIR}/ComputeRates.sh

  if [ ${RunTrajFlg} -eq 1 ]; then
    if [ ${ParNodes} -ge 1 ]; then
      echo "[CoarseAIR]: Calling ComputeTrajsPBS"
      echo " "
      ComputeTrajsPBS
      echo " "
      EXIT
    else
      echo "[CoarseAIR]: Calling ComputeTrajs"
      echo " "
      ComputeTrajs
      echo " "
    fi
  fi

  if [ ${PostFlg} -eq 1 ]; then
    if [ ${ParNodes} -ge 1 ]; then
      echo "[CoarseAIR]: Calling PostTrajectoriesPBS"
      echo " "
      PostTrajectoriesPBS
      echo " "
      EXIT
    else
      echo "[CoarseAIR]: Calling PostTrajectoriesAtNode"
      echo " "
      PostTrajectoriesAtNode
      echo " "
    fi
  fi
}


function LoadMergeAllRates() {
  source ${COARSEAIR_SH_DIR}/ComputeRates.sh
  MergeAllRates
}


function LoadDeriveQuantities() {
  source ${COARSEAIR_SH_DIR}/DeriveQuantities.sh
  DeriveQuantities
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


# --- Generating Levels --------------------------------------------------------------------------#
if [ "${GenLevFlg}" == "1" ]; then
  echo "[CoarseAIR]: Calling LoadGenerateLevels"
  LoadGenerateLevels
  echo " "
fi
#=================================================================================================#    


for Tran in "${Tran_vec[@]}"; do :
  echo '[CoarseAIR] '
  if [ ${TranFlg} -eq 0 ]; then 
    echo "[CoarseAIR]: ===   Translational Energy           = " ${Tran} " ============================ "
    #for Tint in "${Tint_Vec[@]}"; do :
    Tint=0.0
  elif [ ${TranFlg} -eq 1 ]; then 
    echo "[CoarseAIR]: ===   Translational Temperature      = " ${Tran} " ============================ "
    #for Tint in "${Tint_Vec[@]}"; do :
    Tint=${Tran}
  elif [ ${TranFlg} -eq 2 ]; then 
    echo "[CoarseAIR]: ===   Mean Translational Temperature = " ${Tran} " ============================ "
    #for Tint in "${Tint_Vec[@]}"; do :
    Tint=${Tran}
  fi
    # echo "[CoarseAIR]: ===== Internal Temperature           = " ${Tint} " ======================== "
    # echo '[CoarseAIR] '     
      
      
    # --- Preprocessing Levels ------------------------------------------------------------------ #
    if [ ${PrepFlg} -eq 1 ]; then
      echo "[CoarseAIR]: Calling LoadPreprocLevels"
      LoadPreprocLevels
    fi
    #=============================================================================================#


    if [ ${TranFlg} -eq 0 ]; then 
      mkdir -p ${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}
    elif [ ${TranFlg} -eq 1 ]; then 
      mkdir -p ${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}
    elif [ ${TranFlg} -eq 2 ]; then 
      mkdir -p ${COARSEAIR_OUTPUT_DIR}/"EMu"
    fi


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
    #=============================================================================================#


    # --- Computing Rates ----------------------------------------------------------------------- #
    if [ ${RunTrajFlg} -eq 1 ] || [ ${PostFlg} -eq 1 ]; then
      echo "[CoarseAIR]: Calling LoadComputeRates"
      LoadComputeRates
    fi
    #=============================================================================================#

  
    echo "[CoarseAIR]: ======================================================================== "
  echo "[CoarseAIR]: ============================================================================ "
done
wait


if [ ${ParNodes} -gt 1 ]; then
  echo "[CoarseAIR]: File Executed on Multiple Nodes: waiting for Trajectories' PBSs being Done!"
  exit 0
fi


# --- Merging All the Trajectories ---------------------------------------------------------- #
if [ ${MergeAllFlg} -eq 1 ]; then
  echo "[CoarseAIR]: Merging All the Trajectories in One File"
  echo " "
  LoadMergeAllRates
  echo " "
fi
#=============================================================================================#


endTot=`date +%s`
Totruntime=$((endTot-startTot))
echo "[CoarseAIR]: Done with CoarseAIR. Total RunTime = "${Totruntime}"s"
echo " "
###################################################################################################
###################################################################################################

exit 0
