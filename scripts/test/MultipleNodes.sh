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
NNode=2                                                                                           # Nb of Nodes
NProc=4                                                                                           # Nb of Processors
SlncFlg=0                                                                                         # Silence bash file echoes? (0/1 for no/yes)
RmTrajFlg=1                                                                                       # Remove Trajectories Files? (0/1 for no/yes)
#-------------------------------------------------------------------------------------------------#




###################################### DEFINING FUNCTIONS #########################################
function SetData() {
  mkdir -p ${COARSEAIR_OUTPUT_DIR}
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


function LoadComputeRates() {
  source ${COARSEAIR_SH_DIR}/ComputeRates.sh
  if [ ${NNode} -gt 1 ]; then
    ComputeRatesPBS
  else
    ComputeRates
  fi
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


for Ttra in "${Ttra_vec[@]}"; do :
  echo '[CoarseAIR] '
  echo "[CoarseAIR]: ===   Translational Temperature      = " ${Ttra} " ============================ "
  #for Tint in "${Tint_Vec[@]}"; do :
    Tint=${Ttra}
    echo "[CoarseAIR]: ===== Internal Temperature           = " ${Tint} " ======================== "
    echo '[CoarseAIR] '     
    
    
    mkdir -p ${COARSEAIR_OUTPUT_DIR}/"T_"${Ttra%.*}"_"${Tint%.*}
      
      
    # --- Preprocessing Levels ------------------------------------------------------------------ #
    if [ ${PrepFlg} -eq 1 ]; then
      echo "[CoarseAIR]: Calling LoadPreprocLevels"
      LoadPreprocLevels
    fi
    #=============================================================================================#
    

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


if [ ${NNode} -gt 1 ]; then
  echo "[CoarseAIR]: File Executed on Multiple Nodes: waiting for PBS Files being Done!"
  exit 0
fi


# --- Computing Derived Quantities ---------------------------------------------------------- #
if [ ${DerQntFlg} -eq 1 ]; then
  echo "[CoarseAIR]: Calling DeriveQuantities"
  LoadDeriveQuantities
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
