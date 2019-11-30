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

# ----------------------------------------------------------------------------------------------------------------- Clean #
function RmAll {
  
  export COARSEAIR_OUTPUT_DIR=$(pwd)/Test
  export TranFlg=${1}
  export Tran=${2}
  export Tint=${3}
  export MinLevel1=${4}
  export MaxLevel1=${5}
  
  echo "  [Clean]: COARSEAIR_OUTPUT_DIR  = "${COARSEAIR_OUTPUT_DIR}
  echo "  [Clean]: TranFlg            = "${TranFlg}
  echo "  [Clean]: Tran               = "${Tran}
  echo "  [Clean]: Tint               = "${Tint}
  echo "  [Clean]: MinLevel1          = "${MinLevel1}
  echo "  [Clean]: MaxLevel1          = "${MaxLevel1}
  echo "  [Clean]: MinLevel2          = "${MinLevel2}
  echo "  [Clean]: MaxLevel2          = "${MaxLevel2}


  iLevels1=${MinLevel1}
  while [ ${iLevels1} -le ${MaxLevel1} ]
    do
    echo "  [Clean]: --- Molecule 1, Level/Bin " ${iLevels1} " ----------------------------- "

    iLevels2=0  
    #iLevels2=${MinLevel2}
    #while [ ${iLevels2} -le ${MaxLevel2} ]; do
    #  echo "  [Clean]: ----- Molecule 2, Level/Bin = " ${iLevels2} " --------------------- "
    #  echo "  [Clean]"
      
      
      if [ ${TranFlg} -eq 0 ]; then 
        export COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevels1}"_"${iLevels2}
      else
        export COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevels1}"_"${iLevels2}
      fi
      
      
      #if [ -f ${COARSEAIR_BIN_OUTPUT_DIR} ]; then
      #  echo ${COARSEAIR_BIN_OUTPUT_DIR}
        rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}
      #fi
    
    
    #  echo "  [Clean]: ---------------------------------------------------------- "
    #  iLevels2=$((iLevels2+1))
    #done
    
    echo "  [Clean]: ------------------------------------------------------------ "
    echo " "
    iLevels1=$((iLevels1+1))
  done         
      
}
#================================================================================================================================#


# ------------------------------------------------------------------------------------------------------------------------ Clean #
function Clean {
  # Ex: Clean 1 10000.0 10000.0 1 9390 0 0

  export COARSEAIR_OUTPUT_DIR=$(pwd)/Test
  export TranFlg=${1}
  export Tran=${2}
  export Tint=${3}
  export MinLevel1=${4}
  export MaxLevel1=${5}
  
  echo "  [Clean]: COARSEAIR_OUTPUT_DIR  = "${COARSEAIR_OUTPUT_DIR}
  echo "  [Clean]: TranFlg            = "${TranFlg}
  echo "  [Clean]: Tran               = "${Tran}
  echo "  [Clean]: Tint               = "${Tint}
  echo "  [Clean]: MinLevel1          = "${MinLevel1}
  echo "  [Clean]: MaxLevel1          = "${MaxLevel1}
  echo "  [Clean]: MinLevel2          = "${MinLevel2}
  echo "  [Clean]: MaxLevel2          = "${MaxLevel2}


  iLevels1=${MinLevel1}
  while [ ${iLevels1} -le ${MaxLevel1} ]
    do
    echo "  [Clean]: --- Molecule 1, Level/Bin " ${iLevels1} " ----------------------------- "

    iLevels2=0  
    #iLevels2=${MinLevel2}
    #while [ ${iLevels2} -le ${MaxLevel2} ]; do
    #  echo "  [Clean]: ----- Molecule 2, Level/Bin = " ${iLevels2} " --------------------- "
    #  echo "  [Clean]"
      
      
      if [ ${TranFlg} -eq 0 ]; then 
        export COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevels1}"_"${iLevels2}
      else
        export COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevels1}"_"${iLevels2}
      fi
      
      
      if [ -f ${COARSEAIR_BIN_OUTPUT_DIR}/NConvTraj.dat ]; then
        rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/AdjstTrj*
        rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/Node*
        rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/statistics-*
      fi
    
    
    #  echo "  [Clean]: ---------------------------------------------------------- "
    #  iLevels2=$((iLevels2+1))
    #done
    
    echo "  [Clean]: ------------------------------------------------------------ "
    echo " "
    iLevels1=$((iLevels1+1))
  done         
      
}
#================================================================================================================================#
