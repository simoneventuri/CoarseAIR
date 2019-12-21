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


# Ex. of Call: bash RunKONIG.sh O3 10000 $WORKSPACE_PATH/neqplasma_QCT/cvode_code/install $WORKSPACE_PATH/Mars_Database/Run_0D/database/kinetics/ $WORKSPACE_PATH/Mars_Database/Run_0D/ 1 1 1 1

echo '------------------------------------------------------------------------------------------'
echo ' CoarseAIR: Coarse-Grained Quasi-Classical Trajectories 								    '
echo '------------------------------------------------------------------------------------------'
echo ' '
echo '------------------------------------------------------------------------------------------'
echo '   PipeLine for Running External Codes				 								    '
echo '------------------------------------------------------------------------------------------'
echo ' '

export System=${1}
export TTran=${2}
export PathToKONIGFldr={3}
export PathToKinFldr=${4}
export PathToRunFldr=${5}
export DissFlg=${6}
export InelFlg=${7}
export ExchFlg1=${8}
export ExchFlg2=${9}

ExtCode_SH_DIR=${COARSEAIR_SOURCE_DIR}"/scripts/postprocessing/ExtCode_PipeLine/"

echo '------------------------------------------------------'
echo '  Paths:'
echo '------------------------------------------------------'
echo '  ExtCode .sh  directory   = '${ExtCode_SH_DIR}
echo '  KONING install directory = '${PathToKONIGFldr}
echo '  Kinetic Data directory   = '${PathToKinFldr}
echo '  KONING running directory = '${PathToRunFldr}
echo '------------------------------------------------------'
echo ' '

echo '------------------------------------------------------'
echo '  Inputs:'
echo '------------------------------------------------------'
echo '  System                         = '${System}
echo '  Translational Temp.            = '${TTran}
echo '  Writing Dissociation?          = '${DissFlg}
echo '  Writing Inelastic?             = '${InelFlg}
echo '  Firts Exchanges to be Written? = '${ExchFlg1}
echo '  Last Exchanges  to be Written? = '${ExchFlg2}
echo '------------------------------------------------------'
echo ' '


function Load_Initialize_0D{
  source $ExtCode_SH_DIR/Initialize_0D_Database_Function.sh
  Initialize_0D_Database_Function
}


function Call_KONIG{
  cd ${PathToRunFldr}
  export OutputFldr='output_'${System}'_T'${TTran}'K_'${DissFlg}'_'${InelFlg}'_'${ExchFlg1}'_'${ExchFlg2}
  echo "[RunKONIG]: KONIG will be executed in the Folder "${OutputFldr}
  mkdir -p ./${OutputFldr}
  cd ./{OutputFldr}
  scp ${PathToKONIGFldr}/${System}/'T'${TTran}'K/run_' ./
  ./run_ ${OMP_PROC}
}



echo "[RunKONIG]: Calling Load_Initialize_0D"
Load_Initialize_0D
echo " "


echo "[RunKONIG]: Calling Call_KONIG"
Call_KONIG
echo " "