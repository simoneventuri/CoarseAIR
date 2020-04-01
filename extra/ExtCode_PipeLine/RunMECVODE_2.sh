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


# Ex. of Call: bash RunMECVODE.sh O3 10000 $WORKSPACE_PATH/neqplasma_QCT/ME_CVODE $WORKSPACE_PATH/Mars_Database/Run_0D/database/ $WORKSPACE_PATH/Mars_Database/Run_0D/ 1 1 0 -1

echo '------------------------------------------------------------------------------------------'
echo ' CoarseAIR: Coarse-Grained Quasi-Classical Trajectories                     '
echo '------------------------------------------------------------------------------------------'
echo ' '
echo '------------------------------------------------------------------------------------------'
echo '   PipeLine for Running External Codes                            '
echo '------------------------------------------------------------------------------------------'
echo ' '

source ~/.bashrc
module purge
COARSEAIR_UPDATE
COARSEAIR_release
#PLATONORECOMB_gnu_release
PLATO_gnu_release

export System='O3_UMN'
export Molecule='O2'
export FldrName='_OnlyBound'
export Tran_vec=(10000 15000) 
export T0=300
export PathToMECVODEFldr=$WORKSPACE_PATH/neqplasma_QCT/ME_CVODE
export PathToDtbFldr=$WORKSPACE_PATH/Mars_Database/Run_0D/database/
export PathToRunFldr=$WORKSPACE_PATH/Mars_Database/Run_0D/

export DissFlg=2
export InelFlg=1
export ExchFlg1=1
export ExchFlg2=1


ExtCode_SH_DIR=${COARSEAIR_SOURCE_DIR}"/scripts/postprocessing/ExtCode_PipeLine/"

echo '------------------------------------------------------'
echo '  Paths:'
echo '------------------------------------------------------'
echo '  $PLATO_LIB      directory = '${PLATO_LIB}
echo '  ExtCode .sh     directory = '${ExtCode_SH_DIR}
echo '  MeCvode install directory = '${PathToMECVODEFldr}
echo '  MeCvode Dtb     directory = '${PathToDtbFldr}
echo '  MeCvode running directory = '${PathToRunFldr}
echo '------------------------------------------------------'
echo ' '

echo '------------------------------------------------------'
echo '  Inputs:'
echo '------------------------------------------------------'
echo '  System                         = '${System}
echo '  Vector of Translational Temp.s = '${Tran_vec}
echo '  Writing Dissociation?          = '${DissFlg}
echo '  Writing Inelastic?             = '${InelFlg}
echo '  Firts Exchanges to be Written? = '${ExchFlg1}
echo '  Last Exchanges  to be Written? = '${ExchFlg2}
echo '------------------------------------------------------'
echo ' '


function Load_Initialize_0D() {
  source ${ExtCode_SH_DIR}/Initialize_0D_Database_Function.sh
  Initialize_0D_Database
}


function Call_MeCvode() {
  cd ${PathToRunFldr}
  export OutputFldr='output_'${System}${FldrName}'_T'${TTran}'K_'${DissFlg}'_'${InelFlg}'_'${ExchFlg1}'_'${ExchFlg2}
  mkdir -p ./${OutputFldr}
  cd ./${OutputFldr} 
  if [ $DissFlg -eq 0 ]; then
    export ExFldr=${PathToMECVODEFldr}/${System}/'Mars_T'${TTran}'K_Danil_NoDiss'
  else
    export ExFldr=${PathToMECVODEFldr}/${System}/'Mars_T'${TTran}'K_Danil'
  fi  
  echo "[RunMECVODE]: Copying MeCvode Executable from "${ExFldr}/'exec/box_'
  scp ${ExFldr}'/exec/box_' ./
  echo "[RunMECVODE]: MeCvode will be executed in the Folder "$(pwd)
  ./box_ 8
}


for TTran in "${Tran_vec[@]}"; do :
  echo "[RunMECVODE]: Translational Temperature, TTran = "${TTran}

  echo "[RunMECVODE]: Calling Load_Initialize_0D"
  Load_Initialize_0D
  echo " "

  echo "[RunMECCVODE]: Calling Call_MeCvode"
  Call_MeCvode
  echo " "

done