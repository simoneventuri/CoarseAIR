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
echo ' CoarseAIR: Coarse-Grained Quasi-Classical Trajectories                                   '
echo '------------------------------------------------------------------------------------------'
echo ' '
echo '------------------------------------------------------------------------------------------'
echo '   PipeLine for Running External Codes                                                    '
echo '------------------------------------------------------------------------------------------'
echo ' '

source ~/.bashrc
#module purge
COARSEAIR_UPDATE
COARSEAIR_release
PLATO_gnu_release

export System='O4_UMN'
export Molecule_vec=('O2')
export FldrName='_VSM'
export Tran_vec=(10000) #(1500 2500 5000 6000 8000 10000 12000 14000 15000 20000)
export T0=1400
export PathToMECVODEFldr=$WORKSPACE_PATH/neqplasma_QCT/ME_CVODE
export PathToDtbFldr=$WORKSPACE_PATH/Air_Database/Run_0D/database/
export PathToRunFldr=$WORKSPACE_PATH/Air_Database/Run_0D/

export DissFlg=0
export InelFlg=1
export ExchFlg1=1
export ExchFlg2=0


ExtCode_SH_DIR=${COARSEAIR_SOURCE_DIR}"/extra/ExtCode_PipeLine/"

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


ExtCode_SH_DIR=${COARSEAIR_SOURCE_DIR}"/extra/ExtCode_PipeLine/"

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

  if [ $DissFlg -eq 5 ]; then
    PLATONORECOMB_gnu_release
    export ExFldr=${PathToMECVODEFldr}/Generic_NoRecomb/
    echo "[RunMECVODE]: Copying MeCvode Executable from: "${ExFldr}/'exec/box_'
    scp ${ExFldr}'/exec/box_' ./
  else
    PLATO_gnu_release
    export ExFldr=${PathToMECVODEFldr}/Generic/
    echo "[RunMECVODE]: Copying MeCvode Executable from: "${ExFldr}/'exec/box_'
    scp ${ExFldr}'/exec/box_' ./
  fi


  if [ $DissFlg -eq 0 ]; then
    export InputFile=${PathToDtbFldr}'/input/'${System}'/NoDiss/T'${TTran}'K.inp'
  elif [ $InelFlg -eq 0 ] && [ $ExchFlg1 -eq 0 ] && [ $ExchFlg2 -eq 0 ]; then
    export InputFile=${PathToDtbFldr}'/input/'${System}'/OnlyDiss/T'${TTran}'K.inp'
  elif [ $DissFlg -eq 5 ]; then
    export InputFile=${PathToDtbFldr}'/input/'${System}'/Maninder/T'${TTran}'K.inp'
  else
    export InputFile=${PathToDtbFldr}'/input/'${System}'/All/T'${TTran}'K.inp'
  fi  
  echo "[RunMECVODE]: Input File: "${InputFile}
  
  echo "[RunMECVODE]: MeCvode will be executed in the Folder: "$(pwd)
  ./box_ ${OPENBLAS_NUM_THREADS} $T0 ${InputFile}
}


for TTran in "${Tran_vec[@]}"; do :
  echo "[RunMECVODE]: Translational Temperature: TTran = "${TTran}

 
  #########################################################################################################################################################
  #########################################################################################################################################################
  echo "Units=cm^3/s" > $PathToDtbFldr"/kinetics/KineticsTEMP_T"$TTran"K_"$System


  #######################################################################################################3
  #### Adding Dissociation Processes
  #### 
  if [ $DissFlg -gt 0 ]; then 

    echo "  [Initialize_0D_Database]: Adding Dissociation Kinetics for "${System}" to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K_"$System
    cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/DissInel.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K_"$System

  fi


  #######################################################################################################3
  #### Adding Inelastic Processes
  #### 
  iInel=${InelFlg}
  if [ ${iInel} -eq 1 ]; then
    echo "  [Initialize_0D_Database]: Adding Inelastic Kinetics for "${System}" to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K_"$System
    cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Inel.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K_"$System
  fi



  #######################################################################################################3
  #### Adding Exchange Processes
  #### 
  iExch=${ExchFlg1}
  if [ ${iExch} -eq 1 ]; then
    echo "  [Initialize_0D_Database]: Adding Exchange Kinetics for "${System}" to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K_"$System
    cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Exch_Type1.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K_"$System
  fi


  #######################################################################################################3
  #### Adding Molecules Info
  #### 
  for Molecule in "${Molecule_vec[@]}"; do :

    echo "  [Initialize_0D_Database]: Copying Thermo File "$PathToDtbFldr/"/thermo/"${Molecule}"_T"${TTran}"K"
    #cat $PathToDtbFldr/"/thermo/"${System}"/"${Molecule}"_T"${TTran}"K" > $PathToDtbFldr/"/thermo/"${Molecule}"_T"${TTran}"K"
    cat $PathToDtbFldr/"/thermo/"${System}${FldrName}"/"${Molecule}"_T"${TTran}"K" > $PathToDtbFldr/"/thermo/"${Molecule}"_T"${TTran}"K"

    echo "  [Initialize_0D_Database]: Copying Initial Mole Fraction File "$PathToDtbFldr/"/thermo/"${Molecule}"_T"${T0}"K"
    #cat $PathToDtbFldr/"/thermo/"${System}"/"${Molecule}"_InitialMoleFracs_T"${T0}"K.dat" > $PathToDtbFldr/"/thermo/"${Molecule}"_InitialMoleFracs_T"${T0}"K.dat"
    cat $PathToDtbFldr/"/thermo/"${System}${FldrName}"/"${Molecule}"_InitialMoleFracs_T"${T0}"K.dat" > $PathToDtbFldr/"/thermo/"${Molecule}"_InitialMoleFracs_T"${T0}"K.dat"

  done
  #########################################################################################################################################################
  #########################################################################################################################################################


  echo "[RunMECCVODE]: Calling Call_MeCvode"
  Call_MeCvode
  echo " "

  rm -rf $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K_"$System

done
