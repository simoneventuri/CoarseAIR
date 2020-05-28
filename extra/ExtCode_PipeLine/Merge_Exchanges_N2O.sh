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

export System1='N2O_UMN'
export System2='NON_UMN'
export Molecule_vec=(N2 NO)
export FldrName=''
export Tran_vec=(5000 10000 20000) #(2500 5000 7500 10000 12500 15000 20000) 
export T0=300
export PathToMECVODEFldr=$WORKSPACE_PATH/neqplasma_QCT/ME_CVODE
export PathToDtbFldr=$WORKSPACE_PATH/Mars_Database/Run_0D/database/
export PathToRunFldr=$WORKSPACE_PATH/Mars_Database/Run_0D/

export DissFlg=0
export InelFlg=1
export ExchFlg1=1
export ExchFlg2=0

export NBins=0

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



for TTran in "${Tran_vec[@]}"; do :
  echo "[MergeExchanges]: Translational Temperature, TTran = "${TTran}


  echo "[MergeExchanges]:   Copying Exchange of N2+O"
  scp $PathToDtbFldr"/kinetics/"${System1}${FldrName}"/T"$TTran"K/Exch_Type1.dat"    $PathToDtbFldr"/kinetics/"${System1}${FldrName}"/T"$TTran"K/Exch_Comb_Type1.dat"
  scp $PathToDtbFldr"/kinetics/"${System1}${FldrName}"/T"$TTran"K/Exch_Type1.dat"    $PathToDtbFldr"/kinetics/"${System2}${FldrName}"/T"$TTran"K/Exch_Comb_Type1.dat"


  echo "[MergeExchanges]:   Copying Exchange of NO+N"
  cat $PathToDtbFldr"/kinetics/"${System2}${FldrName}"/T"$TTran"K/Exch_Type1.dat" >> $PathToDtbFldr"/kinetics/"${System1}${FldrName}"/T"$TTran"K/Exch_Comb_Type1.dat"
  cat $PathToDtbFldr"/kinetics/"${System2}${FldrName}"/T"$TTran"K/Exch_Type1.dat" >> $PathToDtbFldr"/kinetics/"${System2}${FldrName}"/T"$TTran"K/Exch_Comb_Type1.dat"

  cat $PathToDtbFldr"/kinetics/"${System2}${FldrName}"/T"$TTran"K/Exch_Type2.dat" >> $PathToDtbFldr"/kinetics/"${System2}${FldrName}"/T"$TTran"K/Exch_Comb_Type2.dat"


done