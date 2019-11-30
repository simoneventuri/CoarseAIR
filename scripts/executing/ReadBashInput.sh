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

# ---------------------------------------------------------------------------------------------------------------- ReadBashInput #
function ReadBashInput {

  
  
  #echo "  [ReadBashInput]: COARSEAIR_OUTPUT_DIR     = "${COARSEAIR_OUTPUT_DIR}
  
  ######################################################################################
  #InitializeCommand="coarseair-initialize.x "${SilencingString}
  #
  #cd ${COARSEAIR_OUTPUT_DIR}
  #eval ${InitializeCommand}
  ######################################################################################
  
  echo "  [ReadBashInput]: "
  echo "  [ReadBashInput]:  -----------------------"
  echo "  [ReadBashInput]:  Input Variables: "
  echo "  [ReadBashInput]:  -----------------------"
  iLine=1
  System=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  System=$(echo $System | tr -d " ")
  echo "  [ReadBashInput]:  System = "${System}
  echo "  [ReadBashInput]:  "
  
  iLine=$((iLine+1))
  TranFlg=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  TranFlg=$(echo $TranFlg | tr -d " ")
  echo "  [ReadBashInput]:  TranFlg = "${TranFlg}
  
  
  if [ ${TranFlg} -eq 0 ]; then 
    echo "  [ReadBashInput]:  Translational Energy Specific. CoarseAIR will compute Cross Sections!"
    
    iLine=$((iLine+1))
    NTran=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
    NTran=$(echo $NTran | tr -d " ")
    echo "  [ReadBashInput]:  NTran = "${NTran}
    iEtra=1
    while [ ${iEtra} -le ${NTran} ]; do
      iLine=$((iLine+1))
      xTemp=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      xTemp=$(echo $xTemp | tr -d " ")
      Tran_vec[${iEtra}]=${xTemp}
      iEtra=$((iEtra+1))
    done
    echo "  [ReadBashInput]:  Tran_vec = "${Tran_vec[@]}
    
  else
    echo "  [ReadBashInput]:  Translational Temperature Specific. CoarseAIR will compute Rates!"
    
    iLine=$((iLine+1))
    NTran=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
    NTran=$(echo $NTran | tr -d " ")
    echo "  [ReadBashInput]:  NTran = "${NTran}
    iTtra=1
    while [ ${iTtra} -le ${NTran} ]; do
      iLine=$((iLine+1))
      xTemp=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      xTemp=$(echo $xTemp | tr -d " ")
      Tran_vec[${iTtra}]=${xTemp}
      iTtra=$((iTtra+1))
    done
    echo "  [ReadBashInput]:  Tran_vec = "${Tran_vec[@]}
  fi
  

  iLine=$((iLine+1))
  NTint=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  NTint=$(echo $NTint | tr -d " ")
  echo "  [ReadBashInput]:  NTint = "${NTint}
  if [[ ${NTint} -gt 1 ]]; then
    iTint=1
    while [ ${iTint} -le ${NTint} ]; do
      iLine=$((iLine+1))
      xTemp=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      xTemp=$(echo ${xTemp} | tr -d " ")
      Tint_vec[${iTint}]=${xTemp}
      iTint=$((iTint+1))
    done
    Tint_vec
    echo "  [ReadBashInput]:  Tint_vec = "${Tint_vec[@]}
  fi
  echo "  [ReadBashInput]:  "
  

  iLine=$((iLine+1))
  NMolecules=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  NMolecules=$(echo $NMolecules | tr -d " ")
  echo "  [ReadBashInput]:  NMolecules = "${NMolecules}
  iMolecules=1
  while [ ${iMolecules} -le ${NMolecules} ]; do
    if [[ ${iMolecules} -eq 1 ]]; then
      iLine=$((iLine+1))
      Molecule1=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      Molecule1=$(echo $Molecule1 | tr -d " ")
      echo "  [ReadBashInput]:    Molecule1 = "${Molecule1}
      iLine=$((iLine+1))
      Mthd1=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      Mthd1=$(echo $Mthd1 | tr -d " ")
      echo "  [ReadBashInput]:    Mthd1     = "${Mthd1}
      iLine=$((iLine+1))
      VarA1=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      VarA1=$(echo $VarA1 | tr -d " ")
      echo "  [ReadBashInput]:    VarA1     = "${VarA1}
      iLine=$((iLine+1))
      VarB1=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      VarB1=$(echo $VarB1 | tr -d " ")
      echo "  [ReadBashInput]:    VarB1     = "${VarB1}  
      iLine=$((iLine+1))
      VarC1=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      VarC1=$(echo $VarC1 | tr -d " ")
      echo "  [ReadBashInput]:    VarC1     = "${VarC1}   
    elif [[ ${iMolecules} -eq 2 ]]; then
      iLine=$((iLine+1))
      Molecule2=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      Molecule2=$(echo $Molecule2 | tr -d " ")
      echo "  [ReadBashInput]:    Molecule2 = "${Molecule2}
      iLine=$((iLine+1))
      Mthd2=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      Mthd2=$(echo $Mthd2 | tr -d " ")
      echo "  [ReadBashInput]:    Mthd2     = "${Mthd2}
      iLine=$((iLine+1))
      VarA2=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      VarA2=$(echo $VarA2 | tr -d " ")
      echo "  [ReadBashInput]:    VarA2     = "${VarA2}
      iLine=$((iLine+1))
      VarB2=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      VarB2=$(echo $VarB2 | tr -d " ")
      echo "  [ReadBashInput]:    VarB2     = "${VarB2}  
      iLine=$((iLine+1))
      VarC2=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
      VarC2=$(echo $VarC2 | tr -d " ")
      echo "  [ReadBashInput]:    VarC2     = "${VarC2}        
    fi  
    iMolecules=$((iMolecules+1))  
  done  
  
  if [ ${Mthd1} == "State-Specific" ] || [ ${Mthd1} == "Vib-Specific" ] ; then
    if [ ${VarA1} -lt 0 ]; then
      MinLevel1=1
    else
      MinLevel1=${VarA1}
    fi
    if [ ${VarB1} -gt 0 ]; then
      MaxLevel1=${VarB1}
    fi
    echo "  [ReadBashInput]:    MinLevel1 = "${MinLevel1}
    echo "  [ReadBashInput]:    MaxLevel1 = "${MaxLevel1}
    #echo "  [ReadBashInput]:    NLevels1  = "${NLevels1}
  elif [ ${Mthd1} == "From-File" ] || [ ${Mthd1} == "RoVib-CG" ]; then
    NLevels1=${VarA1}
    MinLevel1=${VarB1}
    MaxLevel1=${VarC1}
  fi
  if [ ${NMolecules} -gt 1 ]; then
    if [ ${Mthd2} == "State-Specific" ] || [ ${Mthd2} == "Vib-Specific" ] ; then
      if [ ${VarA2} -lt 0 ]; then
        MinLevel2=1
      else
        MinLevel2=${VarA2}
      fi
      if [ ${VarB2} -gt 0 ]; then
        MaxLevel2=${VarB2}
      fi
    elif [ ${Mthd1} == "From-File" ] || [ ${Mthd2} == "RoVib-CG" ]; then
      NLevels2=${VarA2}
      MinLevel2=${VarB2}
      MaxLevel2=${VarC2}
      echo "  [ReadBashInput]:    MinLevel2 = "${MinLevel2}
      echo "  [ReadBashInput]:    MaxLevel2 = "${MaxLevel2}
      #echo "  [ReadBashInput]:    NLevels2  = "${NLevels2}
    fi
  else
    VarA2=0
    VarB2=0
    VarC2=0
    Molecule2=None
    MinLevel2=0
    MaxLevel2=0
    NLevels2=0
  fi
  echo "  [ReadBashInput]:  "

  iLine=$((iLine+1))
  GenLevFlg=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  GenLevFlg=$(echo $GenLevFlg | tr -d " ")
  echo "  [ReadBashInput]:  Generatig Energy Levels?           = "${GenLevFlg}
  
  iLine=$((iLine+1))
  PrepFlg=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  PrepFlg=$(echo $PrepFlg | tr -d " ")
  echo "  [ReadBashInput]:  Preprocessing Levels?              = "${PrepFlg}
  
  iLine=$((iLine+1))
  RunTrajFlg=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  RunTrajFlg=$(echo $RunTrajFlg | tr -d " ")
  echo "  [ReadBashInput]:  Running Trajectories?              = "${RunTrajFlg}
  
  iLine=$((iLine+1))
  PostFlg=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  PostFlg=$(echo $PostFlg | tr -d " ")
  echo "  [ReadBashInput]:  Postprocessing Trajectories?       = "${PostFlg}
  
  iLine=$((iLine+1))
  DerQntFlg=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  DerQntFlg=$(echo $DerQntFlg | tr -d " ")
  echo "  [ReadBashInput]:  Computing Derived Quantities?      = "${DerQntFlg}

  iLine=$((iLine+1))
  StochPESFlg=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  StochPESFlg=$(echo $StochPESFlg | tr -d " ")
  echo "  [ReadBashInput]:  Stochastic PES?                    = "${StochPESFlg}
  
  iLine=$((iLine+1))
  NPESs=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  NPESs=$(echo $NPESs | tr -d " ")
  echo "  [ReadBashInput]:  Nb of PESs                         = "${NPESs}
  
  iLine=$((iLine+1))
  iPESStart=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  iPESStart=$(echo $iPESStart | tr -d " ")
  iPESStart=$((iPESStart+1))
  echo "  [ReadBashInput]:  Starting Index for Stochastic PESs = "${iPESStart}
  echo "  [ReadBashInput]:  "

  iLine=$((iLine+1))
  RunExtCodeFlg=$(sed -n ${iLine}'p' ${COARSEAIR_OUTPUT_DIR}/'InputForBash.inp')
  RunExtCodeFlg=$(echo $RunExtCodeFlg | tr -d " ")
  echo "  [ReadBashInput]:  RunExtCodeFlg = "${RunExtCodeFlg}
  
  echo "  [ReadBashInput]:  -----------------------"
    
}
#================================================================================================================================#
