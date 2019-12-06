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
set -e


#. ~/.bash_profile
#COARSEAIR_UPDATE
#COARSEAIR_release


COARSEAIR_WORKING_DIR=${1}
COARSEAIR_OUTPUT_DIR=${2}
COARSEAIR_SH_DIR=${3}
System=${4}
StochPESFlg=${5}
NPESs=${6}
iPESStart=${7}
TranFlg=${8}
Tran=${9}
Tint=${10}
Velocity=${11}
NNode=${12}
iNode=${13}
NProc=${14}
iProc=${15}
Molecule1=${16}
NLevels1=${17}
MinLevel1=${18}
MaxLevel1=${19}
Molecule2=${20}
NLevels2=${21}
MinLevel2=${22}
MaxLevel2=${23}
RmTrajFlg=${24}
BinaryTrajFlg=${25}                                                            



echo "    [PostTrajectoriesAtProc.sh]: COARSEAIR_WORKING_DIR = "${COARSEAIR_WORKING_DIR}
echo "    [PostTrajectoriesAtProc.sh]: COARSEAIR_OUTPUT_DIR  = "${COARSEAIR_OUTPUT_DIR}
echo "    [PostTrajectoriesAtProc.sh]: COARSEAIR_SH_DIR      = "${COARSEAIR_SH_DIR}
echo "    [PostTrajectoriesAtProc.sh]: System                = "${System}
echo "    [PostTrajectoriesAtProc.sh]: StochPESFlg           = "${StochPESFlg}
echo "    [PostTrajectoriesAtProc.sh]: NPESs                 = "${NPESs}
echo "    [PostTrajectoriesAtProc.sh]: iPESStart             = "${iPESStart}
echo "    [PostTrajectoriesAtProc.sh]: TranFlg               = "${TranFlg}
echo "    [PostTrajectoriesAtProc.sh]: Tran                  = "${Tran}
echo "    [PostTrajectoriesAtProc.sh]: Tint                  = "${Tint}
echo "    [PostTrajectoriesAtProc.sh]: Velocity              = "${Velocity}
echo "    [PostTrajectoriesAtProc.sh]: NNode                 = "${NNode}
echo "    [PostTrajectoriesAtProc.sh]: iNode                 = "${iNode}
echo "    [PostTrajectoriesAtProc.sh]: NProc                 = "${NProc}
echo "    [PostTrajectoriesAtProc.sh]: iProc                 = "${iProc}
echo "    [PostTrajectoriesAtProc.sh]: Molecule1             = "${Molecule1}
echo "    [PostTrajectoriesAtProc.sh]: NLevels1              = "${NLevels1}
echo "    [PostTrajectoriesAtProc.sh]: MinLevel1             = "${MinLevel1}
echo "    [PostTrajectoriesAtProc.sh]: MaxLevel1             = "${MaxLevel1}
echo "    [PostTrajectoriesAtProc.sh]: Molecule2             = "${Molecule2}
echo "    [PostTrajectoriesAtProc.sh]: NLevels2              = "${NLevels2}
echo "    [PostTrajectoriesAtProc.sh]: MinLevel2             = "${MinLevel2}
echo "    [PostTrajectoriesAtProc.sh]: MaxLevel2             = "${MaxLevel2}
echo "    [PostTrajectoriesAtProc.sh]: RmTrajFlg             = "${RmTrajFlg}
echo "    [PostTrajectoriesAtProc.sh]: BinaryTrajFlg         = "${BinaryTrajFlg}

source ${COARSEAIR_SH_DIR}/ComputeRates.sh


if [ ${MinLevel1} -eq 0 ]; then 
  echo "    [PostTrajectoriesAtProc.sh]: Reading Levels/Bins from File "${COARSEAIR_INPUT_DIR}/LevelsToRunList.inp
  
  if [ ${NNode} -eq 1 ]; then
    iNode=1
  else
    iNode=${MaxLevel1}
  fi
  NLevelsToRun=$(wc -l < "${COARSEAIR_INPUT_DIR}/LevelsToRunList.inp")
  NLevelsPerNode="$(bc <<< "scale = 10; ${NLevelsToRun} / ${NNode}")"
  NLevelsPerNode="$(echo ${NLevelsPerNode} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
  MinLevelInNode=$(($((iNode-1))*NLevelsPerNode+1))
  MaxLevelInNode=$((iNode*NLevelsPerNode))
  #echo "    [PostTrajectoriesAtProc.sh]: For Node "${iNode}", the first Level to read from file is the "${MinLevelInNode}"-th in the List"
  #echo "    [PostTrajectoriesAtProc.sh]: For Node "${iNode}", the last  Level to read from file is the "${MaxLevelInNode}"-th in the List"

  NLevelsPerProc="$(bc <<< "scale = 10; ${NLevelsPerNode} / ${NProc}")"
  NLevelsPerProc="$(echo ${NLevelsPerProc} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
  MinLevelInProc=$(($(($((iProc-1))*NLevelsPerProc+1))+${MinLevelInNode}-1))
  MaxLevelInProc=$((iProc*NLevelsPerProc+${MinLevelInNode}-1))
  if [ ${MaxLevelInProc} -gt ${MaxLevelInNode} ]; then
    MaxLevelInProc=${MaxLevelInNode}
  fi
  echo "    [PostTrajectoriesAtProc.sh]: For Node "${iNode}", Proc "${iProc}", the first Level to read from file is the "${MinLevelInProc}"-th in the List"
  echo "    [PostTrajectoriesAtProc.sh]: For Node "${iNode}", Proc "${iProc}", the last  Level to read from file is the "${MaxLevelInProc}"-th in the List"
  
  iCount=0
  while IFS= read -r iLevel1
    do
    iCount=$((iCount+1))
    
    if [ ${iCount} -ge ${MinLevelInProc} ] && [ ${iCount} -le ${MaxLevelInNode} ]; then
    
      echo "    [PostTrajectoriesAtProc.sh]: --- Molecule 1, Level/Bin " ${iLevel1} " ----------------------------- "

      #iLevel2=0  
      iLevel2=${MinLevel2}
      while [ ${iLevel2} -le ${MaxLevel2} ]; do
        echo "    [PostTrajectoriesAtProc.sh]: ----- Molecule 2, Level/Bin = " ${iLevel2} " --------------------- "        
        
        if [ ${TranFlg} -eq 0 ]; then 
          COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
        else
          COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
        fi
        

        VelocityFile=${COARSEAIR_OUTPUT_DIR}/Velocity_${Tran}.dat
        if [ -f $exist ]; then
          Velocity=$(sed '2q;d' ${VelocityFile})
          echo "    [PostTrajectoriesAtProc.sh]: Velocity = "${Velocity}
        else
          echo "    [PostTrajectoriesAtProc.sh]: ERROR! Velocity File does exist! CoarseAIR cannot compute Cross Sections!"
          exit1
        fi
        
        echo "    [PostTrajectoriesAtProc.sh]: Calling PostTrajectories"
        PostTrajectories


        if [ ${RmTrajFlg} -eq 1 ] && [ -f ${COARSEAIR_BIN_OUTPUT_DIR}/NConvTraj.dat ]; then
          rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/NConvTraj.dat
          rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/Node*
          rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/statistics*
          rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/*.log
          if [ ${BinaryTrajFlg} -eq 1 ]; then
            rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv*
          fi
        fi

    
        echo "    [PostTrajectoriesAtProc.sh]: ----- Molecule 2, Level/Bin = " ${iLevel2} " ------------------- DONE -- "
        iLevel2=$((iLevel2+1))
      done
    
    echo "    [PostTrajectoriesAtProc.sh]: --- Molecule 1, Level/Bin " ${iLevel1} " --------------------- DONE -- "
    fi
    

  done < "${COARSEAIR_INPUT_DIR}/LevelsToRunList.inp"

else

  iLevel1=${MinLevel1}
  while [ ${iLevel1} -le ${MaxLevel1} ]
    do
    echo "    [PostTrajectoriesAtProc.sh]: --- Molecule 1, Level/Bin  = " ${iLevel1} " --------------------------- "

    #iLevel2=0  
    iLevel2=${MinLevel2}
    while [ ${iLevel2} -le ${MaxLevel2} ]; do
      echo "    [PostTrajectoriesAtProc.sh]: ----- Molecule 2, Level/Bin = " ${iLevel2} " --------------------- "
      
      
      if [ ${TranFlg} -eq 0 ]; then 
        COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
      else
        COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
      fi
      

      VelocityFile=${COARSEAIR_OUTPUT_DIR}/Velocity_${Tran}.dat
      if [ -f $exist ]; then
        Velocity=$(sed '2q;d' ${VelocityFile})
        echo "    [PostTrajectoriesAtProc.sh]: Velocity = "${Velocity}
      else
        echo "    [PostTrajectoriesAtProc.sh]: ERROR! Velocity File does exist! CoarseAIR cannot compute Cross Sections!"
        exit1
      fi
      
      echo "    [PostTrajectoriesAtProc.sh]: Calling PostTrajectories"
      PostTrajectories
      
      if [ ${RmTrajFlg} -eq 1 ] && [ -f ${COARSEAIR_BIN_OUTPUT_DIR}/NConvTraj.dat ]; then
        rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/NConvTraj.dat
        rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/Node*
        rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/statistics*
        rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/*.log
        if [ ${BinaryTrajFlg} -eq 1 ]; then
          rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv*
        fi
      fi


      echo "    [PostTrajectoriesAtProc.sh]: ----- Molecule 2, Level/Bin = " ${iLevel2} " ------------------- DONE -- "
      iLevel2=$((iLevel2+1))
    done
    
    echo "    [PostTrajectoriesAtProc.sh]: --- Molecule 1, Level/Bin = " ${iLevel1} " ------------------- DONE -- "

    iLevel1=$((iLevel1+1))
  done         

fi


exit 0