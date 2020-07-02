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
#set -e


# . ~/.bash_profile
# COARSEAIR_UPDATE
# COARSEAIR_release




########################################################################################################################################################################## <= SplitTrajsPESsHERE
## Splitting Lines in Trajectory File Based on the Collision PES
##
function SplitTrajsPESsHERE
{

  iPES=${iPESStart}
  while [[ ${iPES} -le ${iPESEnd} ]]; do
    PESFile=${COARSEAIR_BIN_OUTPUT_DIR}'/trajectories.csv.'${iPES}
    echo '#    iTraj, iPES,       bmax,        b_i,           j_i,           v_i,         arr_i,           j_f,           v_f,         arr_f' > ${PESFile} 
    iPES=$(($iPES+1))
  done        

  OrigFile=$COARSEAIR_BIN_OUTPUT_DIR/trajectories.csv
  OLDIFS=$IFS
  IFS=','
  [ ! -f ${OrigFile} ] && { echo "$OrigFile file not found"; exit 99; }
  sed 1d ${OrigFile} | while read Idx iPES bmax b_i j_i v_i arr_i j_f v_f arr_f
  do
    jPES=$(($iPES+0))
    if [[ ${jPES} -ge ${iPESStart} ]] && [[ ${jPES} -le ${iPESEnd} ]]; then
      PESFile=${COARSEAIR_BIN_OUTPUT_DIR}'/trajectories.csv.'${jPES}
      echo $Idx','$iPES','$bmax','$b_i','$j_i','$v_i','$arr_i','$j_f','$v_f','$arr_f >> ${PESFile} 
    fi
  done 
  IFS=$OLDIFS

}
########################################################################################################################################################################## <= SplitTrajsPESsHERE



########################################################################################################################################################################## <= PostTrajectoriesHERE
## Computing Trajectory Statistics -> Cross Sections -> Rate Coefficients
##
function PostTrajectoriesHERE
{


  ###########################################################################################################
  ## Defining Executables for Fortran Code
  ##
  TrajectoriesStatsCommand="coarseair-trajectoriesstats.x"
  PostTrajectoriesCommand="coarseair-posttrajectories.x"


  echo "      [PostTrajectoriesHERE]: "
  echo "      [PostTrajectoriesHERE]: --- Molecule 1, Level/Bin " ${iLevel1} " ----------------------------- "
  echo "      [PostTrajectoriesHERE]: ----- Molecule 2, Level/Bin = " ${iLevel2} " --------------------------- "
  echo "      [PostTrajectoriesHERE]:         Postprocessing Trajectories for iLevel1 = "${iLevel1}" and for iLevel2 = "${iLevel2}", @ iProc = "${iProc}". Command: "${PostTrajectoriesCommand}
  

  ###########################################################################################################
  ## Defining the Run Folder Corresponding to the Current Process
  ##
  if [ ${TranFlg} -eq 0 ]; then 
    COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
  elif [ ${TranFlg} -eq 1 ]; then 
    COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
  elif [ ${TranFlg} -eq 2 ]; then 
    echo "    [PostTrajectoriesHERE]: ERROR! Postprocessing for Gaussian Translational Energy NOT IMPLEMENTED YET! "
    stop
  else 
    echo "    [PostTrajectoriesHERE]: ERROR! Wrong Model for the Tranlational Energy! "
    stop
  fi

   

  ###########################################################################################################
  ## Defining Files for:
  ##      - Trajectories Output
  ##      - Nb of Trajectories in the Trajectory File
  ##      - Error File 
  ## 
  #TrajFile=${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.out.${iPES} # FOR COMPATIBILITY WITH CG-QCT CODE
  TrajFile=${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv
  NTrajFile=${COARSEAIR_BIN_OUTPUT_DIR}/'NConvTraj.dat'
  TrajErrorFile=${COARSEAIR_OUTPUT_DIR}/"RatesErrors_Node"${iNode}"_Proc"${iProc}".dat"
  if [ -f ${TrajFile} ]; then                                                              
    NTraj=$(wc -l < ${TrajFile})                                                        
    NTraj=$((NTraj-1))
  else 
    NTraj=0
  fi
  echo ${NTraj} > ${NTrajFile}
  echo "      [PostTrajectoriesHERE]:         Tot Nb of Trajectories: "${NTraj}


  ###########################################################################################################
  ## Deciding Whether to Split the Trajectories Based on Their PESs
  ## 
  if [ ${SplitPESsFlg} -eq 0 ]; then 
    NPESs=0
    iPESEnd=0
    iPES=0
  else
    iPESEnd=$((iPESStart + NPESs - 1))
    # if [[ ${NTraj} -gt 0 ]]; then     
    #   echo "      [PostTrajectoriesHERE]:         Calling SplitTrajsPESs"
      
    #   ########################################################################################################################################################################## <= SplitTrajsPESsHERE
    #   SplitTrajsPESsHERE
    #   ########################################################################################################################################################################## <= SplitTrajsPESsHERE

    # fi
    iPES=${iPESStart}
  fi


  ###########################################################################################################
  ## Looping on PESs
  ## 
  while [ ${iPES} -le ${iPESEnd} ]; do
    start=`date +%s`


    ###########################################################################################################
    ## Checking Whether the Run Folder Exists
    ## 
    if [ -e "$COARSEAIR_BIN_OUTPUT_DIR" ]; then
      cd ${COARSEAIR_BIN_OUTPUT_DIR}
      
      ###########################################################################################################
      ## Checking Whether the Trajectory File Existed and Getting Nb of Trajectories
      ## 
      if [ -e "${NTrajFile}" ]; then                                                                                        
        NTraj=`cat ${NTrajFile}`
        #typeset -i NTraj=$(cat ${NTrajFile})


        if [ ${NTraj} -gt 0 ] || [ ${BinaryTrajFlg} -gt 0 ]; then
          ###########################################################################################################
          ## If the Trajectory File actually contains Trajectories, then Compute Cross Sections and Rate Coefficients
          ## 

          echo "      [PostTrajectoriesHERE]:         Computing Cross Sections from Trajectories. Command:      "${TrajectoriesStatsCommand}
          eval ${TrajectoriesStatsCommand} ${Tran} ${Tint} ${BinaryTrajFlg} ${iPES}

          echo "      [PostTrajectoriesHERE]:         Computing Rate Coefficients from Cross Sections. Command: "${TrajectoriesStatsCommand}
          rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/"Post.log"
          eval ${PostTrajectoriesCommand} ${Tran} ${Tint} ${NTraj} ${Velocity} ${iPES} ${iLevel1} ${iLevel2} 
        
        else                                                                                                             
          ###########################################################################################################
          ## If the Trajectory File is Empty, then Report the Current Process in the Error File
          ## 
          echo "      [PostTrajectoriesHERE]:         No Trajectories for Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}

          if [ -e "${TrajErrorFile}" ]; then
            echo ${iLevel1}","${iLevel2}",1" >> ${TrajErrorFile}
          else
          echo "# List of Levels / Bins that Generated Errors during Rates Computation (iLevel, jLevel, Flag). Flag == 1: No Trajectories; Flag == 2: No Trajectory File; Flag == 3: No Process Folder" > ${TrajErrorFile}
            echo ${iLevel1}","${iLevel2}",1" >> ${TrajErrorFile}
          fi

        fi
        
      else
        ###########################################################################################################
        ## If the Trajectory File Does Not Exist, then Report the Current Process in the Error File
        ## 
        echo "      [PostTrajectoriesHERE]:         No Trajectory File for Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}
      
        if [ -e "${TrajErrorFile}" ]; then                                                                                   
          echo ${iLevel1}","${iLevel2}",2" >> ${TrajErrorFile}
        else
          echo "# List of Levels / Bins that Generated Errors during Rates Computation (iLevel, jLevel, Flag). Flag == 1: No Trajectories; Flag == 2: No Trajectory File; Flag == 3: No Process Folder" > ${TrajErrorFile}
          echo ${iLevel1}","${iLevel2}",2" >> ${TrajErrorFile}
        fi
        
      fi                                                                                                                  
      
    else
      ###########################################################################################################
      ## If the Process Folder Does Not Exist, then Report the Current Process in the Error File
      ##   
      echo "      [PostTrajectoriesHERE]:          No Process Folder for Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}

      if [ -e "$TrajErrorFile" ]; then
        echo ${iLevel1}","${iLevel2}",3" >> ${TrajErrorFile}
      else
          echo "# List of Levels / Bins that Generated Errors during Rates Computation (iLevel, jLevel, Flag). Flag == 1: No Trajectories; Flag == 2: No Trajectory File; Flag == 3: No Process Folder" > ${TrajErrorFile}
        echo ${iLevel1}","${iLevel2}",3" >> ${TrajErrorFile}
      fi
      
    fi
    

    ###########################################################################################################
    ## Cleaning Files for the Current PES
    ##
    if [ ${SplitPESsFlg} -eq 1 ]; then 
      # rm -rf ${TrajFile}
      # rm -rf ${NTrajFile}
      rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/statistics*.{iPES}
      echo "      [PostTrajectoriesHERE]: ------- PES "${iPES} " --------------- DONE - "
    fi
  

    end=`date +%s`
    runtime=$((end-start))

    iPES=$((iPES+1))
  done   


  ###########################################################################################################
  ## Cleaning Files and Folders 
  ## 
  if [ ${RmTrajFlg} -eq 1 ] && [ -f ${COARSEAIR_BIN_OUTPUT_DIR}/NConvTraj.dat ]; then
    #rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/NConvTraj.dat
    rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/Node*
    rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/statistics*
    rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/*.log
    if [ ${BinaryTrajFlg} -eq 1 ]; then
      rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv*
    fi
  fi



  echo "      [PostTrajectoriesHERE]: ----- Molecule 2, Level/Bin = " ${iLevel2} " ------------------- DONE -- "
  echo "      [PostTrajectoriesHERE]: --- Molecule 1, Level/Bin " ${iLevel1} " --------------------- DONE -- "
  echo " "


}
########################################################################################################################################################################## <= PostTrajectoriesHERE




##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################


###########################################################################################################
## From ComputeRates.sh, PostTrajectoriesAtNode
##    bash  ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} ${iProc} ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc}


COARSEAIR_OUTPUT_DIR=${1}
SplitPESsFlg=${2}
NPESs=${3}
iPESStart=${4}
TranFlg=${5}
Tran=${6}
Tint=${7}
Velocity=${8}
iNode=${9}
iProc=${10}
NMolecules=${11}
SymmFlg=${12}
NLevels1=${13}
MinLevel1=${14}
NLevels2=${15}
MinLevel2=${16}
RmTrajFlg=${17}
BinaryTrajFlg=${18}   
MinProcessInNode=${19}
MaxProcessInNode=${20}
NProcessesPerProc=${21}                                   


echo "    [PostTrajectoriesAtProc.sh]: COARSEAIR_OUTPUT_DIR  = "${COARSEAIR_OUTPUT_DIR}
echo "    [PostTrajectoriesAtProc.sh]: SplitTrajsPESs        = "${SplitTrajsPESs}
echo "    [PostTrajectoriesAtProc.sh]: NPESs                 = "${NPESs}
echo "    [PostTrajectoriesAtProc.sh]: iPESStart             = "${iPESStart}
echo "    [PostTrajectoriesAtProc.sh]: TranFlg               = "${TranFlg}
echo "    [PostTrajectoriesAtProc.sh]: Tran                  = "${Tran}
echo "    [PostTrajectoriesAtProc.sh]: Tint                  = "${Tint}
echo "    [PostTrajectoriesAtProc.sh]: Velocity              = "${Velocity}
echo "    [PostTrajectoriesAtProc.sh]: iNode                 = "${iNode}
echo "    [PostTrajectoriesAtProc.sh]: iProc                 = "${iProc}
echo "    [PostTrajectoriesAtProc.sh]: NMolecules            = "${NMolecules}
echo "    [PostTrajectoriesAtProc.sh]: SymmFlg               = "${SymmFlg}
echo "    [PostTrajectoriesAtProc.sh]: NLevels1              = "${NLevels1}
echo "    [PostTrajectoriesAtProc.sh]: MinLevel1             = "${MinLevel1}
echo "    [PostTrajectoriesAtProc.sh]: NLevels2              = "${NLevels2}
echo "    [PostTrajectoriesAtProc.sh]: MinLevel2             = "${MinLevel2}
echo "    [PostTrajectoriesAtProc.sh]: RmTrajFlg             = "${RmTrajFlg}
echo "    [PostTrajectoriesAtProc.sh]: BinaryTrajFlg         = "${BinaryTrajFlg}
echo "    [PostTrajectoriesAtProc.sh]: MinProcessInNode      = "${MinProcessInNode}
echo "    [PostTrajectoriesAtProc.sh]: MaxProcessInNode      = "${MaxProcessInNode}
echo "    [PostTrajectoriesAtProc.sh]: NProcessesPerProc     = "${NProcessesPerProc}


# echo "    [PostTrajectoriesAtProc.sh]: COARSEAIR_WORKING_DIR = "${COARSEAIR_WORKING_DIR}
# echo "    [PostTrajectoriesAtProc.sh]: COARSEAIR_SH_DIR      = "${COARSEAIR_SH_DIR}
# echo "    [PostTrajectoriesAtProc.sh]: System                = "${System}
# echo "    [PostTrajectoriesAtProc.sh]: NNode                 = "${NNode}
# echo "    [PostTrajectoriesAtProc.sh]: NProc                 = "${NProc}
# echo "    [PostTrajectoriesAtProc.sh]: Molecule1             = "${Molecule1}
# echo "    [PostTrajectoriesAtProc.sh]: MaxLevel1             = "${MaxLevel1}
# echo "    [PostTrajectoriesAtProc.sh]: Molecule2             = "${Molecule2}
# echo "    [PostTrajectoriesAtProc.sh]: MaxLevel2             = "${MaxLevel2}
# echo "    [PostTrajectoriesAtProc.sh]: NProcessesPerNode     = "${NProcessesPerNode}

#source ${COARSEAIR_SH_DIR}/ComputeRates.sh



###########################################################################################################
## Deciding Minimum and Maximum Indexes for the Process to compute in this PROCESSOR
##
MinProcessInProc=$(( ${MinProcessInNode} + $((iProc-1))*${NProcessesPerProc}    ))
MaxProcessInProc=$(( ${MinProcessInNode} +     ${iProc}*${NProcessesPerProc} -1 ))
if [ ${MaxProcessInProc} -gt  ${MaxProcessInNode} ]; then
  MaxProcessInProc=${MaxProcessInNode}
fi

echo "    [PostTrajectoriesAtProc.sh] For Node "${iNode}", Processor "${iProc}", the first Process to be computed is the "${MinProcessInProc}"-th in the List"
echo "    [PostTrajectoriesAtProc.sh] For Node "${iNode}", Processor "${iProc}", the last  Process to be computed is the "${MaxProcessInProc}"-th in the List"



###########################################################################################################
## Deciding Whether to Read the List of Processes to run from the ProcessesToRunList.inp File ...
##     ... or to Go in Increasing Order
##
if [ ${MinLevel1} -eq 0 -a ${MinLevel2} -eq 0 ]; then 
  echo "    [PostTrajectoriesAtProc.sh] Reading Levels/Bins from File "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp
  
  
  ###########################################################################################################
  ## Reading the List of Processes and Checking if the Current is OK for this Processor
  ##
  iProcess=0
  while IFS= read -r ProcessLine
    do
    iProcess=$((iProcess+1))
    if [ ${iProcess} -ge ${MinProcessInProc} ] && [ ${iProcess} -le ${MaxProcessInProc} ]; then
      iLevel1=$(echo $ProcessLine | awk -F'-' '{print $1}')
      if [ -z "$iLevel1" ]; then 
        echo "    [PostTrajectoriesAtProc.sh] -> ERROR! Unable to read the "${iProcess}"-th Line of the Processes-To-Run List ("${COARSEAIR_INPUT_DIR}"/ProcessesToRunList.inp)!"
      fi
      iLevel2=$(echo $ProcessLine | awk -F'-' '{print $2}')
      if [ -z "$iLevel2" ]; then 
        iLevel2=0
      fi

      ############################################################################################################################################################################## <= PostTrajectoriesHERE
      PostTrajectoriesHERE
      ############################################################################################################################################################################## <= PostTrajectoriesHERE

    fi
  done < "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp"


else


  ###########################################################################################################
  ## Going in Increasing Order (Paying Attention to System Symmetries e.g., AB_i + AB_j is the same as AB_j + AB_i) and Deciding if the Current Process is OK for this Processor
  ##
  for (( iProcessesTot=${MinProcessInProc}; iProcessesTot<=${MaxProcessInProc}; iProcessesTot++ )); do
    if [ ${NMolecules} -eq 1 ]; then 
      iLevel1=${iProcessesTot}
      iLevel2=0
    else
      if [ ${SymmFlg} -eq 1 ]; then
          NBetw=0
          iLevel1=1
          while [ ${NBetw} -lt ${iProcessesTot} ] && [ $iLevel1 -le ${NLevels1} ]; do  
            NMin=$(( ${NLevels1} - (${iLevel1}+1) + 1 ))
            NMax=${NLevels1}
            NBetw=$(( (${NMax}+1)*(${NMax})/2 - (${NMin}+1)*(${NMin})/2 ))
            iLevel1=$((${iLevel1} + 1))
          done
          iLevel1=$(( ${iLevel1} - 1))
          iLevel2=$(( ${NLevels1} - (${NBetw} - ${iProcessesTot}) ))
      else
        iLevel1=$(( (${iProcessesTot}-1) / ${NLevels2} + 1 ))
        iLevel1=${iLevel1%.*}
        Temp=$(( (${iLevel1} - 1) * ${NLevels2} ))
        iLevel2=$(( ${iProcessesTot} - ${Temp} ))
      fi
    fi

    ############################################################################################################################################################################## <= PostTrajectoriesHERE
    PostTrajectoriesHERE
    ############################################################################################################################################################################## <= PostTrajectoriesHERE

  done


fi


exit 0

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################