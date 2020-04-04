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
NMolecules=${16}
SymmFlg=${17}
Molecule1=${18}
NLevels1=${19}
MinLevel1=${20}
MaxLevel1=${21}
Molecule2=${22}
NLevels2=${23}
MinLevel2=${24}
MaxLevel2=${25}
RmTrajFlg=${26}
BinaryTrajFlg=${27}   
MinProcessInNode=${28}
MaxProcessInNode=${29}
NProcessesPerNode=${30}
NProcessesPerProc=${31}
                                                         

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
echo "    [PostTrajectoriesAtProc.sh]: NMolecules            = "${NMolecules}
echo "    [PostTrajectoriesAtProc.sh]: SymmFlg               = "${SymmFlg}
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
echo "    [PostTrajectoriesAtProc.sh]: MinProcessInNode      = "${MinProcessInNode}
echo "    [PostTrajectoriesAtProc.sh]: MaxProcessInNode      = "${MaxProcessInNode}
echo "    [PostTrajectoriesAtProc.sh]: NProcessesPerNode     = "${NProcessesPerNode}
echo "    [PostTrajectoriesAtProc.sh]: NProcessesPerProc     = "${NProcessesPerProc}

#source ${COARSEAIR_SH_DIR}/ComputeRates.sh

TrajectoriesStatsCommand="coarseair-trajectoriesstats.x"
PostTrajectoriesCommand="coarseair-posttrajectories.x"


MinProcessInProc=$(( ${MinProcessInNode} + $((iProc-1))*${NProcessesPerProc}    ))
MaxProcessInProc=$(( ${MinProcessInNode} +     ${iProc}*${NProcessesPerProc} -1 ))
if [ ${MaxProcessInProc} -gt  ${MaxProcessInNode} ]; then
  MaxProcessInProc=${MaxProcessInNode}
fi

echo "    [PostTrajectoriesAtProc.sh] For Node "${iNode}", Processor "${iProc}", the first Process to be computed is the "${MinProcessInProc}"-th in the List"
echo "    [PostTrajectoriesAtProc.sh] For Node "${iNode}", Processor "${iProc}", the last  Process to be computed is the "${MaxProcessInProc}"-th in the List"



if [ ${MinLevel1} -eq 0 -a ${MinLevel2} -eq 0 ]; then 
  echo "    [PostTrajectoriesAtProc.sh] Reading Levels/Bins from File "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp
  
  
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

      echo "    [PostTrajectoriesAtProc.sh]: "
      echo "    [PostTrajectoriesAtProc.sh]: --- Molecule 1, Level/Bin " ${iLevel1} " ----------------------------- "
      echo "    [PostTrajectoriesAtProc.sh]: ----- Molecule 2, Level/Bin = " ${iLevel2} " --------------------- "        
      
      if [ ${TranFlg} -eq 0 ]; then 
        COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
      else
        COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
      fi
      

      VelocityFile=${COARSEAIR_OUTPUT_DIR}/Velocity_${Tran}.dat
      if [ -f $exist ]; then
        Velocity=$(sed '2q;d' ${VelocityFile})
        echo "    [PostTrajectoriesAtProc.sh]:    -> Velocity = "${Velocity}
      else
        echo "    [PostTrajectoriesAtProc.sh]:    -> ERROR! Velocity File does exist! CoarseAIR cannot compute Cross Sections!"
        exit1
      fi
      

      ############################################################################################################################################################################## <= PostTrajectories
      if [ ${StochPESFlg} -eq 0 ]; then 
        NPESs=0
        iPESEnd=0
        iPES=0
      else
        iPES=${iPESStart}
        iPESEnd=$((iPESStart + NPESs - 1))
      fi


      while [ ${iPES} -le ${iPESEnd} ]
      do

        if [ ${StochPESFlg} -eq 1 ]; then 
          echo "    [PostTrajectoriesAtProc.sh]: --- PES "${iPES}" ----------------------- "
        fi

        iTry=1
        NLinesTry=0
        while [ ${NLinesTry} -le 10 ] && [ ${iTry} -le 3 ]; do
          
          if [ ${StochPESFlg} -eq 1 ]; then 
            scp ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv.${iPES} ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv
          fi
          if [ -f ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv ]; then                                                              #
          #if [ -f ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.out ]; then                                                              # FOR COMPATIBILITY WITH CG-QCT CODE
            NTraj=$(wc -l < ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv)                                                           #
            #NTraj=$(wc -l < ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.out)                                                           # FOR COMPATIBILITY WITH CG-QCT CODE
            NTraj=$((NTraj-1))
          else 
            NTraj=0
          fi
          echo ${NTraj} > ${COARSEAIR_BIN_OUTPUT_DIR}/'NConvTraj.dat'
          echo "    [PostTrajectoriesAtProc.sh]:    -> Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}": Tot Nb of Trajectories for PES "${iPES}": "${NTraj}


          if [ ${NTraj} -gt 0 ] || [ ${BinaryTrajFlg} -gt 0 ]; then              
              ################################################################################################################################ <= TrajectoriesStats         
              cd ${COARSEAIR_BIN_OUTPUT_DIR}
              
              echo "    [PostTrajectoriesAtProc.sh]:    -> Computing Statistics for Trajectories. Command: "${TrajectoriesStatsCommand}
              eval ${TrajectoriesStatsCommand} ${Tran} ${Tint} ${BinaryTrajFlg}
              echo "    [PostTrajectoriesAtProc.sh]:    -> Done with TrajectoriesStats"      
              ################################################################################################################################ <= TrajectoriesStats        
          fi
          

          #################################################################################################################################### <= FromCrossToRates 
          TrajErrorFile=${COARSEAIR_OUTPUT_DIR}/"RatesErrors_Node"${iNode}"_Proc"${iProc}".dat"

          echo "    [PostTrajectoriesAtProc.sh]:    -> Postprocessing Trajectories for iLevel1 = "${iLevel1}" and for iLevel2 = "${iLevel2}", @ iProc = "${iProc}". Command: "${PostTrajectoriesCommand} ${Tran} 
          start=`date +%s`
          if [ -e "$COARSEAIR_BIN_OUTPUT_DIR" ]; then
          
            cd ${COARSEAIR_BIN_OUTPUT_DIR}
            TrajFile=${COARSEAIR_BIN_OUTPUT_DIR}/'NConvTraj.dat'
          
            if [ -e "$TrajFile" ]; then                                                                                        
          
              typeset -i NTraj=$(cat ${COARSEAIR_BIN_OUTPUT_DIR}/'NConvTraj.dat')

              if [ ${NTraj} -gt 0 ] || [ ${BinaryTrajFlg} -gt 0 ]; then
                rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/"Post.log"
                eval ${PostTrajectoriesCommand} ${Tran} ${Tint} ${NTraj} ${Velocity} ${iPES} ${iLevel1} ${iLevel2} 
              else                                                                                                             
                if [ -e "$TrajErrorFile" ]; then
                  echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
                else
                  echo "# List of Levels / Bins that Generated Errors during Rates Computation" > $TrajErrorFile
                  echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
                fi
                echo "    [PostTrajectoriesAtProc.sh]:    -> No Trajectories for Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}
                #rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/*
              fi
              
            else
            
             if [ -e "$TrajErrorFile" ]; then                                                                                   
                echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
              else
                echo "# List of Levels / Bins that Generated Errors during Rates Computation" > $TrajErrorFile
                echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
              fi
              echo "    [PostTrajectoriesAtProc.sh]:    -> No Trajectories for Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}
              
            fi                                                                                                                  
            
          else
          
            if [ -e "$TrajErrorFile" ]; then
              echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
            else
              echo "# List of Levels / Bins that Generated Errors during Rates Computation" > $TrajErrorFile
              echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
            fi
            echo "    [PostTrajectoriesAtProc.sh]:    -> No Trajectories for Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}
            
          fi
          
          end=`date +%s`
          runtime=$((end-start))
          echo "    [PostTrajectoriesAtProc.sh]: Done with FromCrossToRates for iLevel1 = "${iLevel1}" and for iLevel2 = "${iLevel2}", @ iProc = "${iProc}". RunTime = "${runtime}"s"
          ##################################################################################################################################### <= FromCrossToRates 



          if [ ${StochPESFlg} -eq 1 ]; then 
            NLinesTry=$(wc -l < ${COARSEAIR_OUTPUT_DIR}"/"${System}"/"${Molecule1}"/Rates/T_"${Tran%.*}"_"${Tint%.*}"/Bin"$iLevel1".csv."$iPES)
            iTry=$((iTry+1))
          else
            iTry=4
          fi
        done
          
        if [ ${StochPESFlg} -eq 1 ]; then 
          rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv
          echo "    [PostTrajectoriesAtProc.sh]: --- iPES "${iPES} " --------------- DONE -- "
        fi
        iPES=$((iPES+1))
      done   


      if [ ${StochPESFlg} -eq 1 ]; then
        #scp ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv.Orig ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv
        rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/NConvTraj.dat
        rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/statistics*
      fi
      ############################################################################################################################################################################## <= PostTrajectories



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
      echo "    [PostTrajectoriesAtProc.sh]: --- Molecule 1, Level/Bin " ${iLevel1} " --------------------- DONE -- "
      echo " "
    fi
    
  done < "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp"

else

  iProcessesTot=0
  for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
    
    iLevel2Start=0
    MinLevel2Temp=0
    if [ ${NMolecules} -eq 2 ]; then 
      iLevel2Start=1
      MinLevel2Temp=1
    fi
    if [ ${SymmFlg} -eq 1 ]; then
      iLevel2Start=${iLevel1}
      MinLevel2Temp=${MinLevel1}
    fi

    for (( iLevel2=${iLevel2Start}; iLevel2<=${NLevels2}; iLevel2++ )); do
      iProcessesTot=$(( ${iProcessesTot} + 1 ))
      
      if [ ${iProcessesTot} -ge ${MinProcessInProc} ] && [ ${iProcessesTot} -le ${MaxProcessInProc} ]; then

        echo "    [PostTrajectoriesAtProc.sh]: "
        echo "    [PostTrajectoriesAtProc.sh]: --- Molecule 1, Level/Bin  = " ${iLevel1} " --------------------------- "
        echo "    [PostTrajectoriesAtProc.sh]: ----- Molecule 2, Level/Bin = " ${iLevel2} " --------------------- "

            
        if [ ${TranFlg} -eq 0 ]; then 
          COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
        else
          COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
        fi
        

        VelocityFile=${COARSEAIR_OUTPUT_DIR}/Velocity_${Tran}.dat
        if [ -f $exist ]; then
          Velocity=$(sed '2q;d' ${VelocityFile})
          echo "    [PostTrajectoriesAtProc.sh]:    -> Velocity = "${Velocity}
        else
          echo "    [PostTrajectoriesAtProc.sh]:    -> ERROR! Velocity File does exist! CoarseAIR cannot compute Cross Sections!"
          exit1
        fi
        


        ############################################################################################################################################################################## <= PostTrajectories
        if [ ${StochPESFlg} -eq 0 ]; then 
          NPESs=0
          iPESEnd=0
          iPES=0
        else
          iPES=${iPESStart}
          iPESEnd=$((iPESStart + NPESs - 1))
        fi


        while [ ${iPES} -le ${iPESEnd} ]
        do

          if [ ${StochPESFlg} -eq 1 ]; then 
            echo "    [PostTrajectoriesAtProc.sh]: --- PES "${iPES}" ----------------------- "
          fi

          iTry=1
          NLinesTry=0
          while [ ${NLinesTry} -le 10 ] && [ ${iTry} -le 3 ]; do
            
            if [ ${StochPESFlg} -eq 1 ]; then 
              scp ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv.${iPES} ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv
            fi
            if [ -f ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv ]; then                                                              #
            #if [ -f ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.out ]; then                                                              # FOR COMPATIBILITY WITH CG-QCT CODE
              NTraj=$(wc -l < ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv)                                                           #
              #NTraj=$(wc -l < ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.out)                                                           # FOR COMPATIBILITY WITH CG-QCT CODE
              NTraj=$((NTraj-1))
            else 
              NTraj=0
            fi
            echo ${NTraj} > ${COARSEAIR_BIN_OUTPUT_DIR}/'NConvTraj.dat'
            echo "    [PostTrajectoriesAtProc.sh]:    -> Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}": Tot Nb of Trajectories for PES "${iPES}": "${NTraj}
            
            if [ ${NTraj} -gt 0 ] || [ ${BinaryTrajFlg} -gt 0 ]; then
              ################################################################################################################################ <= TrajectoriesStats         
              cd ${COARSEAIR_BIN_OUTPUT_DIR}
              
              echo "    [PostTrajectoriesAtProc.sh]:    -> Computing Statistics for Trajectories. Command: "${TrajectoriesStatsCommand}
              eval ${TrajectoriesStatsCommand} ${Tran} ${Tint} ${BinaryTrajFlg}
              #echo "    [PostTrajectoriesAtProc.sh]:    -> Done with TrajectoriesStats"              
              ################################################################################################################################ <= TrajectoriesStats         
            fi
            

            #################################################################################################################################### <= FromCrossToRates 
            TrajErrorFile=${COARSEAIR_OUTPUT_DIR}/"RatesErrors_Node"${iNode}"_Proc"${iProc}".dat"

            echo "    [PostTrajectoriesAtProc.sh]:    -> Postprocessing Trajectories for iLevel1 = "${iLevel1}" and for iLevel2 = "${iLevel2}", @ iProc = "${iProc}". Command: "${PostTrajectoriesCommand} ${Tran} 
            start=`date +%s`
            if [ -e "$COARSEAIR_BIN_OUTPUT_DIR" ]; then
            
              cd ${COARSEAIR_BIN_OUTPUT_DIR}
              TrajFile=${COARSEAIR_BIN_OUTPUT_DIR}/'NConvTraj.dat'
            
              if [ -e "$TrajFile" ]; then                                                                                        
            
                typeset -i NTraj=$(cat ${COARSEAIR_BIN_OUTPUT_DIR}/'NConvTraj.dat')

                if [ ${NTraj} -gt 0 ] || [ ${BinaryTrajFlg} -gt 0 ]; then
                  rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/"Post.log"
                  eval ${PostTrajectoriesCommand} ${Tran} ${Tint} ${NTraj} ${Velocity} ${iPES} ${iLevel1} ${iLevel2} 
                else                                                                                                             
                  if [ -e "$TrajErrorFile" ]; then
                    echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
                  else
                    echo "# List of Levels / Bins that Generated Errors during Rates Computation" > $TrajErrorFile
                    echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
                  fi
                  echo "    [PostTrajectoriesAtProc.sh]:    -> No Trajectories for Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}
                  #rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/*
                fi
                
              else
              
               if [ -e "$TrajErrorFile" ]; then                                                                                   
                  echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
                else
                  echo "# List of Levels / Bins that Generated Errors during Rates Computation" > $TrajErrorFile
                  echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
                fi
                echo "    [PostTrajectoriesAtProc.sh]:     -> No Trajectories for Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}
                
              fi                                                                                                                  
              
            else
            
              if [ -e "$TrajErrorFile" ]; then
                echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
              else
                echo "# List of Levels / Bins that Generated Errors during Rates Computation" > $TrajErrorFile
                echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
              fi
              echo "    [PostTrajectoriesAtProc.sh]:     -> No Trajectories for Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}
              
            fi
            
            end=`date +%s`
            runtime=$((end-start))
            #echo "    [PostTrajectoriesAtProc.sh]:    -> Done with FromCrossToRates for iLevel1 = "${iLevel1}" and for iLevel2 = "${iLevel2}", @ iProc = "${iProc}". RunTime = "${runtime}"s"
            ##################################################################################################################################### <= FromCrossToRates 

            
            if [ ${StochPESFlg} -eq 1 ]; then 
              NLinesTry=$(wc -l < ${COARSEAIR_OUTPUT_DIR}"/"${System}"/"${Molecule1}"/Rates/T_"${Tran%.*}"_"${Tint%.*}"/Bin"$iLevel1".csv."$iPES)
              iTry=$((iTry+1))
            else
              iTry=4
            fi
          done
            
          if [ ${StochPESFlg} -eq 1 ]; then 
            rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv
            echo "    [PostTrajectoriesAtProc.sh]: --- iPES "${iPES} " --------------- DONE -- "
          fi
          iPES=$((iPES+1))
        done   


        if [ ${StochPESFlg} -eq 1 ]; then
          #scp ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv.Orig ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv
          rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/NConvTraj.dat
          rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/statistics*
        fi
        ############################################################################################################################################################################## <= PostTrajectories


        
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
        echo "    [PostTrajectoriesAtProc.sh]: --- Molecule 1, Level/Bin = " ${iLevel1} " ------------------- DONE -- "
        echo " "

      fi

    done

  done
  

fi


exit 0