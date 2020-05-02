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


# ------------------------------------------------------------------------------------------------------------- ComputeTrajsPBS #
function ComputeTrajsPBS {

  echo "  [ComputeTrajsPBS]: COARSEAIR_OUTPUT_DIR = "${COARSEAIR_OUTPUT_DIR}
  echo "  [ComputeTrajsPBS]: COARSEAIR_SH_DIR     = "${COARSEAIR_SH_DIR}
  echo "  [ComputeTrajsPBS]: NNode                = "${NNode}
  echo "  [ComputeTrajsPBS]: ProcType             = "${ProcType}
  echo "  [ComputeTrajsPBS]: NProc                = "${NProc}
  echo "  [ComputeTrajsPBS]: NMolecules           = "${NMolecules}
  echo "  [ComputeTrajsPBS]: SymmFlg              = "${SymmFlg}
  echo "  [ComputeTrajsPBS]: NLevels1             = "${NLevels1}
  echo "  [ComputeTrajsPBS]: MinLevel1            = "${MinLevel1}
  echo "  [ComputeTrajsPBS]: MaxLevel1            = "${MaxLevel1}
  echo "  [ComputeTrajsPBS]: NLevels2             = "${NLevels2}
  echo "  [ComputeTrajsPBS]: MinLevel2            = "${MinLevel2}
  echo "  [ComputeTrajsPBS]: MaxLevel2            = "${MaxLevel2}
  

  ###########################################################################################################
  ## Copying the .pbs Template locally
  ## 
  cd ${COARSEAIR_OUTPUT_DIR}/..
  if [ "${ProcType}" = "none" ]; then
    scp ${COARSEAIR_SH_DIR}/../launching/RunTrajectories-Format-UIUC.pbs     ./
  else
    scp ${COARSEAIR_SH_DIR}/../launching/RunTrajectories-Format-Pleiades.pbs ./
  fi


  ###########################################################################################################
  ## Finding out how many PROCESSES we are required to compute and how many per NODE 
  ## 
  if [ ${MinLevel1} -eq 0 -a ${MinLevel2} -eq 0 ]; then 

    ###########################################################################################################
    ## Reading the Process List from ProcessesToRunList.inp
    ## 
    echo "  [ComputeTrajsPBS]: Reading Levels/Bins from File "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp
    NProcessesFromFile=$(($(wc -l < "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp")))
    echo "  [ComputeTrajsPBS]: -> Total Nb of Processes to Run = "${NProcessesFromFile}

    NProcessesPerNode="$(bc <<< "scale = 10; ${NProcessesFromFile} / ${NNode}")"
    NProcessesPerNode="$(echo ${NProcessesPerNode} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
    echo "  [ComputeTrajsPBS]: -> Nb of Processes Per Node = "${NProcessesPerNode}

    MaxProcessAll=${NProcessesFromFile}
    MinProcessAll=1

  else

    ###########################################################################################################
    ## Selecting the Processes to Run in Increasing Order
    ## 
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
        if [ ${iLevel1} -eq ${MinLevel1} ] && [ ${iLevel2} -eq ${MinLevel2Temp} ]; then
          MinProcessAll=${iProcessesTot}
        fi
        if [ ${iLevel1} -eq ${MaxLevel1} ] && [ ${iLevel2} -eq ${MaxLevel2} ]; then
          MaxProcessAll=${iProcessesTot}
        fi
      done
    done
    NProcessesAll=$(( ${MaxProcessAll} - ${MinProcessAll} + 1 ))
    echo "  [ComputeTrajsPBS]: -> Total Nb of Processes to Run = "${NProcessesAll}

    NProcessesPerNode="$(bc <<< "scale = 10; ${NProcessesAll} / ${NNode}")"
    NProcessesPerNode="$(echo ${NProcessesPerNode} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
    echo "  [ComputeTrajsPBS]: -> Nb of Processes Per Node = "${NProcessesPerNode}

  fi
  

  ###########################################################################################################
  ## Customizing ONE .pbs file PER NODE and Submitting EACH of them
  ## 
  for (( iNode=1; iNode<=${NNode}; iNode++ )); do

    MinProcessInNode=$(( ${MinProcessAll} + $((iNode-1))*NProcessesPerNode ))
    if [ ${MinProcessInNode} -le ${MaxProcessAll} ]; then 
      MaxProcessInNode=$(( ${MinProcessAll} + ${iNode}*NProcessesPerNode - 1))
      if [ ${MaxProcessInNode} -gt ${MaxProcessAll} ]; then 
        MaxProcessInNode=${MaxProcessAll}
      fi
      echo "  [ComputeTrajsPBS]: For Node "${iNode}", the first Process to be Computed is the "${MinProcessInNode}"-th"
      echo "  [ComputeTrajsPBS]: For Node "${iNode}", the last  Process to be Computed is the "${MaxProcessInNode}"-th"

      echo "  [ComputeTrajsPBS]: -> iNode = "${iNode}"; MinProcessInNode = "${MinProcessInNode}"; MaxProcessInNode = "${MaxProcessInNode}
              
      if [ "${ProcType}" = "none" ]; then
        sed -e '3s/$/1:ppn='${NProc}'/'                           'RunTrajectories-Format-UIUC.pbs'     > 'RunTrajectoriesTEMP-1.pbs'
      else
        sed -e '2s/$/'${Queue}'/'                                 'RunTrajectories-Format-Pleiades.pbs' > 'RunTrajectoriesTEMP-A.pbs'
        sed -e '3s/$/1:ncpus='${NProc}':model='${ProcType}'/'     'RunTrajectoriesTEMP-A.pbs'  > 'RunTrajectoriesTEMP-B.pbs'
        sed -e '5s/$/'${WallTime}':00:00/'                        'RunTrajectoriesTEMP-B.pbs'  > 'RunTrajectoriesTEMP-1.pbs'
      fi
      sed -e '13s/$/'${NProc}'/'                                  'RunTrajectoriesTEMP-1.pbs'  > 'RunTrajectoriesTEMP-2.pbs'

      sed -e '14s/$/'${SlncFlg}'/'                                'RunTrajectoriesTEMP-2.pbs'  > 'RunTrajectoriesTEMP-3.pbs'
      sed -e '15s/$/'${MergeAllFlg}'/'                            'RunTrajectoriesTEMP-3.pbs'  > 'RunTrajectoriesTEMP-4.pbs'
      sed -e '16s/$/'${RmTrajFlg}'/'                              'RunTrajectoriesTEMP-4.pbs'  > 'RunTrajectoriesTEMP-5.pbs'
      sed -e '17s/$/'${BinaryTrajFlg}'/'                          'RunTrajectoriesTEMP-5.pbs'  > 'RunTrajectoriesTEMP-6.pbs'

      sed -e '4s/$/'${MinProcessInNode}'_'${MaxProcessInNode}'/'  'RunTrajectoriesTEMP-6.pbs'  > 'RunTrajectoriesTEMP-7.pbs'
      sed -e '6s/$/'${MinProcessInNode}'_'${MaxProcessInNode}'/'  'RunTrajectoriesTEMP-7.pbs'  > 'RunTrajectoriesTEMP-8.pbs'

      sed -e '148s/$/'${MinProcessInNode}'/'                      'RunTrajectoriesTEMP-8.pbs'  > 'RunTrajectoriesTEMP-9.pbs' 
      sed -e '149s/$/'${MaxProcessInNode}'/'                      'RunTrajectoriesTEMP-9.pbs'  > 'RunTrajectoriesTEMP-10.pbs' 

      sed -e '150s/$/'${Tran}'/'                                  'RunTrajectoriesTEMP-10.pbs' > 'RunTrajectoriesTEMP-11.pbs' 
      sed -e '151s/$/'${Tint}'/'                                  'RunTrajectoriesTEMP-11.pbs' > 'RunTrajectoriesTEMP-12.pbs'
      sed -e '152s/$/'${iNode}'/'                                 'RunTrajectoriesTEMP-12.pbs' > 'RunTrajectories-'${MinProcessInNode}'-'${MaxProcessInNode}'.pbs'

      qsub ./'RunTrajectories-'${MinProcessInNode}'-'${MaxProcessInNode}'.pbs'
      
      rm -rf ./RunTrajectoriesTEMP*    
    fi

  done
  echo "  [ComputeTrajsPBS]: Done with Submitting PBS Files for computing Trajectories. Now we have to wait for ALL NODES being done."

      
}
#================================================================================================================================#


# ----------------------------------------------------------------------------------------------------------------- ComputeTrajs #
function ComputeTrajs {

  echo " "
  echo "  ======================== COMPUTING TRAJECTORIES ========================== "
  echo " "
  echo "  [ComputeTrajs]: COARSEAIR_WORKING_DIR = "${COARSEAIR_WORKING_DIR}
  echo "  [ComputeTrajs]: COARSEAIR_OUTPUT_DIR  = "${COARSEAIR_OUTPUT_DIR}
  echo "  [ComputeTrajs]: COARSEAIR_SH_DIR      = "${COARSEAIR_SH_DIR}
  echo "  [ComputeTrajs]: NNode                 = "${NNode}
  echo "  [ComputeTrajs]: ParNodes              = "${ParNodes}
  echo "  [ComputeTrajs]: iNode                 = "${iNode}
  echo "  [ComputeTrajs]: NProc                 = "${NProc}
  echo "  [ComputeTrajs]: System                = "${System}
  echo "  [ComputeTrajs]: Tran                  = "${Tran}
  echo "  [ComputeTrajs]: Tint                  = "${Tint}
  echo "  [ComputeTrajs]: NMolecules            = "${NMolecules}
  echo "  [ComputeTrajs]: SymmFlg               = "${SymmFlg}
  echo "  [ComputeTrajs]: NLevels1              = "${NLevels1}
  echo "  [ComputeTrajs]: MinLevel1             = "${MinLevel1}
  echo "  [ComputeTrajs]: MaxLevel1             = "${MaxLevel1}
  echo "  [ComputeTrajs]: NLevels2              = "${NLevels2}
  echo "  [ComputeTrajs]: MinLevel2             = "${MinLevel2}
  echo "  [ComputeTrajs]: MaxLevel2             = "${MaxLevel2}
  echo "  [ComputeTrajs]: MinProcessInNode      = "${MinProcessInNode}
  echo "  [ComputeTrajs]: MaxProcessInNode      = "${MaxProcessInNode}
  echo "  [ComputeTrajs]: SplitPESsFlg          = "${SplitPESsFlg}
  echo "  [ComputeTrajs]: NPESs                 = "${NPESs}
  echo "  [ComputeTrajs]: iPESStart             = "${iPESStart}
  echo "  [ComputeTrajs]: RmTrajFlg             = "${RmTrajFlg}

  if [ ${MinLevel1} -eq 0 -a ${MinLevel2} -eq 0 ]; then 

    echo "  [ComputeTrajs]: Reading Levels/Bins from File "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp
    if [ ${ParNodes} -eq 0 ]; then
      iNode=1
      MinProcessInNode=1
      NProcessesFromFile=$(($(wc -l < "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp")))
      echo "  [ComputeTrajs]: -> Total Nb of Processes to Run = "${NProcessesFromFile}
      MaxProcessInNode=${NProcessesFromFile}
    fi
    echo "  [ComputeTrajs]: For Node "${iNode}", the first Process to be read from file is the "${MinProcessInNode}"-th in the List"
    echo "  [ComputeTrajs]: For Node "${iNode}", the last  Process to be read from file is the "${MaxProcessInNode}"-th in the List"
    

    iProcess=0
    while IFS= read -r ProcessLine
      do
      iProcess=$((iProcess+1))
      
      if [ ${iProcess} -ge ${MinProcessInNode} ] && [ ${iProcess} -le ${MaxProcessInNode} ]; then
        iLevel1=$(echo $ProcessLine | awk -F'-' '{print $1}')
        if [ -z "$iLevel1" ]; then 
          echo "  [ComputeTrajs]: -> ERROR! Unable to read the "${iProcess}"-th Line of the Processes-To-Run List ("${COARSEAIR_INPUT_DIR}"/ProcessesToRunList.inp)!"
        fi
        iLevel2=$(echo $ProcessLine | awk -F'-' '{print $2}')
        if [ -z "$iLevel2" ]; then 
          iLevel2=0
        fi

        echo "  [ComputeTrajs]: -> Found Process Nb "${iProcess}": iLevel1 = "${iLevel1}"; iLevel2 = "${iLevel2}
        ComputeTrajsHERE

      fi
      
    done < "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp"
  

  else

    if [ ${ParNodes} -eq 0 ]; then
      iNode=1
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
          if [ ${iLevel1} -eq ${MinLevel1} ] && [ ${iLevel2} -eq ${MinLevel2Temp} ]; then
            MinProcessInNode=${iProcessesTot}
          fi
          if [ ${iLevel1} -eq ${MaxLevel1} ] && [ ${iLevel2} -eq ${MaxLevel2} ]; then
            MaxProcessInNode=${iProcessesTot}
          fi
        done
      done
      NProcessesAll=$(( ${MaxProcessInNode} - ${MinProcessInNode} + 1 ))
      echo "  [ComputeTrajs]: -> Total Nb of Processes to Run = "${NProcessesAll}
    fi
    echo "  [ComputeTrajs]: For Node "${iNode}", the first Process to be computed is the "${MinProcessInNode}"-th"
    echo "  [ComputeTrajs]: For Node "${iNode}", the last  Process to be computed is the "${MaxProcessInNode}"-th"


    iProcessesTot=0
    for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
      iLevel2Start=0
      if [ ${NMolecules} -eq 2 ]; then 
        iLevel2Start=1
      fi
      if [ ${SymmFlg} -eq 1 ]; then
        iLevel2Start=${iLevel1}
      fi
      for (( iLevel2=${iLevel2Start}; iLevel2<=${NLevels2}; iLevel2++ )); do
        iProcessesTot=$((iProcessesTot+1))

        if [ ${iProcessesTot} -ge ${MinProcessInNode} ] && [ ${iProcessesTot} -le ${MaxProcessInNode} ]; then
          ComputeTrajsHERE
        fi

      done
    done

  
  fi
      
}
#================================================================================================================================#


# -------------------------------------------------------------------------------------------------------- RunTrajectoriesAtNode #
function ComputeTrajsHERE {

  echo "    [ComputeTrajsHERE]: --- Molecule 1, Level/Bin " ${iLevel1} " -------------------------------- "
  echo "    [ComputeTrajsHERE]: ----- Molecule 2, Level/Bin = " ${iLevel2} " ------------------------ "
  
  if [ ${TranFlg} -eq 0 ]; then 
    COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
  else
    COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
  fi
      
  echo "    [ComputeTrajsHERE]: Calling RunTrajectoriesAtNode"
  RunTrajectoriesAtNode
  wait

  echo "    [ComputeTrajsHERE]: Calling MergeTrajectories"
  MergeTrajectories
  wait

  echo "    [ComputeTrajsHERE]: ----- Molecule 2, Level/Bin = " ${iLevel2} " ---------------- DONE -- "
  echo "    [ComputeTrajsHERE]: --- Molecule 1, Level/Bin " ${iLevel1} " ------------------------ DONE -- "
  echo " "
}
#================================================================================================================================#


# -------------------------------------------------------------------------------------------------------- RunTrajectoriesAtNode #
function RunTrajectoriesAtNode {
  
  iNode=1
  
  start=`date +%s`
  echo "    [RunTrajectoriesAtNode]: Parallelizing the Computation of Trajectories"
  echo "    [RunTrajectoriesAtNode]: Node "${iNode}
  
  mkdir -p ${COARSEAIR_BIN_OUTPUT_DIR}/Node_${iNode}
  python3 ${COARSEAIR_SH_DIR}'/SplitLevelsList.py' ${COARSEAIR_BIN_OUTPUT_DIR}/Node_${iNode} ${COARSEAIR_OUTPUT_DIR} ${System} ${NMolecules} ${Molecule1} ${NLevels1} ${iLevel1} ${Molecule2} ${NLevels2} ${iLevel2}

  if [[ ${NProc} -eq 1 ]]; then
    RunTrajectoriesAtProc 1
  elif [[ ${NProc} -eq 2 ]]; then
    parallel --xapply -j 2 "sh ${COARSEAIR_SH_DIR}/RunTrajectoriesAtProc.sh '${System}' '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} ${TranFlg} ${Tran} ${Tint} ${iLevel1} ${iLevel2} " ::: {1..2} 
  elif [[ ${NProc} -eq 4 ]]; then
    parallel --xapply -j 4 "sh ${COARSEAIR_SH_DIR}/RunTrajectoriesAtProc.sh '${System}' '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} ${TranFlg} ${Tran} ${Tint} ${iLevel1} ${iLevel2} " ::: {1..4} 
  elif [[ ${NProc} -eq 8 ]]; then
    parallel --xapply -j 8 "sh ${COARSEAIR_SH_DIR}/RunTrajectoriesAtProc.sh '${System}' '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} ${TranFlg} ${Tran} ${Tint} ${iLevel1} ${iLevel2} " ::: {1..8} 
  elif [[ ${NProc} -eq 16 ]]; then
    parallel --xapply -j 16 "sh ${COARSEAIR_SH_DIR}/RunTrajectoriesAtProc.sh '${System}' '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} ${TranFlg} ${Tran} ${Tint} ${iLevel1} ${iLevel2} " ::: {1..16} 
  elif [[ ${NProc} -eq 20 ]]; then
    parallel --xapply -j 20 "sh ${COARSEAIR_SH_DIR}/RunTrajectoriesAtProc.sh '${System}' '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} ${TranFlg} ${Tran} ${Tint} ${iLevel1} ${iLevel2} " ::: {1..20} 
  elif [[ ${NProc} -eq 24 ]]; then
    parallel --xapply -j 24 "sh ${COARSEAIR_SH_DIR}/RunTrajectoriesAtProc.sh '${System}' '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} ${TranFlg} ${Tran} ${Tint} ${iLevel1} ${iLevel2} " ::: {1..24} 
  elif [[ ${NProc} -eq 32 ]]; then
    parallel --xapply -j 32 "sh ${COARSEAIR_SH_DIR}/RunTrajectoriesAtProc.sh '${System}' '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} ${TranFlg} ${Tran} ${Tint} ${iLevel1} ${iLevel2} " ::: {1..32} 
  elif [[ ${NProc} -eq 48 ]]; then
    parallel --xapply -j 48 "sh ${COARSEAIR_SH_DIR}/RunTrajectoriesAtProc.sh '${System}' '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} ${TranFlg} ${Tran} ${Tint} ${iLevel1} ${iLevel2} " ::: {1..48} 
  elif [[ ${NProc} -eq 64 ]]; then
    parallel --xapply -j 64 "sh ${COARSEAIR_SH_DIR}/RunTrajectoriesAtProc.sh '${System}' '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} ${TranFlg} ${Tran} ${Tint} ${iLevel1} ${iLevel2} " ::: {1..64} 
  elif [[ ${NProc} -eq 128 ]]; then
    parallel --xapply -j 128 "sh ${COARSEAIR_SH_DIR}/RunTrajectoriesAtProc.sh '${System}' '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} ${TranFlg} ${Tran} ${Tint} ${iLevel1} ${iLevel2} " ::: {1..128} 
  else
    echo "ERROR: Number of Precessors not Coherent with RunTrajectoriesAtNode! (Check ComputeRates.sh)"
    exit 1
  fi
  
  wait
  end=`date +%s`
  runtime=$((end-start))
  echo "    [RunTrajectoriesAtNode]: Done Running Trajectories. RunTime = "${runtime}"s"
  
}
#================================================================================================================================#


# -------------------------------------------------------------------------------------------------------- RunTrajectoriesAtProc #
function RunTrajectoriesAtProc() {
  
  iProc=${1}
  NProcTot=${NProc}
  
  echo "      [RunTrajectoriesAtProc]: COARSEAIR_OUTPUT_DIR = "${COARSEAIR_OUTPUT_DIR}
  echo "      [RunTrajectoriesAtProc]: NNode                = "${NNode}
  echo "      [RunTrajectoriesAtProc]: iNode                = "${iNode}
  echo "      [RunTrajectoriesAtProc]: iProc                = "${iProc}
  echo "      [RunTrajectoriesAtProc]: NProcTot             = "${NProcTot}
  echo "      [RunTrajectoriesAtProc]: Tran                 = "${Tran}
  echo "      [RunTrajectoriesAtProc]: Tint                 = "${Tint}
  echo "      [RunTrajectoriesAtProc]: iLevel1              = "${iLevel1}
  echo "      [RunTrajectoriesAtProc]: iLevel2              = "${iLevel2}
  
  RunTrajectoriesCommand="coarseair-runtrajectories.x"


  mkdir -p ${COARSEAIR_BIN_OUTPUT_DIR}/Node_${iNode}
  mkdir -p ${COARSEAIR_BIN_OUTPUT_DIR}/Node_${iNode}/Proc_${iProc}
  cd ${COARSEAIR_BIN_OUTPUT_DIR}/Node_${iNode}/Proc_${iProc}
  scp ../levels* ./

  echo "      [RunTrajectoriesAtProc]: Running Trajectories for Node "${iNode}" of "${NNode}", Processor "${iProc}
  eval ${RunTrajectoriesCommand} '${COARSEAIR_BIN_OUTPUT_DIR}' ${Tran} ${Tint} ${NProcTot} ${iProc} ${NNode} ${iNode} ${iLevel1} ${iLevel2}
  echo "      [RunTrajectoriesAtProc]: Done Running Trajectories for Node "${iNode}", Processor "${iProc}
  
}
#================================================================================================================================#


# ------------------------------------------------------------------------------------------------------------ MergeTrajectories #
function MergeTrajectories {

  cd ${COARSEAIR_BIN_OUTPUT_DIR}

  echo "    [MergeTrajectories]: Merging Trajectories Output."
  #iNode=1
  #while [ ${iNode} -le ${NNode} ]; do
  #  echo "    [MergeTrajectories]: Merging for iNode = "${iNode}

    iProc=1
    while [ $iProc -le ${NProc} ]; do  
      echo "    [MergeTrajectories]: Merging for iProc = "${iProc}
      
      if [ -f ./Node_$iNode/Proc_$iProc/trajectories.csv ]; then
      
        if [ -f ./trajectories.csv ]; then
      
          tail -n+2 ./Node_$iNode/Proc_$iProc/trajectories.csv >> ./trajectories.csv
          if [ -f ./Node_$iNode/Proc_$iProc/PaQSOl.out ]; then
            tail -n+2 ./Node_$iNode/Proc_$iProc/PaQSOl.out >> ./PaQSOl.out
            rm -rf ./Node_$iNode/Proc_$iProc/PaQSOl.out
          fi
          
        else
        
          cat ./Node_${iNode}/Proc_${iProc}/trajectories.csv > ./trajectories.csv
          if [ -f ./Node_$iNode/Proc_$iProc/PaQSOl.out ]; then
            cat ./Node_$iNode/Proc_$iProc/PaQSOl.out > ./PaQSOl.out
            rm -rf ./Node_$iNode/Proc_$iProc/PaQSOl.out
          fi
          
        fi
        
        #rm -rf ./Node_$iNode/Proc_$iProc/trajectories.csv

      fi
      
    iProc=$((iProc+1))
    done

    if [ ${RmTrajFlg} -eq 1 ]; then
      rm -rf ./Node_$iNode
      echo "    [MergeTrajectories]: -> Done with merging Trajectories for Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel1}". Now I will remove the Node Folder"
    fi
    
  #iNode=$((iNode+1))
  #done
      
  if [ -f ./trajectories.csv ]; then
    NTraj=$(wc -l < ./trajectories.csv)
    NTraj=$((NTraj-1))
  else 
    NTraj=0
  fi
  echo ${NTraj} > ${COARSEAIR_BIN_OUTPUT_DIR}/'NConvTraj.dat'
  echo "    [MergeTrajectories]: -> Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel1}": Tot Nb of Converged Trajectories: "${NTraj}
  
}
#================================================================================================================================#
#================================================================================================================================#


# # ------------------------------------------------------------------------------------------------------------------------------ #
# # --------------------------------------------------------------------------------------------------------------- SplitTrajsPESs #
# function SplitTrajsPESs {

#   echo " "
#   echo "        ======================== SPLITTING TRAJECTORIES ========================== "
#   echo "        [SplitTrajsPESs]: COARSEAIR_BIN_OUTPUT_DIR = "${COARSEAIR_BIN_OUTPUT_DIR}
#   echo "        [SplitTrajsPESs]: iPESStart                = "${iPESStart}
#   echo "        [SplitTrajsPESs]: iPESEnd                  = "${iPESEnd}
#   echo "        [SplitTrajsPESs]: -> Splitting the Trajectory File Based on PESs"

#   iPES=${iPESStart}
#   while [[ ${iPES} -le ${iPESEnd} ]]; do
#     PESFile=${COARSEAIR_BIN_OUTPUT_DIR}'/trajectories.csv.'${iPES}
#     echo '#    iTraj, iPES,       bmax,        b_i,           j_i,           v_i,         arr_i,           j_f,           v_f,         arr_f' > ${PESFile} 
#     iPES=$(($iPES+1))
#   done        

#   OrigFile=$COARSEAIR_BIN_OUTPUT_DIR/trajectories.csv
#   OLDIFS=$IFS
#   IFS=','
#   [ ! -f ${OrigFile} ] && { echo "$OrigFile file not found"; exit 99; }
#   sed 1d ${OrigFile} | while read Idx iPES bmax b_i j_i v_i arr_i j_f v_f arr_f
#   do
#     jPES=$(($iPES+0))
#     if [[ ${jPES} -ge ${iPESStart} ]] && [[ ${jPES} -le ${iPESEnd} ]]; then
#       PESFile=${COARSEAIR_BIN_OUTPUT_DIR}'/trajectories.csv.'${jPES}
#       echo $Idx','$iPES','$bmax','$b_i','$j_i','$v_i','$arr_i','$j_f','$v_f','$arr_f >> ${PESFile} 
#     fi
#   done 
#   IFS=$OLDIFS

#   echo "        [SplitTrajsPESs]: -> Done Splitting PESs in the Trajectory Files"
#   echo "        ========================================================================== "
#   echo " "

# }
# #================================================================================================================================#
# #================================================================================================================================#


# ------------------------------------------------------------------------------------------------------------------------------ #
#----------------------------------------------------------------------------------------------------------- PostTrajectoriesPBS #
function PostTrajectoriesPBS {

  echo "  [PostTrajectoriesPBS]: COARSEAIR_OUTPUT_DIR = "${COARSEAIR_OUTPUT_DIR}
  echo "  [PostTrajectoriesPBS]: COARSEAIR_SH_DIR     = "${COARSEAIR_SH_DIR}
  echo "  [PostTrajectoriesPBS]: NNode                = "${NNode}
  echo "  [PostTrajectoriesPBS]: ProcType             = "${ProcType}
  echo "  [PostTrajectoriesPBS]: NProc                = "${NProc}
  echo "  [PostTrajectoriesPBS]: NMolecules           = "${NMolecules}
  echo "  [PostTrajectoriesPBS]: MinLevel1            = "${MinLevel1}
  echo "  [PostTrajectoriesPBS]: MaxLevel1            = "${MaxLevel1}
  echo "  [PostTrajectoriesPBS]: MinLevel2            = "${MinLevel2}
  echo "  [PostTrajectoriesPBS]: MaxLevel2            = "${MaxLevel2}
  echo "  [PostTrajectoriesPBS]: BinaryTrajFlg        = "${BinaryTrajFlg}


  ###########################################################################################################
  ## Copying the .pbs Template locally
  ## 
  cd ${COARSEAIR_OUTPUT_DIR}/..
  if [ "${ProcType}" = "none" ]; then
    scp ${COARSEAIR_SH_DIR}/../launching/PostTrajectories-Format-UIUC.pbs     ./
  else
    scp ${COARSEAIR_SH_DIR}/../launching/PostTrajectories-Format-Pleiades.pbs ./
  fi


  ###########################################################################################################
  ## Finding out how many PROCESSES we are required to compute and how many per NODE 
  ## 
  if [ ${MinLevel1} -eq 0 -a ${MinLevel2} -eq 0 ]; then 
    
    ###########################################################################################################
    ## Reading the Process List from ProcessesToRunList.inp
    ## 

    echo "  [PostTrajectoriesPBS]: Reading Levels/Bins from File "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp
    NProcessesFromFile=$(($(wc -l < "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp")))
    echo "  [PostTrajectoriesPBS]: -> Total Nb of Processes to Run = "${NProcessesFromFile}

    NProcessesPerNode="$(bc <<< "scale = 10; ${NProcessesFromFile} / ${NNode}")"
    NProcessesPerNode="$(echo ${NProcessesPerNode} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
    echo "  [PostTrajectoriesPBS]: -> Nb of Processes Per Node = "${NProcessesPerNode}

    MinProcessAll=1
    MaxProcessAll=${NProcessesFromFile}

  else
    
    ###########################################################################################################
    ## Selecting the Processes to Run in Increasing Order
    ## 

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
        if [ ${iLevel1} -eq ${MinLevel1} ] && [ ${iLevel2} -eq ${MinLevel2Temp} ]; then
          MinProcessAll=${iProcessesTot}
        fi
        if [ ${iLevel1} -eq ${MaxLevel1} ] && [ ${iLevel2} -eq ${MaxLevel2} ]; then
          MaxProcessAll=${iProcessesTot}
        fi
      done
    done
    NProcessesAll=$(( ${MaxProcessAll} - ${MinProcessAll} + 1 ))
    echo "  [PostTrajectoriesPBS]: -> Total Nb of Processes to Run = "${NProcessesAll}

    NProcessesPerNode="$(bc <<< "scale = 10; ${NProcessesAll} / ${NNode}")"
    NProcessesPerNode="$(echo ${NProcessesPerNode} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
    echo "  [PostTrajectoriesPBS]: -> Nb of Processes Per Node = "${NProcessesPerNode}

  fi
  

  ###########################################################################################################
  ## Customizing ONE .pbs file PER NODE and Submitting EACH of them
  ## 
  for (( iNode=1; iNode<=${NNode}; iNode++ )); do

    MinProcessInNode=$(( ${MinProcessAll} + $((${iNode}-1))*NProcessesPerNode     ))
    if [ ${MinProcessInNode} -le ${MaxProcessAll} ]; then 
      MaxProcessInNode=$(( ${MinProcessAll} +        ${iNode}*NProcessesPerNode - 1 ))
      if [ ${MaxProcessInNode} -gt ${MaxProcessAll} ]; then 
        MaxProcessInNode=${MaxProcessAll}
      fi
      
      echo "  [PostTrajectoriesPBS]: -> iNode = "${iNode}"; MinProcessInNode = "${MinProcessInNode}"; MaxProcessInNode = "${MaxProcessInNode}
          
      if [ "${ProcType}" = "none" ]; then
        sed -e '3s/$/1:ppn='${NProc}'/'                              'PostTrajectories-Format-UIUC.pbs'     > 'PostTrajectoriesTEMP-1.pbs'
      else
        sed -e '2s/$/'${Queue}'/'                                    'PostTrajectories-Format-Pleiades.pbs' > 'PostTrajectoriesTEMP-A.pbs'
        sed -e '3s/$/1:ncpus='${NProc}':model='${ProcType}'/'        'PostTrajectoriesTEMP-A.pbs'  > 'PostTrajectoriesTEMP-B.pbs'
        sed -e '5s/$/'${WallTime}':00:00/'                           'PostTrajectoriesTEMP-B.pbs'  > 'PostTrajectoriesTEMP-1.pbs'
      fi
      sed -e '13s/$/'${NProc}'/'                                     'PostTrajectoriesTEMP-1.pbs'  > 'PostTrajectoriesTEMP-2.pbs'

      sed -e '14s/$/'${SlncFlg}'/'                                   'PostTrajectoriesTEMP-2.pbs'  > 'PostTrajectoriesTEMP-3.pbs'
      sed -e '15s/$/'${MergeAllFlg}'/'                               'PostTrajectoriesTEMP-3.pbs'  > 'PostTrajectoriesTEMP-4.pbs'
      sed -e '16s/$/'${RmTrajFlg}'/'                                 'PostTrajectoriesTEMP-4.pbs'  > 'PostTrajectoriesTEMP-5.pbs'
      sed -e '17s/$/'${BinaryTrajFlg}'/'                             'PostTrajectoriesTEMP-5.pbs'  > 'PostTrajectoriesTEMP-6.pbs'

      sed -e '4s/$/'${MinProcessInNode}'-'${MaxProcessInNode}'/'     'PostTrajectoriesTEMP-6.pbs'  > 'PostTrajectoriesTEMP-7.pbs'
      sed -e '6s/$/'${MinProcessInNode}'_'${MaxProcessInNode}'/'     'PostTrajectoriesTEMP-7.pbs'  > 'PostTrajectoriesTEMP-8.pbs'

      sed -e '148s/$/'${MinProcessInNode}'/'                         'PostTrajectoriesTEMP-8.pbs'  > 'PostTrajectoriesTEMP-9.pbs' 
      sed -e '149s/$/'${MaxProcessInNode}'/'                         'PostTrajectoriesTEMP-9.pbs'  > 'PostTrajectoriesTEMP-10.pbs' 
      
      sed -e '150s/$/'${Tran}'/'                                     'PostTrajectoriesTEMP-10.pbs' > 'PostTrajectoriesTEMP-11.pbs' 
      sed -e '151s/$/'${Tint}'/'                                     'PostTrajectoriesTEMP-11.pbs' > 'PostTrajectoriesTEMP-12.pbs'

      sed -e '152s/$/'${iNode}'/'                                    'PostTrajectoriesTEMP-12.pbs' > 'PostTrajectories-'${MinProcessInNode}'-'${MaxProcessInNode}'.pbs'

      qsub ./'PostTrajectories-'${MinProcessInNode}'-'${MaxProcessInNode}'.pbs'

      rm -rf ./PostTrajectoriesTEMP*
    fi

  done
  echo "  [PostTrajectoriesPBS]: Done with Submitting PBS Files for Postprocessing Trajectories. Now we have to wait for ALL NODES being done."
      

}
#================================================================================================================================#
#================================================================================================================================#


# ------------------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------- PostTrajectoriesAtNode #
function PostTrajectoriesAtNode {
  
  echo " "
  echo "  ==================== POSTPROCESSING TRAJECTORIES ========================= "
  echo " "
  echo "  [PostTrajectoriesAtNode]: COARSEAIR_WORKING_DIR = "${COARSEAIR_WORKING_DIR}
  echo "  [PostTrajectoriesAtNode]: COARSEAIR_OUTPUT_DIR  = "${COARSEAIR_OUTPUT_DIR}
  echo "  [PostTrajectoriesAtNode]: COARSEAIR_SH_DIR      = "${COARSEAIR_SH_DIR}
  echo "  [PostTrajectoriesAtNode]: NNode                 = "${NNode}
  echo "  [PostTrajectoriesAtNode]: ParNodes              = "${ParNodes}
  echo "  [PostTrajectoriesAtNode]: iNode                 = "${iNode}
  echo "  [PostTrajectoriesAtNode]: NProc                 = "${NProc}
  echo "  [PostTrajectoriesAtNode]: System                = "${System}
  echo "  [PostTrajectoriesAtNode]: Tran                  = "${Tran}
  echo "  [PostTrajectoriesAtNode]: Tint                  = "${Tint}
  echo "  [PostTrajectoriesAtNode]: NMolecules            = "${NMolecules}
  echo "  [PostTrajectoriesAtNode]: SymmFlg               = "${SymmFlg}
  echo "  [PostTrajectoriesAtNode]: Molecule1             = "${Molecule1}
  echo "  [PostTrajectoriesAtNode]: NLevels1              = "${NLevels1}
  echo "  [PostTrajectoriesAtNode]: MinLevel1             = "${MinLevel1}
  echo "  [PostTrajectoriesAtNode]: MaxLevel1             = "${MaxLevel1}
  echo "  [PostTrajectoriesAtNode]: Molecule2             = "${Molecule2}
  echo "  [PostTrajectoriesAtNode]: NLevels2              = "${NLevels2}
  echo "  [PostTrajectoriesAtNode]: MinLevel2             = "${MinLevel2}
  echo "  [PostTrajectoriesAtNode]: MaxLevel2             = "${MaxLevel2}
  echo "  [PostTrajectoriesAtNode]: MinProcessInNode      = "${MinProcessInNode}
  echo "  [PostTrajectoriesAtNode]: MaxProcessInNode      = "${MaxProcessInNode}
  echo "  [PostTrajectoriesAtNode]: SplitPESsFlg          = "${SplitPESsFlg}
  echo "  [PostTrajectoriesAtNode]: NPESs                 = "${NPESs}
  echo "  [PostTrajectoriesAtNode]: iPESStart             = "${iPESStart}
  echo "  [PostTrajectoriesAtNode]: RmTrajFlg             = "${RmTrajFlg}
  echo "  [PostTrajectoriesAtNode]: BinaryTrajFlg         = "${BinaryTrajFlg}

  start=`date +%s`

  VelocityFile=${COARSEAIR_OUTPUT_DIR}/Velocity_${Tran}.dat
  if [ -f $exist ]; then
    Velocity=$(sed '2q;d' ${VelocityFile})
    echo "  [PostTrajectoriesAtNode]: Velocity = "${Velocity}
  else
    echo "  [PostTrajectoriesAtNode]: ERROR! Velocity File does exist! CoarseAIR cannot compute Cross Sections!"
    exit1
  fi


  if [ ${MinLevel1} -eq 0 -a ${MinLevel2} -eq 0 ]; then 

    echo "  [PostTrajectoriesAtNode]: Reading Levels/Bins from File "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp
    if [ ${ParNodes} -eq 0 ]; then
      iNode=1
      MinProcessInNode=1
      NProcessesFromFile=$(($(wc -l < "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp")))
      MaxProcessInNode=${NProcessesFromFile}
    fi
    NProcessesPerNode=$(( ${MaxProcessInNode} - ${MinProcessInNode} + 1 ))
    echo "  [PostTrajectoriesAtNode]: NProcessesPerNode = "${NProcessesPerNode}

    NProcessesPerProc="$(bc <<< "scale = 10; ${NProcessesPerNode} / ${NProc}")"
    NProcessesPerProc="$(echo ${NProcessesPerProc} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
    echo "  [PostTrajectoriesAtNode]: NProcessesPerProc = "${NProcessesPerProc}


  else

    if [ ${ParNodes} -eq 0 ]; then
      iNode=1
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
          if [ ${iLevel1} -eq ${MinLevel1} ] && [ ${iLevel2} -eq ${MinLevel2Temp} ]; then
            MinProcessInNode=${iProcessesTot}
          fi
          if [ ${iLevel1} -eq ${MaxLevel1} ] && [ ${iLevel2} -eq ${MaxLevel2} ]; then
            MaxProcessInNode=${iProcessesTot}
          fi
        done
      done
    fi

    NProcessesPerNode=$(( ${MaxProcessInNode} - ${MinProcessInNode} + 1 ))
    echo "  [PostTrajectoriesAtNode]: NProcessesPerNode = "${NProcessesPerNode}
    echo "  [PostTrajectoriesAtNode]: MinProcessInNode  = "${MinProcessInNode}
    echo "  [PostTrajectoriesAtNode]: MaxProcessInNode  = "${MaxProcessInNode}

 
    NProcessesPerProc="$(bc <<< "scale = 10; ${NProcessesPerNode} / ${NProc}")"
    NProcessesPerProc="$(echo ${NProcessesPerProc} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
    echo "  [PostTrajectoriesAtNode]: -> Parallelizing the Postprocessing on "${NProc}" Processors. Nb of Processes per Processor = "${NProcessesPerProc}
    echo "  [PostTrajectoriesAtNode]: -> Wating for ALL the Processors being done"


  fi


  # for (( iProc=1; iProc<=${NProc}; iProc++ )); do
  #   bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} ${iProc} ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc} &
  # done
                  

  if [[ ${NProc} -eq 1 ]]; then
   bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh  ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} 1 ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc}
  elif [[ ${NProc} -eq 2 ]]; then
    parallel --xapply -j 2 -u "bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} {1} ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc}" ::: {1..2} 
  elif [[ ${NProc} -eq 4 ]]; then
    parallel --xapply -j 4 -u "bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} {1} ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc}" ::: {1..4} 
  elif [[ ${NProc} -eq 8 ]]; then
    parallel --xapply -j 8 -u "bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} {1} ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc}" ::: {1..8} 
  elif [[ ${NProc} -eq 12 ]]; then
    parallel --xapply -j 12 -u "bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} {1} ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc}" ::: {1..12} 
  elif [[ ${NProc} -eq 16 ]]; then
    parallel --xapply -j 16 -u "bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} {1} ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc}" ::: {1..16} 
  elif [[ ${NProc} -eq 20 ]]; then
    parallel --xapply -j 20 -u "bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} {1} ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc}" ::: {1..20} 
  elif [[ ${NProc} -eq 24 ]]; then
    parallel --xapply -j 24 -u "bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} {1} ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc}" ::: {1..24} 
  elif [[ ${NProc} -eq 32 ]]; then
    parallel --xapply -j 32 -u "bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} {1} ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc}" ::: {1..32} 
  elif [[ ${NProc} -eq 64 ]]; then
    parallel --xapply -j 64 -u "bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} {1} ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc}" ::: {1..64} 
  elif [[ ${NProc} -eq 128 ]]; then
    parallel --xapply -j 128 -u "bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_OUTPUT_DIR} ${SplitPESsFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${iNode} {1} ${NMolecules} ${SymmFlg} ${NLevels1} ${MinLevel1} ${NLevels2} ${MinLevel2} ${RmTrajFlg} ${BinaryTrajFlg} ${MinProcessInNode} ${MaxProcessInNode} ${NProcessesPerProc}" ::: {1..128} 
  else
    echo "ERROR: Number of Precessors not Coherent with PostTrajectoriesAtNode! (Check ComputeRates.sh)"
    exit 1
  fi

  wait

  end=`date +%s`
  runtime=$((end-start))
  echo "  [PostTrajectoriesAtNode]: Done Postprocessing Trajectories. RunTime = "${runtime}"s"      
      
}
#================================================================================================================================#



# # ----------------------------------------------------------------------------------------------------------------- ComputeRates #
# function MergeAllRates {

#   for Tran in "${Tran_vec[@]}"; do :
#     echo '  [MergeAllRates] '
#     echo "  [MergeAllRates]: ===   Translational Temperature      = " ${Tran} " ============================ "
#     #for Tint in "${Tint_Vec[@]}"; do :
#       Tint=${Tran}
#       echo "  [MergeAllRates]: ===== Internal Temperature           = " ${Tint} " ======================== "
#       echo '  [MergeAllRates]'     

#       iProcessesTot=0
#       for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
#         iLevel2Start=0
#         if [ ${NMolecules} -eq 2 ]; then 
#           iLevel2Start=1
#         fi
#         if [ ${SymmFlg} -eq 1 ]; then
#           iLevel2Start=${iLevel1}
#         fi
#         for (( iLevel2=${iLevel2Start}; iLevel2<=${NLevels2}; iLevel2++ )); do
#           iProcessesTot=$((iProcessesTot+1))

#           if [ ${iProcessesTot} -ge ${MinProcessInNode} ] && [ ${iProcessesTot} -le ${MaxProcessInNode} ]; then

#             if [ ${TranFlg} -eq 0 ]; then 
#               COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
#             else
#               COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
#             fi
#             cd ${COARSEAIR_BIN_OUTPUT_DIR}


#             if [ -f ../trajectories-Tot.csv ]; then  
#               tail -n+2 ./trajectories.csv >> ../trajectories-Tot.csv
#               if [ -f ./PaQSOl.out ]; then
#                 echo "### Molecule 1, Level / Bin Nb "$iLevel1"; Molecule 2, Level / Bin Nb "$iLevel2 >> ../PaQSOl-Tot.out
#                 tail -n+2 ./PaQSOl.out >> ../PaQSOl-Tot.out
#                 rm -rf ./PaQSOl.out
#               fi
#             else
#               cat ./trajectories.csv > ../trajectories-Tot.csv
#               if [ -f ./PaQSOl.out ]; then
#                 echo "### Molecule 1, Level / Bin Nb "$iLevel1"; Molecule 2, Level / Bin Nb "$iLevel2 >> ../PaQSOl-Tot.out
#                 tail -n+2 ./PaQSOl.out >> ../PaQSOl-Tot.out
#                 rm -rf ./PaQSOl.out
#               fi
#             fi
            
#             rm -rf ./trajectories.csv

#           fi

#         done
#       done

#     #done
#   echo "  [MergeAllRates]: ------------------------------------------------------------ "
#   done
#   wait

# }
# #================================================================================================================================#
