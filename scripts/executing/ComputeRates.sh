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
  echo "  [ComputeTrajsPBS]: NLevels1             = "${NLevels1}
  echo "  [ComputeTrajsPBS]: MinLevel1            = "${MinLevel1}
  echo "  [ComputeTrajsPBS]: MaxLevel1            = "${MaxLevel1}
  echo "  [ComputeTrajsPBS]: NLevels2             = "${NLevels2}
  echo "  [ComputeTrajsPBS]: MinLevel2            = "${MinLevel2}
  echo "  [ComputeTrajsPBS]: MaxLevel2            = "${MaxLevel2}
  

  cd ${COARSEAIR_OUTPUT_DIR}/..
  if [ "${ProcType}" = "none" ]; then
    scp ${COARSEAIR_SH_DIR}/RunTrajectories-Format-UIUC.pbs     ./
  elif [ "${ProcType}" = "test" ]; then
    scp ${COARSEAIR_SH_DIR}/RunTrajectories-Format-Test.pbs     ./
  else
    scp ${COARSEAIR_SH_DIR}/RunTrajectories-Format-Pleiades.pbs ./
  fi


  if [ ${MinLevel1} -eq 0 -a ${MinLevel2} -eq 0 ]; then 
    echo "  [ComputeTrajsPBS]: Reading Levels/Bins from File "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp

    NProcessesFromFile=$(($(wc -l < "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp")+1))
    NProcessesPerNode="$(bc <<< "scale = 10; ${NProcessesFromFile} / ${NNode}")"
    NProcessesPerNode="$(echo ${NProcessesPerNode} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"

    for (( iNode=1; iNode<=${NNode}; iNode++ )); do

      MinProcessInNode=$(($((iNode-1))*NProcessesPerNode+1))
      if [ ${MinProcessInNode} -le ${NProcessesFromFile} ]; then 
        MaxProcessInNode=$((iNode*NProcessesPerNode))
        if [ ${MaxProcessInNode} -gt ${NProcessesFromFile} ]; then 
          MaxProcessInNode=${NProcessesFromFile}
        fi
        Tempp=0
        echo "  [ComputeTrajsPBS]: For Node "${iNode}", the first Process to be read from file is the "${MinProcessInNode}"-th in the List"
        echo "  [ComputeTrajsPBS]: For Node "${iNode}", the last  Process to be read from file is the "${MaxProcessInNode}"-th in the List"

        echo "  [ComputeTrajsPBS]: -> iNode = "${iNode}"; MinLevel1 = 0 and MaxLevel1 = "${MinProcessInNode}"; MinLevel2 = 0 and MaxLevel2 = "${MaxProcessInNode}
                
        if [ "${ProcType}" = "none" ]; then
          sed -e '3s/$/1:ppn='${NProc}'/'                           'RunTrajectories-Format-UIUC.pbs'     > 'RunTrajectoriesTEMP-1.pbs'
        if [ "${ProcType}" = "test" ]; then
          sed -e '3s/$/1:ncpus='${NProc}':model='${ProcType}'/'     'RunTrajectories-Format-Test.pbs'     > 'RunTrajectoriesTEMP-1.pbs'
        else
          sed -e '3s/$/1:ncpus='${NProc}':model='${ProcType}'/'     'RunTrajectories-Format-Pleiades.pbs' > 'RunTrajectoriesTEMP-1.pbs'
        fi
        sed -e '12s/$/'${NProc}'/'                                  'RunTrajectoriesTEMP-1.pbs'  > 'RunTrajectoriesTEMP-2.pbs'
        
        sed -e '4s/$/'${MinProcessInNode}'_'${MaxProcessInNode}'/'  'RunTrajectoriesTEMP-2.pbs'  > 'RunTrajectoriesTEMP-3.pbs'
        sed -e '6s/$/'${MinProcessInNode}'_'${MaxProcessInNode}'/'  'RunTrajectoriesTEMP-3.pbs'  > 'RunTrajectoriesTEMP-4.pbs'

        sed -e '147s/$/'${Tempp}'/'                                 'RunTrajectoriesTEMP-4.pbs'  > 'RunTrajectoriesTEMP-5.pbs'
        sed -e '148s/$/'${MinProcessInNode}'/'                      'RunTrajectoriesTEMP-5.pbs'  > 'RunTrajectoriesTEMP-6.pbs' 
        sed -e '149s/$/'${Tempp}'/'                                 'RunTrajectoriesTEMP-6.pbs'  > 'RunTrajectoriesTEMP-7.pbs'
        sed -e '150s/$/'${MaxProcessInNode}'/'                      'RunTrajectoriesTEMP-7.pbs'  > 'RunTrajectoriesTEMP-8.pbs' 

        sed -e '151s/$/'${Tran}'/'                                  'RunTrajectoriesTEMP-8.pbs'  > 'RunTrajectoriesTEMP-9.pbs' 
        sed -e '152s/$/'${Tint}'/'                                  'RunTrajectoriesTEMP-9.pbs'  > 'RunTrajectoriesTEMP-10.pbs'

        sed -e '153s/$/'${iNode}'/'                                 'RunTrajectoriesTEMP-10.pbs' > 'RunTrajectories-'${MinProcessInNode}'-'${MaxProcessInNode}'.pbs'

        qsub ./'RunTrajectories-'${MinProcessInNode}'-'${MaxProcessInNode}'.pbs'
        
        rm -rf ./RunTrajectoriesTEMP*    
      fi

    done

  else

    MinLevel1Temp=${MinLevel1}
    MaxLevel1Temp=${MaxLevel1}
    MinLevel2Temp=${MinLevel2}
    MaxLevel2Temp=${MaxLevel2}
    NProcessesTot=0
    for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
      if [ ${iLevel1} -ge ${MinLevel1Temp} ] && [ ${iLevel1} -le ${MaxLevel1Temp} ]; then
        if [ ${SymmFlg} -eq 1 ]; then
          MinLevel2Temp=${iLevel1}
        fi
        for (( iLevel2=1; iLevel2<=${NLevels2}; iLevel2++ )); do
          if [ ${iLevel2} -ge ${MinLevel2Temp} ] && [ ${iLevel2} -le ${MaxLevel2Temp} ]; then
            NProcessesTot=$((NProcessesTot+1))
          fi
        done
      fi
    done  
    echo "  [ComputeTrajsPBS]: -> Total Nb of Processes to Run = "${NProcessesTot}


    NProcessesPerNode="$(bc <<< "scale = 10; ${NProcessesTot} / ${NNode}")"
    NProcessesPerNode="$(echo ${NProcessesPerNode} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
    echo "  [ComputeTrajsPBS]: -> Nb of Processes per Node = "${NProcessesPerNode}
    echo "  [ComputeTrajsPBS] "


    MinLevel1Temp=${MinLevel1}
    MaxLevel1Temp=${MaxLevel1}
    MinLevel2Temp=${MinLevel2}
    if [ ${SymmFlg} -eq 1 ]; then
      MinLevel2Temp=${MinLevel1}
    fi
    MaxLevel2Temp=${MaxLevel2}
    MinLevel1Node=${MinLevel1Temp}
    MinLevel2Node=${MinLevel2Temp}
    iProcessesTot=0
    iProcessesPerNode=0
    iNode=0
    for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
      if [ ${iLevel1} -ge ${MinLevel1Temp} ] && [ ${iLevel1} -le ${MaxLevel1Temp} ]; then
        if [ ${SymmFlg} -eq 1 ]; then
          MinLevel2Temp=${iLevel1}
        fi
        for (( iLevel2=1; iLevel2<=${NLevels2}; iLevel2++ )); do
          if [ ${iLevel2} -ge ${MinLevel2Temp} ] && [ ${iLevel2} -le ${MaxLevel2Temp} ]; then
            iProcessesTot=$((iProcessesTot+1))
            iProcessesPerNode=$((iProcessesPerNode+1))

            if [ ${iProcessesPerNode} -eq ${NProcessesPerNode} ] || [ ${iProcessesTot} -eq ${NProcessesTot} ]; then 
              iNode=$((iNode+1))
              MaxLevel1Node=${iLevel1}
              MaxLevel2Node=${iLevel2}

              echo "  [ComputeTrajsPBS]: -> iNode = "${iNode}"; MinLevel1Node = "${MinLevel1Node}"; MaxLevel1Node = "${MaxLevel1Node}"; MinLevel2Node = "${MinLevel2Node}"; MaxLevel2Node = "${MaxLevel2Node}
              
              if [ "${ProcType}" = "none" ]; then
                sed -e '3s/$/1:ppn='${NProc}'/'                        'RunTrajectories-Format-UIUC.pbs'     > 'RunTrajectoriesTEMP-1.pbs'
              else
                sed -e '3s/$/1:ncpus='${NProc}':model='${ProcType}'/'  'RunTrajectories-Format-Pleiades.pbs' > 'RunTrajectoriesTEMP-1.pbs'
              fi
              sed -e '12s/$/'${NProc}'/'                               'RunTrajectoriesTEMP-1.pbs'  > 'RunTrajectoriesTEMP-2.pbs'
              
              sed -e '4s/$/'${MinLevel1Node}'_'${MinLevel2Node}'/'     'RunTrajectoriesTEMP-2.pbs'  > 'RunTrajectoriesTEMP-3.pbs'
              sed -e '6s/$/'${MinLevel1Node}'_'${MinLevel2Node}'/'     'RunTrajectoriesTEMP-3.pbs'  > 'RunTrajectoriesTEMP-4.pbs'
       
              sed -e '147s/$/'${MinLevel1Node}'/'                      'RunTrajectoriesTEMP-4.pbs'  > 'RunTrajectoriesTEMP-5.pbs'
              sed -e '148s/$/'${MaxLevel1Node}'/'                      'RunTrajectoriesTEMP-5.pbs'  > 'RunTrajectoriesTEMP-6.pbs' 
              sed -e '149s/$/'${MinLevel2Node}'/'                      'RunTrajectoriesTEMP-6.pbs'  > 'RunTrajectoriesTEMP-7.pbs'
              sed -e '150s/$/'${MaxLevel2Node}'/'                      'RunTrajectoriesTEMP-7.pbs'  > 'RunTrajectoriesTEMP-8.pbs' 

              sed -e '151s/$/'${Tran}'/'                               'RunTrajectoriesTEMP-8.pbs'  > 'RunTrajectoriesTEMP-9.pbs' 
              sed -e '152s/$/'${Tint}'/'                               'RunTrajectoriesTEMP-9.pbs'  > 'RunTrajectoriesTEMP-10.pbs'

              sed -e '153s/$/'${iNode}'/'                              'RunTrajectoriesTEMP-10.pbs' > 'RunTrajectories-'${MinLevel1Node}'-'${MinLevel2Node}'.pbs'

              qsub ./'RunTrajectories-'${MinLevel1Node}'-'${MinLevel2Node}'.pbs'
              
              rm -rf ./RunTrajectoriesTEMP*


              iProcessesPerNode=0
              MinLevel1Node=${MaxLevel1Node}
              if [ ${MaxLevel2Node} -ge ${MaxLevel2} ]; then
                MinLevel1Node=$((MaxLevel1Node+1))
                if [ ${SymmFlg} -eq 1 ]; then
                  MinLevel2Node=${MinLevel1Node}
                else
                  MinLevel2Node=${MinLevel2}
                fi
              else
                MinLevel2Node=$((MaxLevel2Node+1))
              fi
            fi

          fi
        done
      fi
    done  

  fi
  
  
  #rm -rf ./RunTrajectories-*
  
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
  echo "  [ComputeTrajs]: iNode                 = "${iNode}
  echo "  [ComputeTrajs]: NProc                 = "${NProc}
  echo "  [ComputeTrajs]: System                = "${System}
  echo "  [ComputeTrajs]: Tran                  = "${Tran}
  echo "  [ComputeTrajs]: Tint                  = "${Tint}
  echo "  [ComputeTrajs]: NMolecules            = "${NMolecules}
  echo "  [ComputeTrajs]: NLevels1              = "${NLevels1}
  echo "  [ComputeTrajs]: MinLevel1             = "${MinLevel1}
  echo "  [ComputeTrajs]: MaxLevel1             = "${MaxLevel1}
  echo "  [ComputeTrajs]: NLevels2              = "${NLevels2}
  echo "  [ComputeTrajs]: MinLevel2             = "${MinLevel2}
  echo "  [ComputeTrajs]: MaxLevel2             = "${MaxLevel2}
  echo "  [ComputeTrajs]: StochPESFlg           = "${StochPESFlg}
  echo "  [ComputeTrajs]: NPESs                 = "${NPESs}
  echo "  [ComputeTrajs]: iPESStart             = "${iPESStart}

  if [ ${MinLevel1} -eq 0 -a ${MinLevel2} -eq 0 ]; then 
    echo "  [ComputeTrajs]: Reading Levels/Bins from File "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp

    if [ ${NNode} -gt 1 ]; then
      MinProcessInNode=${MaxLevel1}
      MaxProcessInNode=${MaxLevel2}
    else
      
      iNode=1
      MinProcessInNode=1
      NProcessesFromFile=$(($(wc -l < "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp")+1))
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
        MinLevel1Temp=${MinLevel1}
        MaxLevel1Temp=${MaxLevel1}
        MinLevel2Temp=${MinLevel2}
        MaxLevel2Temp=${MaxLevel2}

        echo "  [ComputeTrajs]: -> Found Process Nb "${iProcess}": iLevel1 = "${iLevel1}"; iLevel2 = "${iLevel2}
        
        echo "  [ComputeTrajs]: --- Molecule 1, Level/Bin " ${iLevel1} " -------------------------------- "
        echo "  [ComputeTrajs]: ----- Molecule 2, Level/Bin = " ${iLevel2} " ------------------------ "
        
        if [ ${TranFlg} -eq 0 ]; then 
          COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
        else
          COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
        fi
            
        echo "  [ComputeTrajs]: Calling RunTrajectoriesAtNode"
        RunTrajectoriesAtNode
        wait

        echo "  [ComputeTrajs]: Calling MergeTrajectories"
        MergeTrajectories
        wait

        echo "  [ComputeTrajs]: ----- Molecule 2, Level/Bin = " ${iLevel2} " ---------------- DONE -- "
        echo "  [ComputeTrajs]: --- Molecule 1, Level/Bin " ${iLevel1} " ------------------------ DONE -- "
        echo " "

      fi
      
    done < "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp"
  

  else
  
    MinLevel1Temp=${MinLevel1}
    MaxLevel1Temp=${MaxLevel1}
    MinLevel2Temp=${MinLevel2}
    MaxLevel2Temp=${MaxLevel2}
    for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
      if [ ${iLevel1} -ge ${MinLevel1Temp} ] && [ ${iLevel1} -le ${MaxLevel1Temp} ]; then
        if [ ${SymmFlg} -eq 1 ]; then
          MinLevel2Temp=${iLevel1}
        fi
        for (( iLevel2=1; iLevel2<=${NLevels2}; iLevel2++ )); do
          if [ ${iLevel2} -ge ${MinLevel2Temp} ] && [ ${iLevel2} -le ${MaxLevel2Temp} ]; then
        
            echo "  [ComputeTrajs]: --- Molecule 1, Level/Bin " ${iLevel1} " -------------------------------- "
            echo "  [ComputeTrajs]: ----- Molecule 2, Level/Bin = " ${iLevel2} " ------------------------ "
            echo "  [ComputeTrajs]"
          
            if [ ${TranFlg} -eq 0 ]; then 
              COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
            else
              COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
            fi
          
          
            echo "  [ComputeTrajs]: Calling RunTrajectoriesAtNode"
            RunTrajectoriesAtNode
            wait

            echo "  [ComputeTrajs]: Calling MergeTrajectories"
            MergeTrajectories
            wait

            echo "  [ComputeTrajs]: ----- Molecule 2, Level/Bin = " ${iLevel2} " ---------------- DONE -- "
            echo "  [ComputeTrajs]: --- Molecule 1, Level/Bin " ${iLevel1} " ------------------------ DONE -- "
            echo " "
          fi

        done
      fi
    done         
  
  fi
      
}
#================================================================================================================================#


# -------------------------------------------------------------------------------------------------------- RunTrajectoriesAtNode #
function RunTrajectoriesAtNode {
  
  iNode=1
  
  start=`date +%s`
  echo "    [RunTrajectoriesAtNode]: Parallelizing the Computation of Trajectories"
  echo "    [RunTrajectoriesAtNode]: Node "${iNode}
  
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


# ------------------------------------------------------------------------------------------------------------------------------ #
# --------------------------------------------------------------------------------------------------------------- SplitTrajsPESs #
function SplitTrajsPESs {

  echo " "
  echo "  ======================== SPLITTING TRAJECTORIES ========================== "
  echo " "
  echo "  [SplitTrajsPESs]: COARSEAIR_OUTPUT_DIR     = "${COARSEAIR_OUTPUT_DIR}
  echo "  [SplitTrajsPESs]: COARSEAIR_SH_DIR         = "${COARSEAIR_SH_DIR}
  echo "  [SplitTrajsPESs]: TranFlg                  = "${TranFlg}
  echo "  [SplitTrajsPESs]: Tran                     = "${Tran}
  echo "  [SplitTrajsPESs]: Tint                     = "${Tint}
  echo "  [SplitTrajsPESs]: iPESStart                = "${iPESStart}
  echo "  [SplitTrajsPESs]: NMolecules               = "${NMolecules}
  echo "  [SplitTrajsPESs]: MinLevel1                = "${MinLevel1}
  echo "  [SplitTrajsPESs]: MaxLevel1                = "${MaxLevel1}
  echo "  [SplitTrajsPESs]: MinLevel2                = "${MinLevel2}
  echo "  [SplitTrajsPESs]: MaxLevel2                = "${MaxLevel2}

  echo "  [SplitTrajsPESs]: -> Splitting PESs in the Trajectory Files"

  MinLevel1Temp=${MinLevel1}
  MaxLevel1Temp=${MaxLevel1}
  MinLevel2Temp=${MinLevel2}
  MaxLevel2Temp=${MaxLevel2}
  for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
    if [ ${iLevel1} -ge ${MinLevel1Temp} ] && [ ${iLevel1} -le ${MaxLevel1Temp} ]; then
      if [ ${SymmFlg} -eq 1 ]; then
        MinLevel2Temp=${iLevel1}
      fi
      for (( iLevel2=1; iLevel2<=${NLevels2}; iLevel2++ )); do
        if [ ${iLevel2} -ge ${MinLevel2Temp} ] && [ ${iLevel2} -le ${MaxLevel2Temp} ]; then

          echo "  [SplitTrajsPESs]: Splitting Trajectories for iLevel1 = "${iLevel1}" and iLevel2 = "${iLevel2}

          if [ ${TranFlg} -eq 0 ]; then 
            COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
          else
            COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
          fi

          python3 ${COARSEAIR_SH_DIR}/SplitTrajsPESs.py ${COARSEAIR_BIN_OUTPUT_DIR} ${iPESStart}

        fi
      done
    fi
  done

  echo "  [SplitTrajsPESs]: -> Done Splitting PESs in the Trajectory Files"
  echo "  ========================================================================== "

}
#================================================================================================================================#
#================================================================================================================================#


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


  cd ${COARSEAIR_OUTPUT_DIR}/..
  if [ "${ProcType}" = "none" ]; then
    scp ${COARSEAIR_SH_DIR}/PostTrajectories-Format-UIUC.pbs     ./
  else
    scp ${COARSEAIR_SH_DIR}/PostTrajectories-Format-Pleiades.pbs ./
  fi


  if [ ${MinLevel1} -eq 0 -a ${MinLevel2} -eq 0 ]; then 
    echo "  [PostTrajectoriesPBS]: Reading Levels/Bins from File "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp

    NProcessesFromFile=$(($(wc -l < "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp")+1))
    NProcessesPerNode="$(bc <<< "scale = 10; ${NProcessesFromFile} / ${NNode}")"
    NProcessesPerNode="$(echo ${NProcessesPerNode} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"


    for (( iNode=1; iNode<=${NNode}; iNode++ )); do

      MinProcessInNode=$(($((iNode-1))*NProcessesPerNode+1))
      if [ ${MinProcessInNode} -le ${NProcessesFromFile} ]; then 
        MaxProcessInNode=$((iNode*NProcessesPerNode))
        if [ ${MaxProcessInNode} -gt ${NProcessesFromFile} ]; then 
          MaxProcessInNode=${NProcessesFromFile}
        fi
        Tempp=0

        echo "  [PostTrajectoriesPBS]: -> iNode = "${iNode}"; MinLevel1 = 0 and MaxLevel1 = "${MinProcessInNode}"; MinLevel2 = 0 and MaxLevel2 = "${MaxProcessInNode}
            
        if [ "${ProcType}" = "none" ]; then
          sed -e '3s/$/1:ppn='${NProc}'/'                              'PostTrajectories-Format-UIUC.pbs'     > 'PostTrajectoriesTEMP-1.pbs'
        else
          sed -e '3s/$/1:ncpus='${NProc}':model='${ProcType}'/'        'PostTrajectories-Format-Pleiades.pbs' > 'PostTrajectoriesTEMP-1.pbs'
        fi
        sed -e '12s/$/'${NProc}'/'                                     'PostTrajectoriesTEMP-1.pbs'  > 'PostTrajectoriesTEMP-2.pbs'
        
        sed -e '4s/$/'${MinProcessInNode}'_'${MaxProcessInNode}'/'     'PostTrajectoriesTEMP-2.pbs'  > 'PostTrajectoriesTEMP-3.pbs'
        sed -e '6s/$/'${MinProcessInNode}'_'${MaxProcessInNode}'/'     'PostTrajectoriesTEMP-3.pbs'  > 'PostTrajectoriesTEMP-4.pbs'

        sed -e '147s/$/'${Tempp}'/'                                    'PostTrajectoriesTEMP-4.pbs'  > 'PostTrajectoriesTEMP-5.pbs'
        sed -e '148s/$/'${MinProcessInNode}'/'                         'PostTrajectoriesTEMP-5.pbs'  > 'PostTrajectoriesTEMP-6.pbs' 

        sed -e '149s/$/'${Tempp}'/'                                    'PostTrajectoriesTEMP-6.pbs'  > 'PostTrajectoriesTEMP-7.pbs'
        sed -e '150s/$/'${MaxProcessInNode}'/'                         'PostTrajectoriesTEMP-7.pbs'  > 'PostTrajectoriesTEMP-8.pbs' 
        
        sed -e '151s/$/'${Tran}'/'                                     'PostTrajectoriesTEMP-8.pbs'  > 'PostTrajectoriesTEMP-9.pbs' 
        sed -e '152s/$/'${Tint}'/'                                     'PostTrajectoriesTEMP-9.pbs'  > 'PostTrajectoriesTEMP-10.pbs'

        sed -e '153s/$/'${iNode}'/'                                    'PostTrajectoriesTEMP-10.pbs' > 'PostTrajectories-'${MinProcessInNode}'-'${MaxProcessInNode}'.pbs'

        qsub ./'PostTrajectories-'${MinProcessInNode}'-'${MaxProcessInNode}'.pbs'
        
        rm -rf ./PostTrajectoriesTEMP*

      fi

    done

  else

    MinLevel1Temp=${MinLevel1}
    MaxLevel1Temp=${MaxLevel1}
    MinLevel2Temp=${MinLevel2}
    MaxLevel2Temp=${MaxLevel2}
    NProcessesTot=0
    for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
      if [ ${iLevel1} -ge ${MinLevel1Temp} ] && [ ${iLevel1} -le ${MaxLevel1Temp} ]; then
        if [ ${SymmFlg} -eq 1 ]; then
          MinLevel2Temp=${iLevel1}
        fi
        for (( iLevel2=1; iLevel2<=${NLevels2}; iLevel2++ )); do
          if [ ${iLevel2} -ge ${MinLevel2Temp} ] && [ ${iLevel2} -le ${MaxLevel2Temp} ]; then
            NProcessesTot=$((NProcessesTot+1))
          fi
        done
      fi
    done  
    echo "  [PostTrajectoriesPBS]: -> Total Nb of Processes to Run = "${NProcessesTot}


    NProcessesPerNode="$(bc <<< "scale = 10; ${NProcessesTot} / ${NNode}")"
    NProcessesPerNode="$(echo ${NProcessesPerNode} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
    echo "  [PostTrajectoriesPBS]: -> Nb of Processes per Node = "${NProcessesPerNode}
    echo "  [PostTrajectoriesPBS] "


    MinLevel1Temp=${MinLevel1}
    MaxLevel1Temp=${MaxLevel1}
    MinLevel2Temp=${MinLevel2}
    if [ ${SymmFlg} -eq 1 ]; then
      MinLevel2Temp=${MinLevel1}
    fi
    MaxLevel2Temp=${MaxLevel2}
    MinLevel1Node=${MinLevel1Temp}
    MinLevel2Node=${MinLevel2Temp}
    iProcessesTot=0
    iProcessesPerNode=0
    iNode=0
    for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
      if [ ${iLevel1} -ge ${MinLevel1Temp} ] && [ ${iLevel1} -le ${MaxLevel1Temp} ]; then
        if [ ${SymmFlg} -eq 1 ]; then
          MinLevel2Temp=${iLevel1}
        fi
        for (( iLevel2=1; iLevel2<=${NLevels2}; iLevel2++ )); do
          if [ ${iLevel2} -ge ${MinLevel2Temp} ] && [ ${iLevel2} -le ${MaxLevel2Temp} ]; then
            iProcessesTot=$((iProcessesTot+1))
            iProcessesPerNode=$((iProcessesPerNode+1))

            if [ ${iProcessesPerNode} -eq ${NProcessesPerNode} ] || [ ${iProcessesTot} -eq ${NProcessesTot} ]; then 
              iNode=$((iNode+1))
              MaxLevel1Node=${iLevel1}
              MaxLevel2Node=${iLevel2}
              
              echo "  [PostTrajectoriesPBS]: -> iNode = "${iNode}"; MinLevel1Node = "${MinLevel1Node}"; MaxLevel1Node = "${MaxLevel1Node}"; MinLevel2Node = "${MinLevel2Node}"; MaxLevel2Node = "${MaxLevel2Node}
          
              if [ "${ProcType}" = "none" ]; then
                sed -e '3s/$/1:ppn='${NProc}'/'                        'PostTrajectories-Format-UIUC.pbs'     > 'PostTrajectoriesTEMP-1.pbs'
              elif [ "${ProcType}" = "test" ]; then
                sed -e '3s/$/1:ncpus='${NProc}':model='${ProcType}'/'  'PostTrajectories-Format-Test.pbs'     > 'PostTrajectoriesTEMP-1.pbs'
              else
                sed -e '3s/$/1:ncpus='${NProc}':model='${ProcType}'/'  'PostTrajectories-Format-Pleiades.pbs' > 'PostTrajectoriesTEMP-1.pbs'
              fi
              sed -e '12s/$/'${NProc}'/'                               'PostTrajectoriesTEMP-1.pbs'  > 'PostTrajectoriesTEMP-2.pbs'
              
              sed -e '4s/$/'${MinLevel1Node}'_'${MinLevel2Node}'/'     'PostTrajectoriesTEMP-2.pbs'  > 'PostTrajectoriesTEMP-3.pbs'
              sed -e '6s/$/'${MinLevel1Node}'_'${MinLevel2Node}'/'     'PostTrajectoriesTEMP-3.pbs'  > 'PostTrajectoriesTEMP-4.pbs'

              sed -e '147s/$/'${MinLevel1Node}'/'                      'PostTrajectoriesTEMP-4.pbs'  > 'PostTrajectoriesTEMP-5.pbs'
              sed -e '148s/$/'${MaxLevel1Node}'/'                      'PostTrajectoriesTEMP-5.pbs'  > 'PostTrajectoriesTEMP-6.pbs' 

              sed -e '149s/$/'${MinLevel2Node}'/'                      'PostTrajectoriesTEMP-6.pbs'  > 'PostTrajectoriesTEMP-7.pbs'
              sed -e '150s/$/'${MinLevel2Node}'/'                      'PostTrajectoriesTEMP-7.pbs'  > 'PostTrajectoriesTEMP-8.pbs' 
              
              sed -e '151s/$/'${Tran}'/'                               'PostTrajectoriesTEMP-8.pbs'  > 'PostTrajectoriesTEMP-9.pbs' 
              sed -e '152s/$/'${Tint}'/'                               'PostTrajectoriesTEMP-9.pbs'  > 'PostTrajectoriesTEMP-10.pbs'

              sed -e '153s/$/'${iNode}'/'                              'PostTrajectoriesTEMP-10.pbs' > 'PostTrajectories-'${MinLevel1Node}'-'${MinLevel2Node}'.pbs'

              qsub ./'PostTrajectories-'${MinLevel1Node}'-'${MinLevel2Node}'.pbs'
              
              rm -rf ./PostTrajectoriesTEMP*

              iProcessesPerNode=0
              MinLevel1Node=${MaxLevel1Node}
              if [ ${MaxLevel2Node} -ge ${MaxLevel2} ]; then
                MinLevel1Node=$((MaxLevel1Node+1))
                if [ ${SymmFlg} -eq 1 ]; then
                  MinLevel2Node=${MinLevel1Node}
                else
                  MinLevel2Node=${MinLevel2}
                fi
              else
                MinLevel2Node=$((MaxLevel2Node+1))
              fi
            fi

          fi
        done
      fi
    done  


  fi
  
  #rm -rf ./PostTrajectories-*
  
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
  echo "  [PostTrajectoriesAtNode]: iNode                 = "${iNode}
  echo "  [PostTrajectoriesAtNode]: NProc                 = "${NProc}
  echo "  [PostTrajectoriesAtNode]: System                = "${System}
  echo "  [PostTrajectoriesAtNode]: Tran                  = "${Tran}
  echo "  [PostTrajectoriesAtNode]: Tint                  = "${Tint}
  echo "  [PostTrajectoriesAtNode]: Molecule1             = "${Molecule1}
  echo "  [PostTrajectoriesAtNode]: NLevels1              = "${NLevels1}
  echo "  [PostTrajectoriesAtNode]: MinLevel1             = "${MinLevel1}
  echo "  [PostTrajectoriesAtNode]: MaxLevel1             = "${MaxLevel1}
  echo "  [PostTrajectoriesAtNode]: Molecule2             = "${Molecule2}
  echo "  [PostTrajectoriesAtNode]: NLevels2              = "${NLevels2}
  echo "  [PostTrajectoriesAtNode]: MinLevel2             = "${MinLevel2}
  echo "  [PostTrajectoriesAtNode]: MaxLevel2             = "${MaxLevel2}
  echo "  [PostTrajectoriesAtNode]: StochPESFlg           = "${StochPESFlg}
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
    

    if [ ${NNode} -gt 1 ]; then
      MinProcessInNode=${MaxLevel1}
      MaxProcessInNode=${MaxLevel2}
    else
      iNode=1
      MinProcessInNode=1
      NProcessesFromFile=$(($(wc -l < "${COARSEAIR_INPUT_DIR}/ProcessesToRunList.inp")+1))
      MaxProcessInNode=${NProcessesFromFile}
    fi

    NProcessesPerNode=$((${MaxProcessInNode} - ${MinProcessInNode} + 1))

    NProcessesPerProc="$(bc <<< "scale = 10; ${NProcessesPerNode} / ${NProc}")"
    NProcessesPerProc="$(echo ${NProcessesPerProc} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"

    for (( iProc=1; iProc<=${NProc}; iProc++ )); do
      MinProcessInProc=$(($((iProc-1))*NProcessesPerProc+1))
      MaxProcessInProc=$((iProc*NProcessesPerProc))
      echo "  [PostTrajectoriesAtNode]: For Node "${iNode}", Proc "${iProc}", the first Process to be read from file is the "${MinProcessInProc}"-th in the List"
      echo "  [PostTrajectoriesAtNode]: For Node "${iNode}", Proc "${iProc}", the last  Process to be read from file is the "${MaxProcessInProc}"-th in the List"

      echo "  [PostTrajectoriesAtNode]: Postprocessing Trajectories @ iProc = "${iProc}"; MinLevel1 = 0 and MaxLevel1 = "${MinProcessInProc}"; MinLevel2 = 0 and MaxLevel2 = "${MaxProcessInProc}
      bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_WORKING_DIR} ${COARSEAIR_OUTPUT_DIR} ${COARSEAIR_SH_DIR} ${System} ${StochPESFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${NNode} ${iNode} ${NProc} ${iProc} ${NMolecules} ${SymmFlg} ${Molecule1} ${NLevels1} 0 ${MinProcessInProc} ${Molecule2} ${NLevels2} 0 ${MaxProcessInProc} ${RmTrajFlg} ${BinaryTrajFlg} &

    done
    wait

  else
    echo "  [PostTrajectoriesAtNode]: -> Postprocessing Levels from "${MinLevel1}" to "${MaxLevel1}" for Molecule 1 and from "${MinLevel2}" to "${MaxLevel2}" for Molecule 2"

    MinLevel1Temp=${MinLevel1}
    MaxLevel1Temp=${MaxLevel1}
    MinLevel2Temp=${MinLevel2}
    MaxLevel2Temp=${MaxLevel2}
    NProcessesTot=0
    for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
      if [ ${iLevel1} -ge ${MinLevel1Temp} ] && [ ${iLevel1} -le ${MaxLevel1Temp} ]; then
        if [ ${SymmFlg} -eq 1 ]; then
          MinLevel2Temp=${iLevel1}
        fi
        for (( iLevel2=1; iLevel2<=${NLevels2}; iLevel2++ )); do
          if [ ${iLevel2} -ge ${MinLevel2Temp} ] && [ ${iLevel2} -le ${MaxLevel2Temp} ]; then
            NProcessesTot=$((NProcessesTot+1))
          fi
        done
      fi
    done  
    echo "  [PostTrajectoriesAtNode]: -> Total Nb of Processes to Run in This Node = "${NProcessesTot}


    NProcessesPerProc="$(bc <<< "scale = 10; ${NProcessesTot} / ${NProc}")"
    NProcessesPerProc="$(echo ${NProcessesPerProc} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
    echo "  [PostTrajectoriesAtNode]: -> Nb of Processes per Processor = "${NProcessesPerProc}
    echo "  [PostTrajectoriesAtNode] "


    MinLevel1Temp=${MinLevel1}
    MaxLevel1Temp=${MaxLevel1}
    MinLevel2Temp=${MinLevel2}
    MaxLevel2Temp=${MaxLevel2}
    MinLevel1Proc=${MinLevel1Temp}
    MinLevel2Proc=${MinLevel2Temp}
    iProcessesTot=0
    iProcessesPerProc=0
    iProc=0
    for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
      if [ ${iLevel1} -ge ${MinLevel1Temp} ] && [ ${iLevel1} -le ${MaxLevel1Temp} ]; then
        if [ ${SymmFlg} -eq 1 ]; then
          MinLevel2Temp=${iLevel1}
        fi
        for (( iLevel2=1; iLevel2<=${NLevels2}; iLevel2++ )); do
          if [ ${iLevel2} -ge ${MinLevel2Temp} ] && [ ${iLevel2} -le ${MaxLevel2Temp} ]; then
            iProcessesTot=$((iProcessesTot+1))
            iProcessesPerProc=$((iProcessesPerProc+1))

            if [ ${iProcessesPerProc} -eq ${NProcessesPerProc} ] || [ ${iProcessesTot} -eq ${NProcessesTot} ]; then 
              iProc=$((iProc+1))
              MaxLevel1Proc=${iLevel1}
              MaxLevel2Proc=${iLevel2}
              
              echo "  [PostTrajectoriesAtNode]: Postprocessing Trajectories @ iProc = "${iProc}"; MinLevel1Proc = "${MinLevel1Proc}" and MaxLevel1Proc = "${MaxLevel1Proc}"; MinLevel2Proc = "${MinLevel2Proc}" and MaxLevel2Proc = "${MaxLevel2Proc}

              bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_WORKING_DIR} ${COARSEAIR_OUTPUT_DIR} ${COARSEAIR_SH_DIR} ${System} ${StochPESFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${NNode} ${iNode} ${NProc} ${iProc} ${NMolecules} ${SymmFlg} ${Molecule1} ${NLevels1} ${MinLevel1Proc} ${MaxLevel1Proc} ${Molecule2} ${NLevels2} ${MinLevel2Proc} ${MaxLevel2Proc} ${RmTrajFlg} ${BinaryTrajFlg} &

              iProcessesPerProc=0
              MinLevel1Proc=${MaxLevel1Proc}
              if [ ${MaxLevel2Proc} -ge ${MaxLevel2} ]; then
                MinLevel1Proc=$((MaxLevel1Proc+1))
                if [ ${SymmFlg} -eq 1 ]; then
                  MinLevel2Proc=${MinLevel1Proc}
                else
                  MinLevel2Proc=${MinLevel2}
                fi
              else
                MinLevel2Proc=$((MaxLevel2Proc+1))
              fi
            fi

          fi
        done
      fi
    done
    wait

  fi

  end=`date +%s`
  runtime=$((end-start))
  echo "  [PostTrajectoriesAtNode]: Done Postprocessing Trajectories. RunTime = "${runtime}"s"      
      
}
#================================================================================================================================#


# ------------------------------------------------------------------------------------------------------------- PostTrajectories #
function PostTrajectories {

  echo "      [PostTrajectories]: COARSEAIR_OUTPUT_DIR     = "${COARSEAIR_OUTPUT_DIR}
  echo "      [PostTrajectories]: COARSEAIR_BIN_OUTPUT_DIR = "${COARSEAIR_BIN_OUTPUT_DIR}
  echo "      [PostTrajectories]: COARSEAIR_SH_DIR         = "${COARSEAIR_SH_DIR}
  echo "      [PostTrajectories]: System                   = "${System}
  echo "      [PostTrajectories]: TranFlg                  = "${TranFlg}
  echo "      [PostTrajectories]: Tran                     = "${Tran}
  echo "      [PostTrajectories]: Tint                     = "${Tint}
  echo "      [PostTrajectories]: StochPESFlg              = "${StochPESFlg}
  echo "      [PostTrajectories]: NPESs                    = "${NPESs}
  echo "      [PostTrajectories]: iNode                    = "${iNode}
  echo "      [PostTrajectories]: iProc                    = "${iProc}
  echo "      [PostTrajectories]: iPESStart                = "${iPESStart}
  echo "      [PostTrajectories]: Molecule1                = "${Molecule1}
  echo "      [PostTrajectories]: iLevel1                  = "${iLevel1}
  echo "      [PostTrajectories]: Molecule2                = "${Molecule2}
  echo "      [PostTrajectories]: iLevel2                  = "${iLevel2}
  echo "      [PostTrajectories]: Velocity                 = "${Velocity}
  echo "      [PostTrajectories]: BinaryTrajFlg            = "${BinaryTrajFlg}

  echo "      [PostTrajectories]: -> Postprocessing Molecule 1, Level "${iLevel1}" and Molecule 2, Level "${iLevel2}


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
      echo "      [PostTrajectories]: --- PES "${iPES}" ----------------------- "
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
      echo "      [PostTrajectories]: -> Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}": Tot Nb of Trajectories for PES "${iPES}": "${NTraj}
      
      if [ ${NTraj} -gt 0 ] || [ ${BinaryTrajFlg} -gt 0 ]; then
        
        echo "      [PostTrajectories]: Calling TrajectoriesStats"
        TrajectoriesStats
        
      fi
      
      echo "      [PostTrajectories]: Calling FromCrossToRates"
      FromCrossToRates  
      wait

      if [ ${StochPESFlg} -eq 1 ]; then 
        NLinesTry=$(wc -l < ${COARSEAIR_OUTPUT_DIR}"/"${System}"/"${Molecule1}"/Rates/T_"${Tran%.*}"_"${Tint%.*}"/Bin"$iLevel1".csv."$iPES)
        iTry=$((iTry+1))
      else
        iTry=4
      fi
    done
      
    if [ ${StochPESFlg} -eq 1 ]; then 
      rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv
      echo "      [PostTrajectories]: --- iPES "${iPES} " --------------- DONE -- "
    fi
    iPES=$((iPES+1))
  done   


  if [ ${StochPESFlg} -eq 1 ]; then
    #scp ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv.Orig ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv
    rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/NConvTraj.dat
    rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/statistics*
  fi

}
#================================================================================================================================#


# ------------------------------------------------------------------------------------------------------------ TrajectoriesStats #
function TrajectoriesStats {

  echo "        [TrajectoriesStats]: COARSEAIR_BIN_OUTPUT_DIR = "${COARSEAIR_BIN_OUTPUT_DIR}
  echo "        [TrajectoriesStats]: Tran                     = "${Tran}
  echo "        [TrajectoriesStats]: Tint                     = "${Tint}
  echo "        [TrajectoriesStats]: BinaryTrajFlg            = "${BinaryTrajFlg}

  TrajectoriesStatsCommand="coarseair-trajectoriesstats.x"

  cd ${COARSEAIR_BIN_OUTPUT_DIR}
  
  echo "        [TrajectoriesStats]: Computing Statistics for Trajectories. Command: "${TrajectoriesStatsCommand}
  eval ${TrajectoriesStatsCommand} ${Tran} ${Tint} ${BinaryTrajFlg}
  echo "        [TrajectoriesStats]: Done with TrajectoriesStats"
  
}
#================================================================================================================================#


# ------------------------------------------------------------------------------------------------------------- FromCrossToRates #
function FromCrossToRates {

  echo "        [FromCrossToRates]: COARSEAIR_BIN_OUTPUT_DIR = "${COARSEAIR_BIN_OUTPUT_DIR}
  echo "        [FromCrossToRates]: Tran                     = "${Tran}
  echo "        [FromCrossToRates]: Tint                     = "${Tint}
  echo "        [FromCrossToRates]: iLevel1                  = "${iLevel1}
  echo "        [FromCrossToRates]: iLevel2                  = "${iLevel2}
  echo "        [FromCrossToRates]: iPES                     = "${iPES}
  echo "        [FromCrossToRates]: iNode                    = "${iNode}
  echo "        [FromCrossToRates]: iProc                    = "${iProc}
  echo "        [FromCrossToRates]: TranFlg                  = "${TranFlg}
  echo "        [FromCrossToRates]: Velocity                 = "${Velocity}
  echo "        [FromCrossToRates]: BinaryTrajFlg            = "${BinaryTrajFlg}


  PostTrajectoriesCommand="coarseair-posttrajectories.x"
  TrajErrorFile=${COARSEAIR_OUTPUT_DIR}/"RatesErrors_Node"${iNode}"_Proc"${iProc}".dat"

  echo "        [FromCrossToRates]: Postprocessing Trajectories for iLevel1 = "${iLevel1}" and for iLevel2 = "${iLevel2}", @ iProc = "${iProc}". Command: "${PostTrajectoriesCommand} ${Tran} 
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
        echo "        [FromCrossToRates]: No Trajectories for Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}
        #rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/*
      fi
      
    else
    
     if [ -e "$TrajErrorFile" ]; then                                                                                   
        echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
      else
        echo "# List of Levels / Bins that Generated Errors during Rates Computation" > $TrajErrorFile
        echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
      fi
      echo "        [FromCrossToRates]: No Trajectories for Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}
      
    fi                                                                                                                  
    
  else
  
    if [ -e "$TrajErrorFile" ]; then
      echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
    else
      echo "# List of Levels / Bins that Generated Errors during Rates Computation" > $TrajErrorFile
      echo ${iLevel1}","${iLevel2} >> $TrajErrorFile
    fi
    echo "        [FromCrossToRates]: No Trajectories for Molecule 1, Level/Bin "${iLevel1}"; Molecule 2, Level/Bin "${iLevel2}
    
  fi
  
  end=`date +%s`
  runtime=$((end-start))
  echo "        [FromCrossToRates]: Done with FromCrossToRates for iLevel1 = "${iLevel1}" and for iLevel2 = "${iLevel2}", @ iProc = "${iProc}". RunTime = "${runtime}"s"

  
}
#================================================================================================================================#


# ----------------------------------------------------------------------------------------------------------------- ComputeRates #
function MergeAllRates {

  for Tran in "${Tran_vec[@]}"; do :
    echo '  [MergeAllRates] '
    echo "  [MergeAllRates]: ===   Translational Temperature      = " ${Tran} " ============================ "
    #for Tint in "${Tint_Vec[@]}"; do :
      Tint=${Tran}
      echo "  [MergeAllRates]: ===== Internal Temperature           = " ${Tint} " ======================== "
      echo '  [MergeAllRates]'     

      MinLevel1Temp=${MinLevel1}
      MaxLevel1Temp=${MaxLevel1}
      MinLevel2Temp=${MinLevel2}
      MaxLevel2Temp=${MaxLevel2}
      for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
        if [ ${iLevel1} -ge ${MinLevel1Temp} ] && [ ${iLevel1} -le ${MaxLevel1Temp} ]; then
          if [ ${SymmFlg} -eq 1 ]; then
            MinLevel2Temp=${iLevel1}
          fi
          for (( iLevel2=1; iLevel2<=${NLevels2}; iLevel2++ )); do
            if [ ${iLevel2} -ge ${MinLevel2Temp} ] && [ ${iLevel2} -le ${MaxLevel2Temp} ]; then
          

              if [ ${TranFlg} -eq 0 ]; then 
                COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
              else
                COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
              fi
              cd ${COARSEAIR_BIN_OUTPUT_DIR}


              if [ -f ../trajectories-Tot.csv ]; then  
                tail -n+2 ./trajectories.csv >> ../trajectories-Tot.csv
                if [ -f ./PaQSOl.out ]; then
                  echo "### Molecule 1, Level / Bin Nb "$iLevel1"; Molecule 2, Level / Bin Nb "$iLevel2 >> ../PaQSOl-Tot.out
                  tail -n+2 ./PaQSOl.out >> ../PaQSOl-Tot.out
                  rm -rf ./PaQSOl.out
                fi
              else
                cat ./trajectories.csv > ../trajectories-Tot.csv
                if [ -f ./PaQSOl.out ]; then
                  echo "### Molecule 1, Level / Bin Nb "$iLevel1"; Molecule 2, Level / Bin Nb "$iLevel2 >> ../PaQSOl-Tot.out
                  tail -n+2 ./PaQSOl.out >> ../PaQSOl-Tot.out
                  rm -rf ./PaQSOl.out
                fi
              fi
              
              rm -rf ./trajectories.csv

            fi
          done
        fi
      done

    #done
  echo "  [MergeAllRates]: ------------------------------------------------------------ "
  done
  wait

}
#================================================================================================================================#
