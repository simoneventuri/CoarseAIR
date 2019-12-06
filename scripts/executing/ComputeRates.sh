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
  echo "  [ComputeTrajsPBS]: MinLevel1            = "${MinLevel1}
  echo "  [ComputeTrajsPBS]: MaxLevel1            = "${MaxLevel1}
  
  
  cd ${COARSEAIR_OUTPUT_DIR}/..
  if [ "${ProcType}" = "none" ]; then
    scp ${COARSEAIR_SH_DIR}/RunTrajectories-Format-UIUC.pbs     ./
  else
    scp ${COARSEAIR_SH_DIR}/RunTrajectories-Format-Pleiades.pbs ./
  fi
  
  NLevelsTemp1=$((MaxLevel1 - $((MinLevel1 - 1)) ))
  echo "  [ComputeTrajsPBS]: -> Computing Trajectories for "${NLevelsTemp1}" Levels/Bins"
  
  NLevelsPerNode="$(bc <<< "scale = 10; ${NLevelsTemp1} / ${NNode}")"
  NLevelsPerNode="$(echo ${NLevelsPerNode} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
  echo "  [ComputeTrajsPBS]: -> Nb Levels/Bins per Node = "${NLevelsPerNode}

  iNode=1
  iLevelMax=$((MinLevel1 - 1))
  while [ ${iNode} -le ${NNode} ]; do 
    
    if [ ${MinLevel1} -eq 0 ]; then
      iLevelMin=0
      iLevelMax=${iNode}
    else 
      iLevelMin=$((iLevelMax + 1))
      iLevelMax=$((iLevelMin + ${NLevelsPerNode} - 1))
      if [ ${iLevelMax} -gt ${MaxLevel1} ]; then
        iLevelMax=${MaxLevel1}
      fi
    fi
    echo "  [ComputeTrajsPBS]: -> iNode = "${iNode}"; iLevelMin = "${iLevelMin}"; iLevelMax = "${iLevelMax}
    
    if [ "${ProcType}" = "none" ]; then
      sed -e '3s/$/1:ppn='${NProc}'/'                       'RunTrajectories-Format-UIUC.pbs'     > 'RunTrajectoriesTEMP-1.pbs'
    else
      sed -e '3s/$/1:ncpus='${NProc}':model='${ProcType}'/' 'RunTrajectories-Format-Pleiades.pbs' > 'RunTrajectoriesTEMP-1.pbs'
    fi
    sed -e '12s/$/'${NProc}'/'                              'RunTrajectoriesTEMP-1.pbs'  > 'RunTrajectoriesTEMP-2.pbs'
    
    sed -e '4s/$/'${iLevelMin}'/'                           'RunTrajectoriesTEMP-2.pbs'  > 'RunTrajectoriesTEMP-3.pbs'
    sed -e '6s/$/'${iLevelMin}'/'                           'RunTrajectoriesTEMP-3.pbs'  > 'RunTrajectoriesTEMP-4.pbs'
    sed -e '147s/$/'${iLevelMin}'/'                         'RunTrajectoriesTEMP-4.pbs'  > 'RunTrajectoriesTEMP-5.pbs'
    sed -e '148s/$/'${iLevelMax}'/'                         'RunTrajectoriesTEMP-5.pbs'  > 'RunTrajectoriesTEMP-6.pbs' 
    
    sed -e '149s/$/'${Tran}'/'                              'RunTrajectoriesTEMP-6.pbs'  > 'RunTrajectoriesTEMP-7.pbs' 
    sed -e '150s/$/'${Tint}'/'                              'RunTrajectoriesTEMP-7.pbs'  > 'RunTrajectoriesTEMP-8.pbs'

    sed -e '151s/$/'${iNode}'/'                             'RunTrajectoriesTEMP-8.pbs'  > 'RunTrajectories-'${iLevelMin}'.pbs'

    qsub ./'RunTrajectories-'${iLevelMin}'.pbs'
    
    rm -rf ./RunTrajectoriesTEMP*

    iNode=$((iNode+1))
  done 
  
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
  echo "  [ComputeTrajs]: NLevels1              = "${NLevels1}
  echo "  [ComputeTrajs]: MinLevel1             = "${MinLevel1}
  echo "  [ComputeTrajs]: MaxLevel1             = "${MaxLevel1}
  echo "  [ComputeTrajs]: NLevels2              = "${NLevels2}
  echo "  [ComputeTrajs]: MinLevel2             = "${MinLevel2}
  echo "  [ComputeTrajs]: MaxLevel2             = "${MaxLevel2}
  echo "  [ComputeTrajs]: StochPESFlg           = "${StochPESFlg}
  echo "  [ComputeTrajs]: NPESs                 = "${NPESs}
  echo "  [ComputeTrajs]: iPESStart             = "${iPESStart}
  
  if [ ${MinLevel1} -eq 0 ]; then 
    echo "  [ComputeTrajs]: Reading Levels/Bins from File "${COARSEAIR_INPUT_DIR}/LevelsToRunList.inp
    
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
    echo "  [ComputeTrajs]: For Node "${iNode}", the first Level to read from file is the "${MinLevelInNode}"-th in the List"
    echo "  [ComputeTrajs]: For Node "${iNode}", the last  Level to read from file is the "${MaxLevelInNode}"-th in the List"
    
    iCount=0
    while IFS= read -r iLevel1
      do
      iCount=$((iCount+1))
      
      if [ ${iCount} -ge ${MinLevelInNode} ] && [ ${iCount} -le ${MaxLevelInNode} ]; then
      
        echo "  [ComputeTrajs]: --- Molecule 1, Level/Bin " ${iLevel1} " -------------------------------- "

        #iLevel2=0  
        iLevel2=${MinLevel2}
        while [ ${iLevel2} -le ${MaxLevel2} ]; do
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
          iLevel2=$((iLevel2+1))
        done
      
      fi
      
      echo "  [ComputeTrajs]: --- Molecule 1, Level/Bin " ${iLevel1} " ------------------------ DONE -- "
    done < "${COARSEAIR_INPUT_DIR}/LevelsToRunList.inp"
  
  else
  
    iLevel1=${MinLevel1}
    while [ ${iLevel1} -le ${MaxLevel1} ]
      do
      echo "  [ComputeTrajs]: --- Molecule 1, Level/Bin " ${iLevel1} " -------------------------------- "

      #iLevel2=0  
      iLevel2=${MinLevel2}
      while [ ${iLevel2} -le ${MaxLevel2} ]; do
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
        iLevel2=$((iLevel2+1))
      done
      
      echo "  [ComputeTrajs]: --- Molecule 1, Level/Bin " ${iLevel1} " ------------------------ DONE -- "
      iLevel1=$((iLevel1+1))
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
  echo "  [SplitTrajsPESs]: MinLevel1                = "${MinLevel1}
  echo "  [SplitTrajsPESs]: MaxLevel1                = "${MaxLevel1}
  echo "  [SplitTrajsPESs]: MinLevel2                = "${MinLevel2}
  echo "  [SplitTrajsPESs]: MaxLevel2                = "${MaxLevel2}

  echo "  [SplitTrajsPESs]: -> Splitting PESs in the Trajectory Files"

  iLevel1=${MinLevel1}
  while [ $iLevel1 -le ${MaxLevel1} ]; do  

    iLevel2=${MinLevel2}
    while [ $iLevel2 -le ${MaxLevel2} ]; do  

      echo "  [SplitTrajsPESs]: Splitting Trajectories for iLevel1 = "${iLevel1}" and iLevel2 = "${iLevel2}

      if [ ${TranFlg} -eq 0 ]; then 
        COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
      else
        COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}
      fi

      python3 ${COARSEAIR_SH_DIR}/SplitTrajsPESs.py ${COARSEAIR_BIN_OUTPUT_DIR} ${iPESStart}

      iLevel2=$((iLevel2+1))
    done
    
    iLevel1=$((iLevel1+1))
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
  echo "  [PostTrajectoriesPBS]: MinLevel1            = "${MinLevel1}
  echo "  [PostTrajectoriesPBS]: MaxLevel1            = "${MaxLevel1}
  echo "  [PostTrajectoriesPBS]: BinaryTrajFlg        = "${BinaryTrajFlg}

  
  cd ${COARSEAIR_OUTPUT_DIR}/..
  if [ "${ProcType}" = "none" ]; then
    scp ${COARSEAIR_SH_DIR}/PostTrajectories-Format-UIUC.pbs     ./
  else
    scp ${COARSEAIR_SH_DIR}/PostTrajectories-Format-Pleiades.pbs ./
  fi
  
  NLevelsTemp1=$((MaxLevel1 - $((MinLevel1 - 1)) ))
  echo "  [PostTrajectoriesPBS]: -> Postprocessing Trajectories for "${NLevelsTemp1}" Levels/Bins"
  
  NLevelsPerNode="$(bc <<< "scale = 10; ${NLevelsTemp1} / ${NNode}")"
  NLevelsPerNode="$(echo ${NLevelsPerNode} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
  echo "  [PostTrajectoriesPBS]: -> Nb Levels/Bins per Node = "${NLevelsPerNode}

  iNode=1
  iLevelMax=$((MinLevel1 - 1))
  while [ ${iNode} -le ${NNode} ]; do 
    
    if [ ${MinLevel1} -eq 0 ]; then
      iLevelMin=0
      iLevelMax=${iNode}
    else 
      iLevelMin=$((iLevelMax + 1))
      iLevelMax=$((iLevelMin + ${NLevelsPerNode} - 1))
      if [ ${iLevelMax} -gt ${MaxLevel1} ]; then
        iLevelMax=${MaxLevel1}
      fi
    fi
    echo "  [PostTrajectoriesPBS]: -> iNode = "${iNode}"; iLevelMin = "${iLevelMin}"; iLevelMax = "${iLevelMax}
    
    if [ "${ProcType}" = "none" ]; then
      sed -e '3s/$/1:ppn='${NProc}'/'                       'PostTrajectories-Format-UIUC.pbs'     > 'PostTrajectoriesTEMP-1.pbs'
    else
      sed -e '3s/$/1:ncpus='${NProc}':model='${ProcType}'/' 'PostTrajectories-Format-Pleiades.pbs' > 'PostTrajectoriesTEMP-1.pbs'
    fi
    sed -e '12s/$/'${NProc}'/'                              'PostTrajectoriesTEMP-1.pbs'  > 'PostTrajectoriesTEMP-2.pbs'
    
    sed -e '4s/$/'${iLevelMin}'/'                           'PostTrajectoriesTEMP-2.pbs'  > 'PostTrajectoriesTEMP-3.pbs'
    sed -e '6s/$/'${iLevelMin}'/'                           'PostTrajectoriesTEMP-3.pbs'  > 'PostTrajectoriesTEMP-4.pbs'
    sed -e '147s/$/'${iLevelMin}'/'                         'PostTrajectoriesTEMP-4.pbs'  > 'PostTrajectoriesTEMP-5.pbs'
    sed -e '148s/$/'${iLevelMax}'/'                         'PostTrajectoriesTEMP-5.pbs'  > 'PostTrajectoriesTEMP-6.pbs' 
    
    sed -e '149s/$/'${Tran}'/'                              'PostTrajectoriesTEMP-6.pbs'  > 'PostTrajectoriesTEMP-7.pbs' 
    sed -e '150s/$/'${Tint}'/'                              'PostTrajectoriesTEMP-7.pbs'  > 'PostTrajectoriesTEMP-8.pbs'

    sed -e '151s/$/'${iNode}'/'                             'PostTrajectoriesTEMP-8.pbs'  > 'PostTrajectories-'${iLevelMin}'.pbs'

    qsub ./'PostTrajectories-'${iLevelMin}'.pbs'
    
    rm -rf ./PostTrajectoriesTEMP*

    iNode=$((iNode+1))
  done 
  
  #rm -rf ./PostTrajectories-*
  
  echo "  [PostTrajectoriesPBS]: Done with Submitting PBS Files for Postprocessing Trajectories. Now we have to wait for ALL NODES being done."
      
}
#================================================================================================================================#


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

  echo "  [PostTrajectoriesAtNode]: -> Postprocessing Levels from "${MinLevel1}" to "${MaxLevel1}" for Molecule 1 and from "${MinLevel2}" to "${MaxLevel2}" for Molecule 2"

  VelocityFile=${COARSEAIR_OUTPUT_DIR}/Velocity_${Tran}.dat
  if [ -f $exist ]; then
    Velocity=$(sed '2q;d' ${VelocityFile})
    echo "  [PostTrajectoriesAtNode]: Velocity = "${Velocity}
  else
    echo "  [PostTrajectoriesAtNode]: ERROR! Velocity File does exist! CoarseAIR cannot compute Cross Sections!"
    exit1
  fi

  
  NLevelsTemp1=$((MaxLevel1-MinLevel1+1))
  NLevelsPerProc1="$(bc <<< "scale = 10; ${NLevelsTemp1} / ${NProc}")"
  NLevelsPerProc1="$(echo ${NLevelsPerProc1} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
  echo "  [PostTrajectoriesAtNode]: -> Molecule 1, Nb Levels/Bins per Processor = "${NLevelsPerProc1}

  NLevelsTemp2=$((MaxLevel2-MinLevel2+1))
  NLevelsPerProc2="$(bc <<< "scale = 10; ${NLevelsTemp2} / ${NProc}")"
  NLevelsPerProc2="$(echo ${NLevelsPerProc2} | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}')"
  echo "  [PostTrajectoriesAtNode]: -> Molecule 2, Nb Levels/Bins per Processor = "${NLevelsPerProc2}


  start=`date +%s`
 
  iProc=1
  MaxLevel1Tot=$MaxLevel1
  MaxLevel2Tot=$MaxLevel2
  while [ $iProc -le ${NProc} ] && [ $MinLevel1 -le ${MaxLevel1Tot} ]; do  
    MaxLevel1=$((MinLevel1+NLevelsPerProc1-1)) 
    MaxLevel2=$((MinLevel2+NLevelsPerProc2-1))
    #MaxLevel2=0

    echo "  [PostTrajectoriesAtNode]: Postprocessing Trajectories @ iProc = "${iProc}"; MinLevel1 = "${MinLevel1}" and MaxLevel1 = "${MaxLevel1}"; MinLevel2 = "${MinLevel2}" and MaxLevel2 = "${MaxLevel2}
    if [ ${MaxLevel1} -ge ${MaxLevel1Tot} ]; then
      MaxLevel1=${MaxLevel1Tot}
    fi
    if [ ${MaxLevel2} -ge ${MaxLevel2Tot} ]; then
      MaxLevel2=${MaxLevel2Tot}
    fi

    bash ${COARSEAIR_SH_DIR}/PostTrajectoriesAtProc.sh ${COARSEAIR_WORKING_DIR} ${COARSEAIR_OUTPUT_DIR} ${COARSEAIR_SH_DIR} ${System} ${StochPESFlg} ${NPESs} ${iPESStart} ${TranFlg} ${Tran} ${Tint} ${Velocity} ${NNode} ${iNode} ${NProc} ${iProc} ${Molecule1} ${NLevels1} ${MinLevel1} ${MaxLevel1} ${Molecule2} ${NLevels2} ${MinLevel2} ${MaxLevel2} ${RmTrajFlg} ${BinaryTrajFlg} &
    
    MinLevel1=$((MaxLevel1+1))
    #MinLevel2=$((MaxLevel2+1)) 
    MinLevel2=0
    iProc=$((iProc+1))
  done
  wait

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
      if [ -f ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv ]; then
        NTraj=$(wc -l < ${COARSEAIR_BIN_OUTPUT_DIR}/trajectories.csv)
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

      iLevel1=${MinLevel1}
      while [ ${iLevel1} -le ${MaxLevel1} ]
        do
        #echo "  [MergeAllRates]: --- Molecule 1, Level/Bin " ${iLevel1} " ----------------------------- "

        #iLevel2=0  
        iLevel2=${MinLevel2}
        while [ ${iLevel2} -le ${MaxLevel2} ]; do
          #echo "  [ComputeRates]: ----- Molecule 2, Level/Bin = " ${iLevel2} " --------------------- "
          #echo "  [ComputeRates]"
          
          
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
        
          #echo "  [ComputeRates]: ---------------------------------------------------------- "
          iLevel2=$((iLevel2+1))
        done
        
        #echo "  [MergeAllRates]: ------------------------------------------------------------ "
        #echo " "
        iLevel1=$((iLevel1+1))
      done         
  
  echo "  [MergeAllRates]: ------------------------------------------------------------ "
  done
  wait

}
#================================================================================================================================#
