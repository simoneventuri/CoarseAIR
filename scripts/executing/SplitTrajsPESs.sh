#!/bin/bash
#===============================================================================================================

iPESStart=201

COARSEAIR_SH_DIR="/home1/sventuri/WORKSPACE/CoarseAIR/coarseair/scripts/executing/"

COARSEAIR_OUTPUT_DIR=$(pwd)/Test/

Tran=10000.0
Tint=10000.0

iLevel1=2
iLevel2=0
MaxLevel=200
while [ $iLevel1 -le ${MaxLevel1} ]; do  

  echo " Splitting Trajectories for iLvevel1 = "${iLevel1}

  COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevel1}"_"${iLevel2}

  python3 ${COARSEAIR_SH_DIR}/SplitTrajsPESs.py ${COARSEAIR_BIN_OUTPUT_DIR} ${iPESStart}
  
  iLevel1=$((iLevel1+1))
done

exit 0