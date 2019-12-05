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

# --------------------------------------------------------------------------------------------------------------- GenerateLevels #
function GenerateLevels {
  
  echo "  [GenerateLevels]: COARSEAIR_WORKING_DIR = "${COARSEAIR_WORKING_DIR}
  echo "  [GenerateLevels]: COARSEAIR_OUTPUT_DIR  = "${COARSEAIR_OUTPUT_DIR}
  echo "  [GenerateLevels]: COARSEAIR_SH_DIR      = "${COARSEAIR_SH_DIR}
  echo "  [GenerateLevels]: NNode                 = "${NNode}
  echo "  [GenerateLevels]: NProc                 = "${NProc}
  echo "  [GenerateLevels]: System                = "${System}
  
  echo "  [GenerateLevels]: Calling ComputeLevelsAtNode"        
  ComputeLevelsAtNode
  
  echo "  [GenerateLevels]: Calling MergeLevels"        
  MergeLevels
  
}
#================================================================================================================================#


# ---------------------------------------------------------------------------------------------------------- ComputeLevelsAtNode #
function ComputeLevelsAtNode {

  iNode=1
  
  start=`date +%s`
  echo "    [ComputeLevelsAtNode]: Generating Energy Levels. Command: "${CompLevCommand}
  echo "    [ComputeLevelsAtNode]: Node "${iNode}
  
  if [[ ${NProc} -eq 1 ]]; then
    ComputeLevelsAtProc 1
  elif [[ ${NProc} -eq 2 ]]; then
	  parallel --xapply -j 2 "sh ${COARSEAIR_SH_DIR}/ComputeLevelsAtProc.sh '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_BIN_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} " ::: {1..2}
  elif [[ ${NProc} -eq 4 ]]; then
	  parallel --xapply -j 4 "sh ${COARSEAIR_SH_DIR}/ComputeLevelsAtProc.sh '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_BIN_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} " ::: {1..4}
  elif [[ ${NProc} -eq 8 ]]; then
	  parallel --xapply -j 8 "sh ${COARSEAIR_SH_DIR}/ComputeLevelsAtProc.sh '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_BIN_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} " ::: {1..8}
  elif [[ ${NProc} -eq 16 ]]; then 
    parallel --xapply -j 16 "sh ${COARSEAIR_SH_DIR}/ComputeLevelsAtProc.sh '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_BIN_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} " ::: {1..16}
  elif [[ ${NProc} -eq 20 ]]; then        
    parallel --xapply -j 32 "sh ${COARSEAIR_SH_DIR}/ComputeLevelsAtProc.sh '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_BIN_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} " ::: {1..20}
  elif [[ ${NProc} -eq 32 ]]; then  	    
    parallel --xapply -j 32 "sh ${COARSEAIR_SH_DIR}/ComputeLevelsAtProc.sh '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_BIN_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} " ::: {1..32}
  elif [[ ${NProc} -eq 40 ]]; then        
    parallel --xapply -j 32 "sh ${COARSEAIR_SH_DIR}/ComputeLevelsAtProc.sh '${COARSEAIR_OUTPUT_DIR}' '${COARSEAIR_BIN_OUTPUT_DIR}' '${COARSEAIR_WORKING_DIR}' ${NNode} ${iNode} ${NProc} {1} " ::: {1..40}
  else
    echo "ERROR: Number of Precessors not Coherent with ComputeLevelsAtNode (Check GenerateLevels.sh)"
    exit 1
  fi
  
  wait
  end=`date +%s`
  runtime=$((end-start))
  echo "    [ComputeLevelsAtNode]: Done with CompLevAtNodes. RunTime = "${runtime}"s"
}
#================================================================================================================================#


# ---------------------------------------------------------------------------------------------------------- ComputeLevelsAtProc #
function ComputeLevelsAtProc {

  iNode=1
  iProc=${1}
  
  echo "      [ComputeLevelsAtProc]: COARSEAIR_OUTPUT_DIR     = "${COARSEAIR_OUTPUT_DIR}
  echo "      [ComputeLevelsAtProc]: COARSEAIR_BIN_OUTPUT_DIR = "${COARSEAIR_BIN_OUTPUT_DIR}
  echo "      [ComputeLevelsAtProc]: NNode                    = "${NNode}
  echo "      [ComputeLevelsAtProc]: iNode                    = "${iNode}
  echo "      [ComputeLevelsAtProc]: NProc                    = "${NProc}
  echo "      [ComputeLevelsAtProc]: iProc                    = "${iProc}


  ComputeLevelsCommand="coarseair-computelevels.x"

  echo "      [ComputeLevelsAtProc]: Generating Energy Levels for Node "${iNode}", Processor "${iProc}

  cd ${COARSEAIR_OUTPUT_DIR}

  eval ${ComputeLevelsCommand} ${COARSEAIR_BIN_OUTPUT_DIR} ${NProc} ${iProc} ${NNode} ${iNode}

  echo "      [ComputeLevelsAtProc]: Done Generating Energy Levels for Node "${iNode}", Processor "${iProc}
}
#================================================================================================================================#


# ------------------------------------------------------------------------------------------------------------------ MergeLevels #
function MergeLevels {
  
  MergeLevelsCommand="coarseair-mergelevels.x"
    
  cd ${COARSEAIR_OUTPUT_DIR}

  echo "    [MergeLevels]: Merging Energy Levels. Command: ${MergeLevelsCommand}"
  eval ${MergeLevelsCommand} ${NNode} ${NProc}
  echo "    [MergeLevels]: Done with MergeLevels"
  
}
#================================================================================================================================#
