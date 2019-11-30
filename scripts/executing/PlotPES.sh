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

# ---------------------------------------------------------------------------------------------------------------------- PlotPES #
function PlotPES {

  PlotPESCommand="coarseair-plotpes.x"

  # ---  PATHS ------------------------------------------------------------------------------------ #
  export RUN_DIR=$(pwd)
  export COARSEAIR_WORKING_DIR=$(pwd)
  export COARSEAIR_INPUT_DIR=$(pwd)/input
  export COARSEAIR_INPUT_FILE=$(pwd)/input/CoarseAIR.inp
  export COARSEAIR_OUTPUT_DIR=$(pwd)/Test
  export COARSEAIR_BIN_OUTPUT_DIR=$(pwd)/Test
  #-------------------------------------------------------------------------------------------------#

  startTot=`date +%s`
  
  echo "  [PlotPES]: COARSEAIR_WORKING_DIR = "${COARSEAIR_WORKING_DIR}
  echo "  [PlotPES]: COARSEAIR_OUTPUT_DIR  = "${COARSEAIR_OUTPUT_DIR}
  echo "  [PlotPES]: COARSEAIR_INPUT_DIR   = "${COARSEAIR_INPUT_DIR}
  echo "  [PlotPES]"

  echo "  [PlotPES]: Computing and Plotting Potential Energy Surface. Command: ${PlotPESCommand}"
  cd ${COARSEAIR_OUTPUT_DIR}
  eval ${PlotPESCommand} ${COARSEAIR_OUTPUT_DIR}
  cd ..
  
  endTot=`date +%s`
  Totruntime=$((endTot-startTot))
  echo "  [PlotPES]: Done with PlotPES. Total RunTime = "${Totruntime}"s"
  echo " "

}
#================================================================================================================================#