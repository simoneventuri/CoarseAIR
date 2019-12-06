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

# ---  PARAMETERS ------------------------------------------------------------------------------- #
MinLevel1=1
MaxLevel1=100
#=================================================================================================#


mkdir '../../../bnn_results/T_10000_10000'
iLevels1=${MinLevel1}
while [ ${iLevels1} -le ${MaxLevel1} ]
  do
  echo "  [CopyResults]: --- Molecule 1, Level/Bin " ${iLevels1} " ----------------------------- "
  
  mkdir '../../../bnn_results/T_10000_10000/Bins_'${iLevels1}'_0'
  scp  './Bins_'${iLevels1}'_0/trajectories.csv' '../../../bnn_results/T_10000_10000/Bins_'${iLevels1}'_0'
  
  echo "  [CopyResults]: ------------------------------------------------------------ "
  echo " "
  iLevels1=$((iLevels1+1))
done  
