##==============================================================================================================
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
##==============================================================================================================

# Ex.: python3 ../coarseair/scripts/executing/SelectLevelsToRun.py /home/venturi/WORKSPACE/CoarseAIR/N4_VS/Test/ 1001

import os
import sys
import numpy as np
import pandas

##==============================================================================================================
if (len(sys.argv) > 1):
    COARSEAIR_BIN_OUTPUT_DIR = sys.argv[1]
    print("[SelectLevelsToRun.py]: COARSEAIR_BIN_OUTPUT_DIR = ", COARSEAIR_BIN_OUTPUT_DIR)
    MinNTraj                 = int(sys.argv[2])
    print("[SelectLevelsToRun.py]: MinNTraj                 = ", MinNTraj)
else:
    print("[SelectLevelsToRun.py]: ERROR! This Program requires TWO ARGUMENTS: COARSEAIR_BIN_OUTPUT_DIR and MinNTraj" )
##--------------------------------------------------------------------------------------------------------------

PathToFile = COARSEAIR_BIN_OUTPUT_DIR + '/Overall_NConvTraj.csv'
print('[SelectLevelsToRun.py]:   Reading File: ', PathToFile)

Data = pandas.read_csv(PathToFile, header=None, skiprows=1)
Data = Data.apply(pandas.to_numeric, errors='coerce')

iLevelVec = np.array(Data.values[:,0],  dtype=np.int64  )
jLevelVec = np.array(Data.values[:,1],  dtype=np.int64  )
NTrajVec  = np.array(Data.values[:,2],  dtype=np.int64  )


LevelsFile = COARSEAIR_BIN_OUTPUT_DIR + '/ProcessesToRunList.inp'
print('[SelectLevelsToRun.py]: Writing List of Levels to Run in: ' + LevelsFile )
csvlevels  = open(LevelsFile, 'w')
iProc = -1
for NTraj in NTrajVec:
    iProc  = iProc + 1
    iLevel = iLevelVec[iProc]
    jLevel = jLevelVec[iProc]
    if (NTraj < MinNTraj):
        print('[SelectLevelsToRun.py]: Adding (iLevel, jLevel) = (' + str(iLevel) + ',' + str(jLevel) + '). It has only ' + str(NTraj) + ' Trajectories Converged.' )
        Line = '%d-%d\n' % ((iLevel, jLevel))
        csvlevels.write(Line)