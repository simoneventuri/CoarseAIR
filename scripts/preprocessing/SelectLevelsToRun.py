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

# Ex.: python3 ../coarseair/scripts/executing/SelectLevelsToRun.py /home/venturi/WORKSPACE/CoarseAIR/N4_VS/Test/ 1001 1

import os
import sys
import numpy as np
import pandas

##==============================================================================================================
if (len(sys.argv) > 1):
    COARSEAIR_OUTPUT_DIR = sys.argv[1]
    print("[SelectLevelsToRun.py]: COARSEAIR_OUTPUT_DIR = ", COARSEAIR_OUTPUT_DIR)
    MinNQoI              = int(sys.argv[2])
    print("[SelectLevelsToRun.py]: MinNQoI              = ", MinNQoI)
    QoIFlg               = int(sys.argv[3])
    print("[SelectLevelsToRun.py]: QoIFlg               = ", QoIFlg)
else:
    print("[SelectLevelsToRun.py]: ERROR! This Program requires 3 ARGUMENTS: COARSEAIR_OUTPUT_DIR;   MinNQoI;   QoIFlg" )
##--------------------------------------------------------------------------------------------------------------

if (QoIFlg == 1):
    print('[SelectLevelsToRun.py]:   QoI: Number of Trajectories')
    PathToFile = COARSEAIR_OUTPUT_DIR + '/Overall_NConvTraj.csv'
elif (QoIFlg == 2):
    print('[SelectLevelsToRun.py]:   QoI: Number of Rates per File')
    PathToFile = COARSEAIR_OUTPUT_DIR + '/Overall_NRates.csv'
print('[SelectLevelsToRun.py]:     Reading File: ', PathToFile)

Data = pandas.read_csv(PathToFile, header=None, skiprows=1)
Data = Data.apply(pandas.to_numeric, errors='coerce')

iLevelVec = np.array(Data.values[:,0],  dtype=np.int64  )
jLevelVec = np.array(Data.values[:,1],  dtype=np.int64  )
NQoIVec   = np.array(Data.values[:,2],  dtype=np.int64  )


LevelsFile = COARSEAIR_OUTPUT_DIR + '/ProcessesToRunList.inp'
print('[SelectLevelsToRun.py]: Writing List of Levels to Run in: ' + LevelsFile )
csvlevels  = open(LevelsFile, 'w')
iProc = -1
for NQoI in NQoIVec:
    iProc  = iProc + 1
    iLevel = iLevelVec[iProc]
    jLevel = jLevelVec[iProc]
    if (NQoI < MinNQoI):
        print('[SelectLevelsToRun.py]: Adding (iLevel, jLevel) = (' + str(iLevel) + ',' + str(jLevel) + '). It has only ' + str(NQoI) + ' QoIs.' )
        Line = '%d-%d\n' % ((iLevel, jLevel))
        csvlevels.write(Line)