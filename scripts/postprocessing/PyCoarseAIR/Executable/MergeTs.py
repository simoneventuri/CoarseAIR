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
import os
import sys
import numpy as np

## Paths from where to Import Python Modules
SrcDir = '/Users/sventuri/WORKSPACE//CoarseAIR/coarseair' #os.environ['COARSEAIR_SOURCE_DIR'] 
sys.path.insert(1, SrcDir + '/scripts/postprocessing/PyCoarseAIR/Objects/')
sys.path.insert(2, SrcDir + '/scripts/postprocessing/PyCoarseAIR/ChemicalSystems')
sys.path.insert(2, SrcDir + '/scripts/postprocessing/PyCoarseAIR/Reading')
import ChemicalSystems
from System            import system
from Temperatures      import temperatures

from Reading           import Read_qnsEnBin, Read_RatesCGQCT
##--------------------------------------------------------------------------------------------------------------



################################################################################################################

Temp                     = temperatures()
## Vector of Temperatures
Temp.TraVec              = np.array([1500, 3000, 6000, 8000, 10000, 12000, 14000])
Temp.NTra                = Temp.TraVec.size
Temp.IntVec              = Temp.TraVec
Temp.NInt                = Temp.IntVec.size


## Path to CoarseAIR Output Folder inside the Run Folder
PathToOutput             = '/Users/sventuri/WORKSPACE/CG-QCT/run_O3_All/Test/'
## System Name
SystName                 = 'O3'
UploadSystem             = getattr(ChemicalSystems, SystName + '_Upload')
Syst                     = UploadSystem( Temp ) ### Uploading System Properties

## Kinetic Method for Each of the Molecules (STS/CGM)
Syst.Molecule[0].KinMthd = 'STS'
## Nb of Levels/Bins per Molecule
Syst.Molecule[0].NBins   = 6115
Syst.PathToFolder        = PathToOutput + Syst.Name
##--------------------------------------------------------------------------------------------------------------


Syst = Read_qnsEnBin(Syst)


for iBin in range(Syst.Molecule[0].NBins):

	for iTra in range(Temp.NTra):
		TTra = Temp.TraVec[iTra]
		TInt = TTra

		[Processes, Rates, RatesSD] = Read_RatesCGQCT(Syst, TTra, TInt, iBin)

	for jProc i
