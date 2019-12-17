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
#WORKSPACE = os.environ['WORKSPACE_PATH']
#SrcDir    = '/home/venturi/WORKSPACE//CoarseAIR/coarseair/' #os.environ['COARSEAIR_SOURCE_DIR'] 
SrcDir    = '/Users/sventuri/WORKSPACE//CoarseAIR/coarseair/' #os.environ['COARSEAIR_SOURCE_DIR'] 
sys.path.insert(0, SrcDir + '/scripts/postprocessing/PyCoarseAIR/Objects/')
sys.path.insert(0, SrcDir + '/scripts/postprocessing/PyCoarseAIR/ChemicalSystems/')
sys.path.insert(0, SrcDir + '/scripts/postprocessing/PyCoarseAIR/Reading/')
sys.path.insert(0, SrcDir + '/scripts/postprocessing/PyCoarseAIR/Writing/')
sys.path.insert(0, SrcDir + '/scripts/postprocessing/PyCoarseAIR/Computing/')
print(sys.argv)

if (len(sys.argv) > 1):
    InputFile = sys.argv[1]
    print("Calling PyCoarseAIR with Input File = ", InputFile)
    sys.path.insert(0, InputFile)
else:
    print("Calling PyCoarseAIR without Input File")
    sys.path.insert(0, SrcDir + '/scripts/postprocessing/PyCoarseAIR/InputData/')

import ChemicalSystems
from System            import system
from Temperatures      import temperatures

from Reading           import Read_PartFuncsAndEnergies, Read_qnsEnBin, Read_Rates_CGQCT
from Writing           import Write_PartFuncsAndEnergies, Write_Kinetics_FromOverall
from Computing         import Compute_ThermalDissRates_FromOverall
from CoarseAIR         import InputData
##--------------------------------------------------------------------------------------------------------------



################################################################################################################
Temp                     = temperatures()
## Vector of Temperatures
Temp.TranVec             = InputData.TranVec
Temp.NTran               = Temp.TranVec.size
Temp.IntVec              = Temp.TranVec
Temp.NInt                = Temp.IntVec.size

Temp.iTVec               = InputData.iTVec


## Path to CoarseAIR Output Folder inside the Run Folder
OutputFldr               = InputData.OutputFldr
## System Name
SystName                 = InputData.SystName
UploadSystem             = getattr(ChemicalSystems, SystName + '_Upload')
Syst                     = UploadSystem( Temp ) ### Uploading System Properties

## Kinetic Method for Each of the Molecules (STS/CGM)
for iMol in range(Syst.NMolecules):
    Syst.Molecule[iMol].KinMthd = InputData.KinMthd[iMol]
    ## Nb of Levels/Bins per Molecule
    Syst.Molecule[iMol].NBins   = InputData.NBins[iMol]
NProcTot = 1
for iP in range(3):
	Syst.Pair[iP].NBins    = Syst.Molecule[Syst.Pair[iP].ToMol].NBins
	NProcTot               = NProcTot + Syst.Pair[iP].NBins
	Syst.Pair[iP].NProcTot = NProcTot
Syst.PathToFolder = OutputFldr + Syst.Name
##--------------------------------------------------------------------------------------------------------------




Syst = Read_qnsEnBin(Syst)
Syst = Read_PartFuncsAndEnergies(Syst, Temp)

Syst = Read_Rates_CGQCT(Syst, Temp)

#Write_PartFuncsAndEnergies(Syst, Temp, InputData)
#Write_Kinetics_FromOverall(Syst, Temp, InputData)

#Compute_ThermalDissRates_FromOverall(Syst, Temp, InputData)