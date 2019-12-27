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
import numpy as np
import pandas
import csv
import sys
import os

from System            import system
from Temperatures      import temperatures

import ChemicalSystems

def Initialize_Data(InputData):

    Temp                     = temperatures()
    ## Vector of Temperatures
    Temp.TranVec             = InputData.TranVec
    Temp.NTran               = Temp.TranVec.size
    Temp.IntVec              = Temp.TranVec
    Temp.NInt                = Temp.IntVec.size

    Temp.iTVec               = InputData.iTVec


    ## Path to CoarseAIR Output Folder inside the Run Folder
    OutputFldr               = InputData.FinalFldr
    ## System Name
    SystNameLong             = InputData.SystNameLong
    UploadSystem             = getattr(ChemicalSystems, SystNameLong + '_Upload')
    Syst                     = UploadSystem( Temp ) ### Uploading System Properties

    
    for iMol in range(Syst.NMolecules):
        if hasattr(InputData, 'KinMthd'):
            ## Kinetic Method for Each of the Molecules (STS/CGM)
            Syst.Molecule[iMol].KinMthd = InputData.KinMthd[iMol]
        if hasattr(InputData, 'NBins'):
            ## Nb of Levels/Bins per Molecule
            Syst.Molecule[iMol].NBins   = InputData.NBins[iMol]
            

    NProcTot = 1
    for iP in range(3):
        Syst.Pair[iP].NBins    = Syst.Molecule[Syst.Pair[iP].ToMol].NBins
        NProcTot               = NProcTot + Syst.Pair[iP].NBins
        Syst.Pair[iP].NProcTot = NProcTot
    Syst.PathToFolder = OutputFldr + Syst.Name
    Syst.PathToHDF5   = InputData.HDF5.ReadFldr


    ## Creating Output Folder
    InputData.FinalFldr = InputData.FinalFldr + '/' + SystNameLong + '/'
    if not os.path.exists(InputData.FinalFldr):
        os.makedirs(InputData.FinalFldr)

    return Syst, Temp