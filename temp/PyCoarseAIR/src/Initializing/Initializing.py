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
from Molecule          import groupedmolecule

import ChemicalSystems


def Initialize_Data(InputData):

    Temp                     = temperatures()
    ## Vector of Temperatures
    Temp.TranVec             = InputData.TranVec
    Temp.NTran               = Temp.TranVec.size
    Temp.IntVec              = Temp.TranVec
    Temp.NInt                = Temp.IntVec.size
    Temp.T0                  = InputData.T0
    Temp.iTVec               = InputData.iTVec


    ## Path to CoarseAIR Output Folder inside the Run Folder
    OutputFldr               = InputData.QCTOutFldr
    ## System Name
    SystNameLong             = InputData.SystNameLong
    UploadSystem             = getattr(ChemicalSystems, SystNameLong + '_Upload')
    Syst                     = UploadSystem( Temp ) ### Uploading System Properties


    Syst.PathToFolder = OutputFldr + Syst.Name
    Syst.PathToHDF5   = InputData.HDF5.ReadFldr
    if not os.path.exists(Syst.PathToHDF5):
        os.makedirs(Syst.PathToHDF5)

    Syst.Arr.MinRate  = InputData.Kin.MinRate
    Syst.Arr.MinNbTs  = InputData.Kin.MinNbTs


    ## Creating Output Folders
    InputData.FinalFldr = InputData.FinalFldr + '/' + SystNameLong + '/'
    if not os.path.exists(InputData.FinalFldr):
        os.makedirs(InputData.FinalFldr)
    InputData.FinalFldr = InputData.FinalFldr + '/' + InputData.ME.ProcCode + '/'
    if not os.path.exists(InputData.FinalFldr):
        os.makedirs(InputData.FinalFldr)


    Syst.EqNStatesIn  = np.zeros((Syst.NMolecules), dtype=np.int64)
    Syst.EqNStatesOut = np.zeros((Syst.NMolecules), dtype=np.int64)
    print("  [Initializing.py - Initialize_Data]: Initializing Molecules")    
    for iMol in range(Syst.NMolecules):

        print("  [Initializing.py - Initialize_Data]:   Initializing Molecule Nb " + str(iMol))    
        Syst.Molecule[iMol].Initialize( InputData, Syst, Temp, iMol)
        print("  [Initializing.py - Initialize_Data]:     Done Initializing Molecule Nb " + str(iMol))    

        Syst.EqNStatesIn[iMol]  = Syst.Molecule[iMol].EqNStatesIn
        Syst.EqNStatesOut[iMol] = Syst.Molecule[iMol].EqNStatesOut


    return Syst, Temp
