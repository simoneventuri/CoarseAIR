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
import sys
import os

import numpy as np


from System            import system
from Temperatures      import temperatures
from Molecule          import groupedmolecule

import ChemicalSystems

def Initialize_Data(InputData):

    ## Initializing Vectors of Temperatures
    print("  [Initializing.py - Initialize_Data]: Initializing Vectors of Temperatures")    
    Temp                     = temperatures()
    Temp.TranVec             = InputData.TranVec
    Temp.NTran               = Temp.TranVec.size
    Temp.IntVec              = Temp.TranVec
    Temp.NInt                = Temp.IntVec.size
    Temp.T0                  = InputData.T0
    Temp.iTVec               = InputData.iTVec


    ## Initializing System
    print("  [Initializing.py - Initialize_Data]: Initializing System")    
    SystNameLong             = InputData.SystNameLong
    UploadSystem             = getattr(ChemicalSystems, SystNameLong + '_Upload')
    Syst                     = UploadSystem( Temp ) ### Uploading System Properties


    ## Initializng Paths
    print("  [Initializing.py - Initialize_Data]: Initializng Paths")    
    Syst.PathToReadFldr = InputData.DtbReadFldr + Syst.Name                                                 ## i.e.: CoarseAIR Output Folder / Run Folder / System Folder        
    Syst.PathToHDF5     = InputData.HDF5.DtbFldr
    if not os.path.exists(Syst.PathToHDF5):
        os.makedirs(Syst.PathToHDF5)
    Syst.HDF5Exist_Flg = True
    PathToFile         = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    if not os.path.isfile(PathToFile):
        Syst.HDF5Exist_Flg = False
        print("  [Initializing.py - Initialize_Data]: WARNING: The HDF5 File " + PathToFile + " corresponding to the System " + Syst.NameLong + " does not exist. I will create a new one." )    


    ## Creating Output Folders
    print("  [Initializing.py - Initialize_Data]: Creating Output Folders")    
    InputData.OutputWriteFldr = InputData.OutputWriteFldr + '/' + SystNameLong + '/'
    if not os.path.exists(InputData.OutputWriteFldr):
        os.makedirs(InputData.OutputWriteFldr)
    InputData.OutputWriteFldr = InputData.OutputWriteFldr + '/' + InputData.ME.ProcCode + '/'
    if not os.path.exists(InputData.OutputWriteFldr):
        os.makedirs(InputData.OutputWriteFldr)


    ##
    Syst.EqNStatesIn  = np.zeros((Syst.NMolecules), dtype=np.int64)
    Syst.EqNStatesOut = np.zeros((Syst.NMolecules), dtype=np.int64)
    print("  [Initializing.py - Initialize_Data]: Initializing Molecules and Counting Nb of Levels / Groups per Molecule")    
    for iMol in range(Syst.NMolecules):

        print("  [Initializing.py - Initialize_Data]:   Initializing Molecule Nb " + str(iMol))    
        Syst.Molecule[iMol].Initialize( InputData, Syst, Temp, iMol)
        print("  [Initializing.py - Initialize_Data]:     Done Initializing Molecule Nb " + str(iMol))    

        Syst.EqNStatesIn[iMol]  = Syst.Molecule[iMol].EqNStatesIn
        Syst.EqNStatesOut[iMol] = Syst.Molecule[iMol].EqNStatesOut


    ## 
    print("  [Initializing.py - Initialize_Data]: Initializing Quantities for Arrhenius Fitting")    
    Syst.Arr.MinRate  = InputData.Kin.MinRate
    Syst.Arr.MinNbTs  = InputData.Kin.MinNbTs


    return Syst, Temp