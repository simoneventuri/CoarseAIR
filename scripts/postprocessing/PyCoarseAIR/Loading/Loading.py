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
import h5py


def Load_Levels_HDF5(Syst, iMol):

    PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    f          = h5py.File(PathToFile, 'a')

    TempStr = Syst.Molecule[iMol].Name   
    grp     = f[TempStr]

    Data                          = grp["Levelvqn"]
    Syst.Molecule[iMol].Levelvqn  = Data[...]
    Data                          = grp["Leveljqn"]
    Syst.Molecule[iMol].Leveljqn  = Data[...]
    Data                          = grp["LevelEEh"]
    Syst.Molecule[iMol].LevelEEh  = Data[...]
    Data                          = grp["LevelEeV"]
    Syst.Molecule[iMol].LevelEeV  = Data[...]
    Data                          = grp["LevelEgam"]
    Syst.Molecule[iMol].LevelEgam = Data[...]
    Data                          = grp["LevelrMin"]
    Syst.Molecule[iMol].LevelrMin = Data[...]
    Data                          = grp["LevelrMax"]
    Syst.Molecule[iMol].LevelrMax = Data[...]
    Data                          = grp["LevelVMin"]
    Syst.Molecule[iMol].LevelVMin = Data[...]
    Data                          = grp["LevelVMax"]
    Syst.Molecule[iMol].LevelVMax = Data[...]
    Data                          = grp["LevelTau"]
    Syst.Molecule[iMol].LevelTau  = Data[...]
    Data                          = grp["LevelrIn"]
    Syst.Molecule[iMol].LevelrIn  = Data[...]
    Data                          = grp["LevelrOut"]
    Syst.Molecule[iMol].LevelrOut = Data[...]

    f.close()



def Load_qnsEnBin_HDF5(Syst, iMol):

    PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    f          = h5py.File(PathToFile, 'a')

    TempStr = Syst.Molecule[iMol].Name   
    grp     = f[TempStr]

    # Data                           = grp["Levelvqn"]
    # Syst.Molecule[iMol].Levelvqn   = Data[...]
    # Data                           = grp["Leveljqn"]
    # Syst.Molecule[iMol].Leveljqn   = Data[...]
    # Data                           = grp["LevelEEh"]
    # Syst.Molecule[iMol].LevelEEh   = Data[...]
    Data                           = grp["Levelg"]
    Syst.Molecule[iMol].Levelg     = Data[...]
    Data                           = grp["LevelToBin"]
    Syst.Molecule[iMol].LevelToBin = Data[...]

    f.close()



def Load_PartFuncsAndEnergiesAtT_HDF5(Syst, iMol, iT, TTra, TInt):

    PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    f          = h5py.File(PathToFile, 'a')

    TStr    = 'T_' + str(int(TTra)) + '_' + str(int(TInt))
    TempStr = TStr + '/' + Syst.Molecule[iMol].Name   
    grp     = f[TempStr]

    Data                                 = grp["LevelEeV"]
    Syst.Molecule[iMol].T[iT-1].LevelEeV = Data[...]
    Data                                 = grp["LevelQ"]
    Syst.Molecule[iMol].T[iT-1].Q        = Data[...] 

    f.close()



def Load_RatesAtT_HDF5(Syst, TTra, TInt, iT):

    PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    f          = h5py.File(PathToFile, "r")

    TStr = 'T_' + str(int(TTra)) + '_' + str(int(TInt)) + '/Rates/'
    grp  = f[TStr]

    Data                       = grp["Diss"]
    Syst.T[iT-1].Proc[0].Rates = Data[...]
    Data                       = grp["Inel"]
    Syst.T[iT-1].Proc[1].Rates = Data[...]
    # for iProc in range(2, 4):
    #     ExchStr                        = "Rates_Exch_" + str(iProc-1)
    #     Data                           = grp[ExchStr]
    #     Syst.T[iT-1].Proc[iProc].Rates = Data[...]

    for iProc in range(2, Syst.NProcTypes):
        ExchStr                              = "Exch_" + str(iProc-1)
        Data                                 = grp[ExchStr]
        Syst.T[iT-1].ProcExch[iProc-2].Rates = Data[...]

    # Data                          = grp["Overall_Diss"] 
    # Syst.T[iT-1].ProcTot[0].Rates = Data[...]
    # # Data                          = grp["Overall_Inel"]
    # # Syst.T[iT-1].ProcTot[1].Rates = Data[...]
    # for iProc in range(2, Syst.NProcTypes):
    #     ExchStr                           = "Overall_Exch_" + str(iProc-1)
    #     Data                              = grp[ExchStr]
    #     Syst.T[iT-1].ProcTot[iProc].Rates = Data[...]

    f.close()

    return Syst



def Load_DissRatesAtT_HDF5(Syst, TTra, TInt, iT):

    PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    f          = h5py.File(PathToFile, "r")

    TStr = 'T_' + str(int(TTra)) + '_' + str(int(TInt)) + '/Rates/'
    grp  = f[TStr]

    Data                       = grp["Diss"]
    Syst.T[iT-1].Proc[0].Rates = Data[...]

    f.close()

    return Syst