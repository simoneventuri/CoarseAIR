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



def Save_Levels_HDF5(Syst, iMol):

    PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    f          = h5py.File(PathToFile, 'a')

    TempStr1   = Syst.Molecule[iMol].Name
    TempStr2   = TempStr1 + '/LevelrMin'
    if TempStr2 in f.keys():
        grp       = f[TempStr1]
        Data      = grp["Levelvqn"]
        Data[...] = Syst.Molecule[iMol].Levelvqn
        Data      = grp["Leveljqn"]
        Data[...] = Syst.Molecule[iMol].Leveljqn
        Data      = grp["LevelEEh"]
        Data[...] = Syst.Molecule[iMol].LevelEEh
        Data      = grp["LevelEeV"]
        Data[...] = Syst.Molecule[iMol].LevelEeV
        Data      = grp["LevelEgam"]
        Data[...] = Syst.Molecule[iMol].LevelEgam
        Data      = grp["LevelrMin"]
        Data[...] = Syst.Molecule[iMol].LevelrMin
        Data      = grp["LevelrMax"]
        Data[...] = Syst.Molecule[iMol].LevelrMax
        Data      = grp["LevelVMin"]
        Data[...] = Syst.Molecule[iMol].LevelVMin
        Data      = grp["LevelVMax"]
        Data[...] = Syst.Molecule[iMol].LevelVMax
        Data      = grp["LevelTau"]
        Data[...] = Syst.Molecule[iMol].LevelTau
        Data      = grp["LevelrIn"]
        Data[...] = Syst.Molecule[iMol].LevelrIn
        Data      = grp["LevelrOut"]
        Data[...] = Syst.Molecule[iMol].LevelrOut

    else:
        if (not (TempStr1 in f.keys())):
            grp       = f.create_group(TempStr1)
        else:
            grp       = f[TempStr1]
        Levelvqn  = grp.create_dataset("Levelvqn",  data=Syst.Molecule[iMol].Levelvqn,  compression="gzip", compression_opts=9)
        Leveljqn  = grp.create_dataset("Leveljqn",  data=Syst.Molecule[iMol].Leveljqn,  compression="gzip", compression_opts=9)
        LevelEEh  = grp.create_dataset("LevelEEh",  data=Syst.Molecule[iMol].LevelEEh,  compression="gzip", compression_opts=9)
        LevelEeV  = grp.create_dataset("LevelEeV",  data=Syst.Molecule[iMol].LevelEeV,  compression="gzip", compression_opts=9)
        LevelEgam = grp.create_dataset("LevelEgam", data=Syst.Molecule[iMol].LevelEgam, compression="gzip", compression_opts=9)
        LevelrMin = grp.create_dataset("LevelrMin", data=Syst.Molecule[iMol].LevelrMin, compression="gzip", compression_opts=9)
        LevelrMax = grp.create_dataset("LevelrMax", data=Syst.Molecule[iMol].LevelrMax, compression="gzip", compression_opts=9)
        LevelVMin = grp.create_dataset("LevelVMin", data=Syst.Molecule[iMol].LevelVMin, compression="gzip", compression_opts=9)
        LevelVMax = grp.create_dataset("LevelVMax", data=Syst.Molecule[iMol].LevelVMax, compression="gzip", compression_opts=9)
        LevelTau  = grp.create_dataset("LevelTau",  data=Syst.Molecule[iMol].LevelTau,  compression="gzip", compression_opts=9)
        LevelrIn  = grp.create_dataset("LevelrIn",  data=Syst.Molecule[iMol].LevelrIn,  compression="gzip", compression_opts=9)
        LevelrOut = grp.create_dataset("LevelrOut", data=Syst.Molecule[iMol].LevelrOut, compression="gzip", compression_opts=9)

    f.close()


def Save_qnsEnBin_HDF5(Syst, iMol):

    PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    f          = h5py.File(PathToFile, 'a')

    TempStr1   = Syst.Molecule[iMol].Name
    TempStr2   = TempStr1 + '/Levelg'
    if TempStr2 in f.keys():

        grp       = f[TempStr1]
        # Data      = grp["Levelvqn"]
        # Data[...] = Syst.Molecule[iMol].Levelvqn
        # Data      = grp["Leveljqn"]
        # Data[...] = Syst.Molecule[iMol].Leveljqn
        # Data      = grp["LevelEEh"]
        # Data[...] = Syst.Molecule[iMol].LevelEEh
        Data      = grp["Levelg"]
        Data[...] = Syst.Molecule[iMol].Levelg
        Data      = grp["LevelToBin"]
        Data[...] = Syst.Molecule[iMol].LevelToBin

    else:
        if (not (TempStr1 in f.keys())):
            grp       = f.create_group(TempStr1)
        else:
            grp       = f[TempStr1]
        #Levelvqn   = grp.create_dataset("Levelvqn",   data=Syst.Molecule[iMol].Levelvqn,   compression="gzip", compression_opts=9)
        #Leveljqn   = grp.create_dataset("Leveljqn",   data=Syst.Molecule[iMol].Leveljqn,   compression="gzip", compression_opts=9)
        #LevelEEh   = grp.create_dataset("LevelEEh",   data=Syst.Molecule[iMol].LevelEEh,   compression="gzip", compression_opts=9)
        Levelg     = grp.create_dataset("Levelg",     data=Syst.Molecule[iMol].Levelg,     compression="gzip", compression_opts=9)
        LevelToBin = grp.create_dataset("LevelToBin", data=Syst.Molecule[iMol].LevelToBin, compression="gzip", compression_opts=9)

    f.close()



def Save_PartFuncsAndEnergiesAtT_HDF5(Syst, iMol, iT, TTra, TInt):

    PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    f          = h5py.File(PathToFile, 'a')

    TStr     = 'T_' + str(int(TTra)) + '_' + str(int(TInt))
    TempStr1 = TStr + '/' + Syst.Molecule[iMol].Name 
    TempStr2 = TempStr1 + '/LevelQ'     
    if TempStr2 in f.keys():
        grp       = f[TempStr1]
        Data      = grp["LevelEeV"]
        Data[...] = Syst.Molecule[iMol].T[iT-1].LevelEeV
        Data      = grp["LevelQ"]
        Data[...] = Syst.Molecule[iMol].T[iT-1].Q
    
    else:
        if (not (TempStr1 in f.keys())):
            grp       = f.create_group(TempStr1)
        else:
            grp       = f[TempStr1]
        LevelEeV  = grp.create_dataset("LevelEeV", data=Syst.Molecule[iMol].T[iT-1].LevelEeV, compression="gzip", compression_opts=9)
        LevelQ    = grp.create_dataset("LevelQ",   data=Syst.Molecule[iMol].T[iT-1].Q,        compression="gzip", compression_opts=9)

    f.close()



def Save_RatesAtT_HDF5(Syst, iT, TTra, TInt):

    PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    f          = h5py.File(PathToFile, 'a')

    TStr = 'T_' + str(int(TTra)) + '_' + str(int(TInt)) + '/Rates/'       
    if TStr in f.keys():
        grp               = f[TStr]

        Data      = grp["Diss"]
        Data[...] = Syst.T[iT-1].Proc[0].Rates
        Data      = grp["Inel"]
        Data[...] = Syst.T[iT-1].Proc[1].Rates
        # for iProc in range(2, 4):
        #     ExchStr       = "Rates_Exch_" + str(iProc-1)
        #     Data          = grp[ExchStr]
        #     Data[...]     = Syst.T[iT-1].Proc[iProc].Rates

        for iProc in range(2, Syst.NProcTypes):
            ExchStr   = "Exch_" + str(iProc-1)
            Data      = grp[ExchStr]
            Data[...] = Syst.T[iT-1].ProcExch[iProc-2].Rates

        # Data       = grp["Overall_Diss"]
        # Data[...]  = Syst.T[iT-1].ProcTot[0].Rates
        # # Data       = grp["Rates_Overall_Inel"]
        # # Data[...]  = Syst.T[iT-1].ProcTot[1].Rates
        # for iProc in range(2, Syst.NProcTypes):
        #     ExchStr      = "Overall_Exch_" + str(iProc-1)
        #     Data         = grp[ExchStr] 
        #     Data[...]    = Syst.T[iT-1].ProcTot[iProc].Rates
    
    else:
        grp           = f.create_group(TStr)

        Proc0         = grp.create_dataset("Diss", data=Syst.T[iT-1].Proc[0].Rates, compression="gzip", compression_opts=9)
        Proc1         = grp.create_dataset("Inel", data=Syst.T[iT-1].Proc[1].Rates, compression="gzip", compression_opts=9)
        # for iProc in range(2, 4):
        #     ExchStr   = "Rates_Exch_" + str(iProc-1)
        #     Proci     = grp.create_dataset(ExchStr, data=Syst.T[iT-1].Proc[iProc].Rates, compression="gzip", compression_opts=9)

        for iProc in range(2, Syst.NProcTypes):
            ExchStr   = "Exch_" + str(iProc-1)
            ProcExchi = grp.create_dataset(ExchStr, data=Syst.T[iT-1].ProcExch[iProc-2].Rates, compression="gzip", compression_opts=9)

        # ProcTot0      = grp.create_dataset("Overall_Diss", data=Syst.T[iT-1].ProcTot[0].Rates, compression="gzip", compression_opts=9)
        # ProcTot1      = grp.create_dataset("Overall_Inel", data=Syst.T[iT-1].ProcTot[1].Rates, compression="gzip", compression_opts=9) 
        # for iProc in range(2, Syst.NProcTypes):
        #     ExchStr   = "Overall_Exch_" + str(iProc-1)
        #     ProcToti  = grp.create_dataset(ExchStr, data=Syst.T[iT-1].ProcTot[iProc].Rates, compression="gzip", compression_opts=9)

    f.close()