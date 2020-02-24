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
import numpy as np
import pandas
import h5py
from os import path

sys.path.insert(0, '../Saving/')
sys.path.insert(0, '../Loading/')
sys.path.insert(0, '../Computing/')
sys.path.insert(0, '../Parameters/')
sys.path.insert(0, '../Writing/')
sys.path.insert(0, '../Plotting/')

from Saving        import Save_Levels_HDF5, Save_qnsEnBin_HDF5, Save_PartFuncsAndEnergiesAtT_HDF5, Save_RatesAtT_HDF5
from Loading       import Load_Levels_HDF5, Load_qnsEnBin_HDF5, Load_PartFuncsAndEnergiesAtT_HDF5, Load_RatesAtT_HDF5, Load_DissRatesAtT_HDF5
from Computing     import Compute_Correction_To_DissRates, Compute_Rates_Overall, Compute_Rates_Thermal, Compute_BackwardRates, Compute_PrefJumps, Compute_WindAvrg_Matrix, Compute_Correction_To_GroupedDissRates
from Parameters    import *
from Writing       import Write_Rates_Thermal, Write_Kinetics, Write_GroupedKinetics, Write_PrefJumps
from Plotting      import Plot_Rates_Thermal

def Read_Levels(Syst, InputData):
    print('\n  [Read_Levels]: Reading Level Properties')

    for iMol in range(Syst.NMolecules):

        PathToFile  = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        f           = h5py.File(PathToFile, 'a')
        TempStr     = Syst.Molecule[iMol].Name + '/LevelrMin'
        PresentFlg  = TempStr in f.keys()
        f.close()
        if (PresentFlg): print('  [Read_Levels]: Found Data for Molecule Nb ' + str(iMol) + ' (' + Syst.Molecule[iMol].Name + ') in the HDF5 File')

        if ( (PresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):
            print('  [Read_Levels]: Loading Data for Molecule Nb ' + str(iMol) + ' (' + Syst.Molecule[iMol].Name + ') from HDF5 File')
            Load_Levels_HDF5(Syst, iMol)

        else:
            PathToFile = Syst.PathToFolder + '/' + Syst.Molecule[iMol].Name + '/levels_cut.inp'
            print('  [Read_Levels]: Reading File: ', PathToFile)
            
            Data = pandas.read_csv(PathToFile, header=None, skiprows=15, delimiter=r"\s+")
            Data = Data.apply(pandas.to_numeric, errors='coerce')

            Syst.Molecule[iMol].Levelvqn   = np.array(Data.values[:,0],  dtype=np.int16  )
            Syst.Molecule[iMol].Leveljqn   = np.array(Data.values[:,1],  dtype=np.int16  )
            Syst.Molecule[iMol].LevelEEh   = np.array(Data.values[:,2],  dtype=np.float64)
            Syst.Molecule[iMol].LevelEeV   = Syst.Molecule[iMol].LevelEEh * Hartree_To_eV 
            Syst.Molecule[iMol].LevelEgam  = np.array(Data.values[:,3],  dtype=np.float64)
            Syst.Molecule[iMol].LevelrMin  = np.array(Data.values[:,4],  dtype=np.float64)
            Syst.Molecule[iMol].LevelrMax  = np.array(Data.values[:,5],  dtype=np.float64)
            Syst.Molecule[iMol].LevelVMin  = np.array(Data.values[:,6],  dtype=np.float64)
            Syst.Molecule[iMol].LevelVMax  = np.array(Data.values[:,7],  dtype=np.float64)
            Syst.Molecule[iMol].LevelTau   = np.array(Data.values[:,8],  dtype=np.float64)
            Syst.Molecule[iMol].LevelrIn   = np.array(Data.values[:,9],  dtype=np.float64)
            Syst.Molecule[iMol].LevelrOut  = np.array(Data.values[:,10], dtype=np.float64)

            if (InputData.HDF5.Save_Flg):
                print('  [Read_Levels]: Saving Data for Molecule Nb ' + str(iMol) + ' (' + Syst.Molecule[iMol].Name + ') in the HDF5 File')
                Save_Levels_HDF5(Syst, iMol)

        Syst.Molecule[iMol].Compute_ERot()

    return Syst



def Read_qnsEnBin(Syst, InputData):
    print('\n  [Read_qnsEnBin]: Reading Level Q.N.s, Energies, Degeneracies and ToBin Mapping')

    for iMol in range(Syst.NMolecules):

        PathToFile  = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        f           = h5py.File(PathToFile, 'a')
        TempStr     = Syst.Molecule[iMol].Name + '/Levelg'
        PresentFlg  = TempStr in f.keys()
        f.close()
        if (PresentFlg): print('  [Read_qnsEnBin]: Found Data for Molecule Nb ' + str(iMol) + ' (' + Syst.Molecule[iMol].Name + ') in the HDF5 File')

        if ( (PresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):
            print('  [Read_qnsEnBin]: Loading Data for Molecule Nb ' + str(iMol) + ' (' + Syst.Molecule[iMol].Name + ') from HDF5 File')
            Load_qnsEnBin_HDF5(Syst, iMol)

        else:
            PathToFile = Syst.PathToFolder + '/' + Syst.Molecule[iMol].Name + '/' + Syst.Molecule[iMol].Name + '_' + str(Syst.Molecule[iMol].NBins) + '/qnsEnBin.dat'
            print('  [Read_qnsEnBin]: Reading File: ', PathToFile)
            
            Data = pandas.read_csv(PathToFile, header=None, skiprows=1, delimiter=r"\s+")
            Data = Data.apply(pandas.to_numeric, errors='coerce')

            Syst.Molecule[iMol].Levelvqn   = np.array(Data.values[:,1], dtype=np.int16  )
            Syst.Molecule[iMol].Leveljqn   = np.array(Data.values[:,2], dtype=np.int16  )
            Syst.Molecule[iMol].LevelEEh   = np.array(Data.values[:,3], dtype=np.float64)
            Syst.Molecule[iMol].Levelg     = np.array(Data.values[:,4], dtype=np.float64)
            Syst.Molecule[iMol].LevelToBin = np.array(Data.values[:,5], dtype=np.int32  )

            if (InputData.HDF5.Save_Flg):
                print('  [Read_qnsEnBin]: Saving Data for Molecule Nb ' + str(iMol) + ' (' + Syst.Molecule[iMol].Name + ') in the HDF5 File')
                Save_qnsEnBin_HDF5(Syst, iMol)

    return Syst



def Read_PartFuncsAndEnergies(Syst, Temp, InputData):
    print('\n  [Read_PartFuncsAndEnergies]: Reading Level Partition Functions and Energies')

    for iMol in range(Syst.NMolecules):

        for iT in Temp.iTVec:
            TTra = Temp.TranVec[iT-1]
            TInt = TTra

            PathToFile  = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
            f           = h5py.File(PathToFile, 'a')
            TStr        = 'T_' + str(int(TTra)) + '_' + str(int(TInt))
            TempStr     = TStr + '/' + Syst.Molecule[iMol].Name + '/LevelQ'
            TPresentFlg = TempStr in f.keys()
            f.close()
            if (TPresentFlg): print('  [Read_PartFuncsAndEnergies]: Found Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File')

            if ( (TPresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):
                print('  [Read_PartFuncsAndEnergies]: Loading Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) from HDF5 File')
                Load_PartFuncsAndEnergiesAtT_HDF5(Syst, iMol, iT, TTra, TInt)

            else:
                PathToFile = Syst.PathToFolder + '/' + Syst.Molecule[iMol].Name + '/' + Syst.Molecule[iMol].Name + '_' + str(Syst.Molecule[iMol].NBins) + '/T' + str(int(Temp.TranVec[iT-1])) + '.dat'
                print('  [Read_PartFuncsAndEnergies]: Reading File: ', PathToFile)

                Data = pandas.read_csv(PathToFile, header=None, skiprows=1, delimiter=r"\s+")
                Data = Data.apply(pandas.to_numeric, errors='coerce')

                Syst.Molecule[iMol].T[iT-1].QRatio   = np.array(Data.values[:,0])
                Syst.Molecule[iMol].T[iT-1].Q        = np.array(Data.values[:,1])
                Syst.Molecule[iMol].T[iT-1].LevelEeV = np.array(Data.values[:,2])


                if (InputData.HDF5.Save_Flg):
                    print('  [Read_Rates_CGQCT]: Saving Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File')
                    Save_PartFuncsAndEnergiesAtT_HDF5(Syst, iMol, iT, TTra, TInt)

            Syst.Molecule[iMol].T[iT-1].LevelQExp = Syst.Molecule[iMol].Levelg * np.exp( - Syst.Molecule[iMol].T[iT-1].LevelEeV * Ue / (TTra * UKb) )
            Syst.Molecule[iMol].T[iT-1].LevelQExp = Syst.Molecule[iMol].T[iT-1].LevelQExp / np.sum(Syst.Molecule[iMol].T[iT-1].LevelQExp)

    return Syst



def Read_RatesFile_CGQCT(Syst, TTra, TInt, iLevel):
    #sed -i 's/D-/E-/g' *

    PathToFile = Syst.PathToFolder + '/' + Syst.Molecule[0].Name + '/Rates/T_' + str(int(TTra)) + '_' + str(int(TInt)) + '/Bin' + str(iLevel+1) + '.dat'
    #print(PathToFile)
    if (path.isfile(PathToFile)):
        count   = 0
        TempFlg = False
        with open(PathToFile, 'r') as f:
            for line in f:
                count +=1
                if (count>5):
                    TempFlg = True
                    break
        f.close()

        if (TempFlg):
            Data  = pandas.read_csv(PathToFile, header=None, skiprows=5, delimiter=r"\s+")
            Data  = Data.apply(pandas.to_numeric, errors='coerce')

            ProcessesTemp = np.array(Data[1].values, dtype=np.int64)
            RatesTemp     = np.array(Data[2].values, dtype=np.float64)
            RatesSDTemp   = np.array(Data[3].values, dtype=np.float64)
        else:
            print('      [Read_RatesFile_CGQCT]: Rate File Does Not Contain Rates: ' + PathToFile )
            ProcessesTemp = np.zeros(1, dtype=int)+1
            RatesTemp     = np.zeros(1)
            RatesSDTemp   = np.zeros(1)
    else:
        print('      [Read_RatesFile_CGQCT]: Rate File Does Not Exist: ' + PathToFile )
        ProcessesTemp = np.zeros(1, dtype=int)+1
        RatesTemp     = np.zeros(1)
        RatesSDTemp   = np.zeros(1)

    return ProcessesTemp, RatesTemp, RatesSDTemp



def Read_RatesAtT_CGQCT(Syst, TTra, TInt, iT):

    Syst.T[iT-1].Proc[0].Rates = np.zeros((Syst.Molecule[0].NBins, 3))
    Syst.T[iT-1].Proc[1].Rates = np.zeros((Syst.Molecule[0].NBins, Syst.Molecule[0].NBins))
    for iProc in range(2, 4):
        Syst.T[iT-1].Proc[iProc].Rates       = np.zeros((Syst.Molecule[0].NBins, Syst.Molecule[Syst.Pair[iProc-1].ToMol].NBins))
    for iProc in range(2, Syst.NProcTypes):
        Syst.T[iT-1].ProcExch[iProc-2].Rates = np.zeros((Syst.Molecule[0].NBins, Syst.Molecule[Syst.ExchtoMol[iProc-2]].NBins))

    i  = 0
    ii = 10
    for iBins in range(Syst.Molecule[0].NBins):
        if (int(iBins/Syst.Molecule[0].NBins*100) == int(i*ii)):
            print('    [Read_RatesAtT_CGQCT]: Read ' + str(i*ii) + '% of the Rate Files')
            i = i+1 
        RatesTempAll                            = np.zeros(Syst.Pair[-1].NProcTot+1)
        [ProcessesTemp, RatesTemp, RatesSDTemp] = Read_RatesFile_CGQCT(Syst, TTra, TInt, iBins)
        RatesTempAll[ProcessesTemp[:]-1] = RatesTemp[:]

        RatesSplitted                       = np.split( RatesTempAll, np.array([1, Syst.Pair[0].NProcTot, Syst.Pair[1].NProcTot, Syst.Pair[2].NProcTot]) )
        Syst.T[iT-1].Proc[0].Rates[iBins,0] = RatesSplitted[0]
        Syst.T[iT-1].Proc[1].Rates[iBins,:] = RatesSplitted[1]
        for iProc in range(2, 4):
            Syst.T[iT-1].Proc[iProc].Rates[iBins,:] = RatesSplitted[iProc]
    
    for iProc in range(2, 4):
        for jProc in range(2, Syst.NProcTypes):
            if (Syst.Pair[iProc-1].ToMol == Syst.ExchtoMol[jProc-2] ):
                Syst.T[iT-1].ProcExch[jProc-2].Rates = Syst.T[iT-1].ProcExch[jProc-2].Rates + Syst.T[iT-1].Proc[iProc].Rates

    return Syst



def Read_Rates_CGQCT(Syst, Temp, InputData):

    print('\n  [Read_Rates_CGQCT]: Uploading Rates')

    if (InputData.Kin.WindAvrgFlg) and (InputData.Kin.Write_Flg):
        print('\n  [Read_Rates_CGQCT]: Computing Data for Window-Averaging the Kinetics')
        Syst = Compute_WindAvrg_Matrix(Syst, InputData)
    
    for iT in Temp.iTVec:
        TTra = Temp.TranVec[iT-1]
        TInt = TTra

        PathToFile  = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        f           = h5py.File(PathToFile, 'a')
        TStr        = 'T_' + str(int(TTra)) + '_' + str(int(TInt)) + '/Rates/'
        TPresentFlg = TStr in f.keys()
        f.close()
        if (TPresentFlg): print('  [Read_Rates_CGQCT]: Found Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File')

        if ( (TPresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):
            print('  [Read_Rates_CGQCT]: Loading Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) from HDF5 File')
        
            Load_RatesAtT_HDF5(Syst, TTra, TInt, iT)

        else:
            print('  [Read_Rates_CGQCT]: Reading Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) from set of .dat Files')

            Syst = Read_RatesAtT_CGQCT(Syst, TTra, TInt, iT)

            if (InputData.HDF5.Save_Flg):
                print('  [Read_Rates_CGQCT]: Saving Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File\n')
                Save_RatesAtT_HDF5(Syst, iT, TTra, TInt)
        

        if ((InputData.Rates.PrefJumps_Flg) or (InputData.Kin.Groups.Flg)):
            print('  [Read_Rates_CGQCT]: Computing Backweard Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)')
            Syst = Compute_BackwardRates(Syst, iT)


        if (InputData.Kin.Groups.Flg):
            print('  [Read_Rates_CGQCT]: Computing Grouped Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)')
            NbLevels = np.array( [Syst.Molecule[0].NBins]            + [Syst.Molecule[Syst.ExchtoMol[iExch-2]].NBins            for iExch in range(2, Syst.NProcTypes)] ) 
            NbGroups = np.array( [Syst.Molecule[0].Grouped.NbGroups] + [Syst.Molecule[Syst.ExchtoMol[iExch-2]].Grouped.NbGroups for iExch in range(2, Syst.NProcTypes)] )
            Syst.Molecule[0].Grouped.Compute_GroupRates(Syst, iT-1, NbLevels, NbGroups)


        if (InputData.Kin.CorrFactor != 1.0):	
                print('  [Read_Rates_CGQCT]: Correcting Dissociation Rates (Corr Factor = ' + str(InputData.Kin.CorrFactor) + ') for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
                Syst = Compute_Correction_To_DissRates(InputData, Syst, iT)
                if (InputData.Kin.Groups.Flg):
                    print('  [Read_Rates_CGQCT]: Correcting Grouped Dissociation Rates (Corr Factor = ' + str(InputData.Kin.CorrFactor) + ') for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
                    Syst = Compute_Correction_To_GroupedDissRates(InputData, Syst, iT)


        print('  [Read_Rates_CGQCT]: Computing Overall Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
        Syst = Compute_Rates_Overall(Syst, iT)


        print('  [Read_Rates_CGQCT]: Computing Thermal Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
        Syst = Compute_Rates_Thermal(Syst, iT)


        if (InputData.Rates.PrefJumps_Flg):
            print('  [Read_Rates_CGQCT]: Computing Preferred Jumps for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)')
            Syst = Compute_PrefJumps(InputData, Syst, iT)
            print('  [Read_Rates_CGQCT]: Writing   Preferred Jumps for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)')
            Write_PrefJumps(Syst, Temp, InputData, iT)
            print('\n')


        if (InputData.Kin.Write_Flg):
            print('  [Read_Rates_CGQCT]: Writing Kinetics File for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
            Write_Kinetics(Syst, Temp, InputData, iT)


        if  (InputData.Kin.Groups.Flg) and (InputData.Kin.Groups.Write_Flg):
            print('  [Read_Rates_CGQCT]: Writing Grouped Kinetics File for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
            NbGroups = np.array( [Syst.Molecule[0].Grouped.NbGroups] + [Syst.Molecule[Syst.ExchtoMol[iExch-2]].Grouped.NbGroups for iExch in range(2, Syst.NProcTypes)] )
            Write_GroupedKinetics(Syst, Temp, InputData, iT, NbGroups)


        if (InputData.DelRateMat_Flg):
            for iProc in range(4):
                del Syst.T[iT-1].Proc[iProc].Rates
            for iProc in range(2, Syst.NProcTypes):
                del Syst.T[iT-1].ProcExch[iProc-2].Rates

    print('  [Read_Rates_CGQCT]: Saving Thermal Rates\n')
    Write_Rates_Thermal(Syst, Temp, InputData)

    print('  [Read_Rates_CGQCT]: Plotting Thermal Rates\n')
    Plot_Rates_Thermal(Syst, Temp, InputData)

    return Syst