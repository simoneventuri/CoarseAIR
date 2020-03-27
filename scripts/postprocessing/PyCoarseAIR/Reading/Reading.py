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

from Saving        import Save_Levels_HDF5, Save_qnsEnBin_HDF5, Save_PartFuncsAndEnergiesAtT_HDF5, Save_RatesAtT_HDF5_3Atoms, Save_RatesAtT_HDF5_4Atoms
from Loading       import Load_Levels_HDF5, Load_qnsEnBin_HDF5, Load_PartFuncsAndEnergiesAtT_HDF5, Load_RatesAtT_HDF5_3Atoms, Load_RatesAtT_HDF5_4Atoms, Load_DissRatesAtT_HDF5
from Computing     import Compute_Correction_To_DissRates, Compute_Rates_Overall, Compute_Rates_Thermal, Compute_BackwardRates, Compute_PrefJumps, Compute_WindAvrg_Matrix, Compute_Correction_To_GroupedDissRates
from Parameters    import *
from Writing       import Write_Rates_Thermal, Write_Kinetics_3Atoms, Write_Kinetics_4Atoms, Write_GroupedKinetics, Write_PrefJumps
from Plotting      import Plot_Rates_Thermal

def Read_Levels(Syst, InputData):
    print('\n  [Read_Levels]: Reading Level Properties')

    for iMol in range(Syst.NMolecules):

        PathToFile  = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        f           = h5py.File(PathToFile, 'a')
        TempStr     = Syst.Molecule[iMol].Name + '/LevelrMin'
        print('  [Read_Levels]: Found Data for Molecule Nb ' + str(iMol) + ' (' + Syst.Molecule[iMol].Name + ') in the HDF5 File')

        PresentFlg  = TempStr in f.keys()
        f.close()
        if (PresentFlg): print('  [Read_Levels]: Found Data for Molecule Nb ' + str(iMol) + ' (' + Syst.Molecule[iMol].Name + ') in the HDF5 File')

        if ( (PresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):
            print('  [Read_Levels]: Loading Data for Molecule Nb ' + str(iMol) + ' (' + Syst.Molecule[iMol].Name + ') from HDF5 File')
            Syst = Load_Levels_HDF5(Syst, iMol)

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



def Read_qnsEnBin(Syst, Temp, InputData):
    print('\n  [Read_qnsEnBin]: Reading Level Q.N.s, Energies, Degeneracies and ToBin Mapping')

    for iMol in range(Syst.NMolecules):

        if ( (InputData.Kin.CGQCTFlg) and (InputData.Kin.Groups.Flg[iMol]) ):
            Syst.Molecule[iMol].EqNbBins = Syst.Molecule[iMol].Grouped.NbGroups
        else:
            Syst.Molecule[iMol].EqNbBins = Syst.Molecule[iMol].NBins          
        print('  [Read_qnsEnBin]:   Nb Levels/Groups: ' + str(Syst.Molecule[iMol].EqNbBins))

        PathToFile  = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        f           = h5py.File(PathToFile, 'a')
        TempStr     = Syst.Molecule[iMol].Name + '/Levelg'
        PresentFlg  = TempStr in f.keys()
        f.close()
        if (PresentFlg): print('  [Read_qnsEnBin]:   Found Data for Molecule Nb ' + str(iMol) + ' (' + Syst.Molecule[iMol].Name + ') in the HDF5 File')

        if ( (PresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):
            print('  [Read_qnsEnBin]:   Loading Data for Molecule Nb ' + str(iMol) + ' (' + Syst.Molecule[iMol].Name + ') from HDF5 File')
            Syst = Load_qnsEnBin_HDF5(Syst, iMol)

            for iT in Temp.iTVec:
                TTra = Temp.TranVec[iT-1]
                TInt = TTra
                
                print('  [Read_qnsEnBin]:   Loading Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File')
                Syst = Load_PartFuncsAndEnergiesAtT_HDF5(Syst, iMol, iT, TTra, TInt)

        else:
            if (InputData.OldVersionFlg):
                PathToFile = Syst.PathToFolder + '/' + Syst.Molecule[iMol].Name + '/' + Syst.Molecule[iMol].Name + '_' + str(Syst.Molecule[iMol].EqNbBins) + '/qnsEnBin.dat'
                Data = pandas.read_csv(PathToFile, header=None, skiprows=1, delimiter=r"\s+")
            else:
                PathToFile = Syst.PathToFolder + '/' + Syst.Molecule[iMol].Name + '/Bins_' + str(Syst.Molecule[iMol].EqNbBins) + '/QNsEnBin.csv'
                Data = pandas.read_csv(PathToFile, header=None, skiprows=1)
            print('  [Read_qnsEnBin]:   Reading File: ', PathToFile)
            
            Data = Data.apply(pandas.to_numeric, errors='coerce')

            Syst.Molecule[iMol].Levelvqn   = np.array(Data.values[:,1], dtype=np.int16  )
            Syst.Molecule[iMol].Leveljqn   = np.array(Data.values[:,2], dtype=np.int16  )
            Syst.Molecule[iMol].LevelEEh   = np.array(Data.values[:,3], dtype=np.float64)
            Syst.Molecule[iMol].Levelg     = np.array(Data.values[:,4], dtype=np.float64)
            Syst.Molecule[iMol].LevelToBin = np.array(Data.values[:,5], dtype=np.int32  )

            if (InputData.HDF5.Save_Flg):
                print('  [Read_qnsEnBin]:   Saving Data for Molecule Nb ' + str(iMol) + ' (' + Syst.Molecule[iMol].Name + ') in the HDF5 File')
                Save_qnsEnBin_HDF5(Syst, iMol)

            for iT in Temp.iTVec:
                TTra = Temp.TranVec[iT-1]
                TInt = TTra

                Syst.Molecule[iMol].T[iT-1].LevelEeV  = Syst.Molecule[iMol].LevelEEh * Hartree_To_eV
                Syst.Molecule[iMol].T[iT-1].Q         = Syst.Molecule[iMol].Levelg
                Syst.Molecule[iMol].T[iT-1].QRatio    = Syst.Molecule[iMol].T[iT-1].Q / np.sum(Syst.Molecule[iMol].T[iT-1].Q)

                Syst.Molecule[iMol].T[iT-1].LevelQExp = Syst.Molecule[iMol].Levelg * np.exp( - Syst.Molecule[iMol].T[iT-1].LevelEeV * Ue / (TTra * UKb) )
                Syst.Molecule[iMol].T[iT-1].LevelQExp = Syst.Molecule[iMol].T[iT-1].LevelQExp / np.sum(Syst.Molecule[iMol].T[iT-1].LevelQExp)

                if (InputData.HDF5.Save_Flg):
                    print('  [Read_qnsEnBin]:   Saving Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File')
                    Save_PartFuncsAndEnergiesAtT_HDF5(Syst, iMol, iT, TTra, TInt)

    return Syst



# def Read_PartFuncsAndEnergies(Syst, Temp, InputData):
#     print('\n  [Read_PartFuncsAndEnergies]: Reading Level Partition Functions and Energies')

#     for iMol in range(Syst.NMolecules):

#         if (InputData.Kin.Groups.Flg[iMol]):
#             NBins = Syst.Molecule[iMol].Grouped.NbGroups
#         else:
#             NBins = Syst.Molecule[iMol].NBins
#         print('\n  [Read_PartFuncsAndEnergies]: Nb Levels/Groups: ' + str(NBins))

#         for iT in Temp.iTVec:
#             TTra = Temp.TranVec[iT-1]
#             TInt = TTra

#             PathToFile  = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
#             f           = h5py.File(PathToFile, 'a')
#             TStr        = 'T_' + str(int(TTra)) + '_' + str(int(TInt))
#             TempStr     = TStr + '/' + Syst.Molecule[iMol].Name + '/LevelQ'
#             TPresentFlg = TempStr in f.keys()
#             f.close()
#             if (TPresentFlg): print('  [Read_PartFuncsAndEnergies]: Found Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File')

#             if ( (TPresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):
#                 print('  [Read_PartFuncsAndEnergies]: Loading Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) from HDF5 File')
#                 Load_PartFuncsAndEnergiesAtT_HDF5(Syst, iMol, iT, TTra, TInt)

#             else:
#                 if (InputData.OldVersionFlg):
#                     PathToFile = Syst.PathToFolder + '/' + Syst.Molecule[iMol].Name + '/' + Syst.Molecule[iMol].Name + '_' + str(NBins) + '/T' + str(int(Temp.TranVec[iT-1])) + '.dat'
#                     Data = pandas.read_csv(PathToFile, header=None, skiprows=1, delimiter=r"\s+")
#                 else:
#                     PathToFile = Syst.PathToFolder + '/' + Syst.Molecule[iMol].Name + '/Bins_' + str(NBins) + '/T' + str(int(Temp.TranVec[iT-1])) + '.csv'
#                     Data = pandas.read_csv(PathToFile, header=None, skiprows=1)
#                 print('  [Read_PartFuncsAndEnergies]: Reading File: ', PathToFile)

#                 Data = Data.apply(pandas.to_numeric, errors='coerce')

#                 Syst.Molecule[iMol].T[iT-1].QRatio   = np.array(Data.values[:,0])
#                 Syst.Molecule[iMol].T[iT-1].Q        = np.array(Data.values[:,1])
#                 Syst.Molecule[iMol].T[iT-1].LevelEeV = np.array(Data.values[:,2])


#                 if (InputData.HDF5.Save_Flg):
#                     print('  [Read_PartFuncsAndEnergies]: Saving Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File')
#                     Save_PartFuncsAndEnergiesAtT_HDF5(Syst, iMol, iT, TTra, TInt)

#             Syst.Molecule[iMol].T[iT-1].LevelQExp = Syst.Molecule[iMol].Levelg * np.exp( - Syst.Molecule[iMol].T[iT-1].LevelEeV * Ue / (TTra * UKb) )
#             Syst.Molecule[iMol].T[iT-1].LevelQExp = Syst.Molecule[iMol].T[iT-1].LevelQExp / np.sum(Syst.Molecule[iMol].T[iT-1].LevelQExp)

#     return Syst



def Read_RatesFile_OldVersion(Syst, TTra, TInt, iLevel):
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
            print('      [Read_RatesFile_OldVersion]: Rate File Does Not Contain Rates: ' + PathToFile )
            ProcessesTemp = np.zeros(1, dtype=int)+1
            RatesTemp     = np.zeros(1)
            RatesSDTemp   = np.zeros(1)
    else:
        print('      [Read_RatesFile_OldVersion]: Rate File Does Not Exist: ' + PathToFile )
        ProcessesTemp = np.zeros(1, dtype=int)+1
        RatesTemp     = np.zeros(1)
        RatesSDTemp   = np.zeros(1)

    return ProcessesTemp, RatesTemp, RatesSDTemp



def Read_RatesAtT_OldVersion(Syst, TTra, TInt, iT):

    NProcTot = 1
    for iP in range(3):
        Syst.Pair[iP].NBins    = Syst.Molecule[Syst.Pair[iP].ToMol].NBins
        NProcTot               = NProcTot + Syst.Pair[iP].NBins
        Syst.Pair[iP].NProcTot = NProcTot + AddProc

    Syst.T[iT-1].Proc[0].Rates = np.zeros((Syst.Molecule[0].NBins, 3))
    Syst.T[iT-1].Proc[1].Rates = np.zeros((Syst.Molecule[0].NBins, Syst.Molecule[0].NBins))
    for iProc in range(2, 4):
        Syst.T[iT-1].Proc[iProc].Rates       = np.zeros((Syst.Molecule[0].NBins, Syst.Molecule[Syst.Pair[iProc-1].ToMol].NBins))
    for iExch in range(2, Syst.NProcTypes):
        Syst.T[iT-1].ProcExch[iExch-2].Rates = np.zeros((Syst.Molecule[0].NBins, Syst.Molecule[Syst.ExchtoMol[iExch-2]].NBins))

    i  = 0
    ii = 10
    for iBins in range(Syst.Molecule[0].NBins):
        if (int(iBins/Syst.Molecule[0].NBins*100) == int(i*ii)):
            print('    [Read_RatesAtT_OldVersion]: Read ' + str(i*ii) + '% of the Rate Files')
            i = i+1 
        RatesTempAll                            = np.zeros(Syst.Pair[-1].NProcTot+1)
        [ProcessesTemp, RatesTemp, RatesSDTemp] = Read_RatesFile_OldVersion(Syst, TTra, TInt, iBins)
        RatesTempAll[ProcessesTemp[:]-1]        = RatesTemp[:]

        RatesSplitted                       = np.split( RatesTempAll, np.array([1, Syst.Pair[0].NProcTot, Syst.Pair[1].NProcTot, Syst.Pair[2].NProcTot]) )
        Syst.T[iT-1].Proc[0].Rates[iBins,0] = RatesSplitted[0]
        Syst.T[iT-1].Proc[1].Rates[iBins,:] = RatesSplitted[1]
        for iProc in range(2, 4):
            Syst.T[iT-1].Proc[iProc].Rates[iBins,:] = RatesSplitted[iProc]
    
    for iProc in range(2, 4):
        for iExch in range(2, Syst.NProcTypes):
            if (Syst.Pair[iProc-1].ToMol == Syst.ExchtoMol[iExch-2] ):
                Syst.T[iT-1].ProcExch[iExch-2].Rates = Syst.T[iT-1].ProcExch[iExch-2].Rates + Syst.T[iT-1].Proc[iProc].Rates

    return Syst



def Read_Rates_OldVersion(Syst, Temp, InputData):

    print('\n  [Read_Rates_OldVersion]: Uploading Rates')

    if (InputData.Kin.WindAvrgFlg) and (InputData.Kin.Write_Flg):
        print('\n  [Read_Rates_OldVersion]: Computing Data for Window-Averaging the Kinetics')
        Syst = Compute_WindAvrg_Matrix(Syst, InputData)
    
    for iT in Temp.iTVec:
        TTra = Temp.TranVec[iT-1]
        TInt = TTra

        PathToFile  = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        f           = h5py.File(PathToFile, 'a')
        TStr        = 'T_' + str(int(TTra)) + '_' + str(int(TInt)) + '/Rates/'
        TPresentFlg = TStr in f.keys()
        f.close()
        if (TPresentFlg): print('  [Read_Rates_OldVersion]: Found Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File')

        if ( (TPresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):
            print('  [Read_Rates_OldVersion]: Loading Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) from HDF5 File')
        
            Syst = Load_RatesAtT_HDF5_3Atoms(Syst, TTra, TInt, iT)

        else:
            print('  [Read_Rates_OldVersion]: Reading Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) from set of .dat Files')

            Syst = Read_RatesAtT_OldVersion(Syst, TTra, TInt, iT)

            if (InputData.HDF5.Save_Flg):
                print('  [Read_Rates_OldVersion]: Saving Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File\n')
                Save_RatesAtT_HDF5_3Atoms(Syst, iT, TTra, TInt)
        

        if ((InputData.Rates.PrefJumps_Flg) or (InputData.Kin.Groups.Flg[0])):
            print('  [Read_Rates_OldVersion]: Computing Backweard Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)')
            Syst = Compute_BackwardRates(Syst, iT)


        if (InputData.Kin.Groups.Flg[0]):
            print('  [Read_Rates_OldVersion]: Computing Grouped Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)')
            NbLevels = np.array( [Syst.Molecule[0].NBins]            + [Syst.Molecule[Syst.ExchtoMol[iExch-2]].NBins            for iExch in range(2, Syst.NProcTypes)] ) 
            NbGroups = np.array( [Syst.Molecule[0].Grouped.NbGroups] + [Syst.Molecule[Syst.ExchtoMol[iExch-2]].Grouped.NbGroups for iExch in range(2, Syst.NProcTypes)] )
            Syst.Molecule[0].Grouped.Compute_GroupRates(Syst, iT-1, NbLevels, NbGroups)


        if (InputData.Kin.CorrFactor != 1.0):	
                print('  [Read_Rates_OldVersion]: Correcting Dissociation Rates (Corr Factor = ' + str(InputData.Kin.CorrFactor) + ') for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
                Syst = Compute_Correction_To_DissRates(InputData, Syst, iT)
                if (InputData.Kin.Groups.Flg[0]):
                    print('  [Read_Rates_OldVersion]: Correcting Grouped Dissociation Rates (Corr Factor = ' + str(InputData.Kin.CorrFactor) + ') for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
                    Syst = Compute_Correction_To_GroupedDissRates(InputData, Syst, iT)


        print('  [Read_Rates_OldVersion]: Computing Overall Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
        Syst = Compute_Rates_Overall(InputData, Syst, iT)


        print('  [Read_Rates_OldVersion]: Computing Thermal Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
        Syst = Compute_Rates_Thermal(Syst, iT)


        if (InputData.Rates.PrefJumps_Flg):
            print('  [Read_Rates_OldVersion]: Computing Preferred Jumps for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)')
            Syst = Compute_PrefJumps(InputData, Syst, iT)
            print('  [Read_Rates_OldVersion]: Writing   Preferred Jumps for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)')
            Write_PrefJumps(Syst, Temp, InputData, iT)
            print('\n')


        if (InputData.Kin.Write_Flg):
            print('  [Read_Rates_OldVersion]: Writing Kinetics File for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
            Write_Kinetics(Syst, Temp, InputData, iT)


        if  (InputData.Kin.Groups.Flg[0]) and (InputData.Kin.Groups.Write_Flg):
            print('  [Read_Rates_OldVersion]: Writing Grouped Kinetics File for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
            NbGroups = np.array( [Syst.Molecule[0].Grouped.NbGroups] + [Syst.Molecule[Syst.ExchtoMol[iExch-2]].Grouped.NbGroups for iExch in range(2, Syst.NProcTypes)] )
            Write_GroupedKinetics(Syst, Temp, InputData, iT, NbGroups)


        if (InputData.DelRateMat_Flg):
            for iProc in range(4):
                del Syst.T[iT-1].Proc[iProc].Rates
            for iProc in range(2, Syst.NProcTypes):
                del Syst.T[iT-1].ProcExch[iProc-2].Rates

    print('  [Read_Rates_OldVersion]: Saving Thermal Rates\n')
    Write_Rates_Thermal(Syst, Temp, InputData)

    print('  [Read_Rates_OldVersion]: Plotting Thermal Rates\n')
    Plot_Rates_Thermal(Syst, Temp, InputData)

    return Syst


#########################################################################################################################################################
#########################################################################################################################################################
def Read_RatesFile(Syst, TTra, TInt, iProc):
    #sed -i 's/D-/E-/g' *

    PathToFile = Syst.PathToFolder + '/Rates/T_' + str(int(TTra)) + '_' + str(int(TInt)) + '/Proc' + str(iProc+1) + '.csv'
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
            Data  = pandas.read_csv(PathToFile, header=None, skiprows=5)
            Data  = Data.apply(pandas.to_numeric, errors='coerce')

            ProcessesTemp = np.array(Data[0].values, dtype=np.int64)
            RatesTemp     = np.array(Data[1].values, dtype=np.float64)
            RatesSDTemp   = np.array(Data[2].values, dtype=np.float64)
        else:
            print('      [Read_RatesFile]: Rate File Does Not Contain Rates: ' + PathToFile )
            ProcessesTemp = np.zeros(1, dtype=int)+1
            RatesTemp     = np.zeros(1)
            RatesSDTemp   = np.zeros(1)
    else:
        print('      [Read_RatesFile]: Rate File Does Not Exist: ' + PathToFile )
        ProcessesTemp = np.zeros(1, dtype=int)+1
        RatesTemp     = np.zeros(1)
        RatesSDTemp   = np.zeros(1)

    return ProcessesTemp, RatesTemp, RatesSDTemp



def Read_RatesAtT_3Atoms(InputData, Syst, TTra, TInt, iT):

    for iMol in range(Syst.NMolecules):
        if ( (InputData.Kin.CGQCTFlg) and (InputData.Kin.Groups.Flg[iMol]) ):
            Syst.Molecule[iMol].EqNbBins = Syst.Molecule[iMol].Grouped.NbGroups
        else:
            Syst.Molecule[iMol].EqNbBins = Syst.Molecule[iMol].NBins          

    NBins0 = Syst.Molecule[0].EqNbBins
    print('    [Read_RatesAtT_3Atoms]: Nb of Levels Initial Molecules = ' + str(NBins0) )

    NProcTot = 1
    print('    [Read_RatesAtT_3Atoms]: Nb of Levels Per Final Pair: ')
    for iP in range(3):
        iMol                   = Syst.Pair[iP].ToMol
        Syst.Pair[iP].NBins    = Syst.Molecule[iMol].EqNbBins  
        NProcTot               = NProcTot + Syst.Pair[iP].NBins 
        Syst.Pair[iP].NProcTot = NProcTot
        print('    [Read_RatesAtT_3Atoms]:   Pair ' + str(iP) + ' = ' + str(Syst.Pair[iP].NBins) )


    Syst.T[iT-1].Proc[0].Rates = np.zeros((NBins0, 3))
    Syst.T[iT-1].Proc[1].Rates = np.zeros((NBins0, Syst.Pair[0].NBins))
    for iProc in range(2, 4):
        Syst.T[iT-1].Proc[iProc].Rates       = np.zeros((NBins0, Syst.Pair[iP-1].NBins))
    for iExch in range(2, Syst.NProcTypes):
        Syst.T[iT-1].ProcExch[iExch-2].Rates = np.zeros((NBins0, Syst.Molecule[Syst.ExchtoMol[iExch-2]].EqNbBins))

    i  = 0
    ii = 10
    for iBins in range(Syst.Molecule[0].NBins):
        if (int(iBins/Syst.Molecule[0].NBins*100) == int(i*ii)):
            print('    [Read_RatesAtT_3Atoms]: Read ' + str(i*ii) + '% of the Rate Files')
            i = i+1 
        RatesTempAll                            = np.zeros(Syst.Pair[-1].NProcTot+1)
        [ProcessesTemp, RatesTemp, RatesSDTemp] = Read_RatesFile(Syst, TTra, TInt, iBins)
        RatesTempAll[ProcessesTemp[:]-1]        = RatesTemp[:]

        RatesSplitted                       = np.split( RatesTempAll, np.array([1, Syst.Pair[0].NProcTot, Syst.Pair[1].NProcTot, Syst.Pair[2].NProcTot]) )
        #Syst.T[iT-1].Proc[0].Rates[iBins,0] = RatesSplitted[0]
        Syst.T[iT-1].Proc[0].Rates[iBins,0] = RatesSplitted[1][0:0]
        Syst.T[iT-1].Proc[1].Rates[iBins,:] = RatesSplitted[1][1:]
        for iProc in range(2, 4):
            Syst.T[iT-1].Proc[0].Rates[iBins,iProc-1] = RatesSplitted[iProc][0:0]
            Syst.T[iT-1].Proc[iProc].Rates[iBins,:]   = RatesSplitted[iProc][1:]
    
    for iProc in range(2, 4):
        for iExch in range(2, Syst.NProcTypes):
            if (Syst.Pair[iProc-1].ToMol == Syst.ExchtoMol[iExch-2] ):
                Syst.T[iT-1].ProcExch[iExch-2].Rates = Syst.T[iT-1].ProcExch[iExch-2].Rates + Syst.T[iT-1].Proc[iProc].Rates

    return Syst


def Read_RatesAtT_4Atoms(InputData, Syst, TTra, TInt, iT):

    iPVec    = np.array([0, 1, 2])
    iPOppVec = np.array([5, 4, 3])

    for iMol in range(Syst.NMolecules):
        if ( (InputData.Kin.CGQCTFlg) and (InputData.Kin.Groups.Flg[iMol]) ):
            Syst.Molecule[iMol].EqNbBins = Syst.Molecule[iMol].Grouped.NbGroups
        else:
            Syst.Molecule[iMol].EqNbBins = Syst.Molecule[iMol].NBins          

    NBins0_1 = Syst.Molecule[Syst.Pair[0].ToMol].EqNbBins
    NBins0_2 = Syst.Molecule[Syst.Pair[5].ToMol].EqNbBins
    print('    [Read_RatesAtT_4Atoms]: Nb of Levels in Initial Molecule 1 = ' + str(NBins0_1) )
    print('    [Read_RatesAtT_4Atoms]: Nb of Levels in Initial Molecule 2 = ' + str(NBins0_2) )

    maxNBins = 0
    print('    [Read_RatesAtT_4Atoms]: Nb of Levels Per Final Pair: ')
    for iP in range(6):
        maxNBins               = np.maximum(Syst.Molecule[iMol].EqNbBins, maxNBins)
        iMol                   = Syst.Pair[iP].ToMol
        Syst.Pair[iP].NBins    = Syst.Molecule[iMol].EqNbBins  
        print('    [Read_RatesAtT_4Atoms]:   Pair ' + str(iP) + ' = ' + str(Syst.Pair[iP].NBins) )

    NProcTot = 1
    for iP in range(3):
        NProcTot               = NProcTot + (Syst.Pair[iP].NBins + 1) * (Syst.Pair[iPOppVec[iP]].NBins + 1) 
        Syst.Pair[iP].NProcTot = NProcTot


    Syst.T[iT-1].DissRates     = np.zeros((NBins0_1, NBins0_2, 4))
    Syst.T[iT-1].Proc[0].Rates = np.zeros((NBins0_1, NBins0_2, maxNBins, 6))
    for iP in range(1, 4):
        print('    [Read_RatesAtT_4Atoms]: Pair ' + str(iP) + '; Rate Matrix shape = (' + str(NBins0_1) + '; ' + str(NBins0_2) + '; ' + str(Syst.Pair[iP-1].NBins) + '; ' + str(Syst.Pair[iPOppVec[iP-1]].NBins) + ')')
        Syst.T[iT-1].Proc[iP].Rates          = np.zeros((NBins0_1, NBins0_2, Syst.Pair[iP-1].NBins, Syst.Pair[iPOppVec[iP-1]].NBins))
    for iExch in range(2, Syst.NProcTypes):
        Syst.T[iT-1].ProcExch[iExch-2].Rates = np.zeros((NBins0_1, NBins0_2, Syst.Molecule[Syst.ExchtoMol[iExch-2,0]].EqNbBins, Syst.Molecule[Syst.ExchtoMol[iExch-2,1]].EqNbBins))


    i  = 0
    ii = 10
    for iBins in range(NBins0_1):
        jBinsStart = 0 
        if (Syst.SymmFlg):
            jBinsStart = iBins
        for jBins in range(NBins0_2):
            i = i+1 
            if (int((iBins*jBins)/(NBins0_1*NBins0_2)*100) == int(i*ii)):
                print('    [Read_RatesAtT_4Atoms]: Read ' + str(i*ii) + '% of the Rate Files')

            if (jBins >= jBinsStart):
                RatesTempAll                            = np.zeros(Syst.Pair[2].NProcTot)
                [ProcessesTemp, RatesTemp, RatesSDTemp] = Read_RatesFile(Syst, TTra, TInt, i-1)
                RatesTempAll[ProcessesTemp[:]-1]        = RatesTemp[:]

                iProc = -1

                iProc = iProc + 1
                Syst.T[iT-1].DissRates[iBins, jBins, 0] = Syst.T[iT-1].DissRates[iBins, jBins, 0] + RatesTempAll[iProc]

                for iP in range(1, 4):
                    for kBins in range(Syst.Pair[iP].NBins+1):
                        for lBins in range(Syst.Pair[iPOppVec[iP-1]].NBins+1):
                            iProc = iProc + 1

                            if ( (kBins == 0) and (lBins == 0) ): 
                                Syst.T[iT-1].DissRates[iBins, jBins, iP]                          = Syst.T[iT-1].DissRates[iBins, jBins, iP]                          + RatesTempAll[iProc]
                            elif (kBins == 0):
                                Syst.T[iT-1].Proc[0].Rates[iBins, jBins, lBins-1, iP-1]           = Syst.T[iT-1].Proc[0].Rates[iBins, jBins, lBins-1, iP-1]           + RatesTempAll[iProc]
                            elif (lBins == 0):
                                Syst.T[iT-1].Proc[0].Rates[iBins, jBins, kBins-1, iPOppVec[iP-1]] = Syst.T[iT-1].Proc[0].Rates[iBins, jBins, kBins-1, iPOppVec[iP-1]] + RatesTempAll[iProc]
                            else:
                                iTemp1  = kBins-1
                                iTemp2  = lBins-1
                                SymmFlg = 1.0
                                if (Syst.SymmFlg):
                                    if (iTemp1 > iTemp2):
                                        iTemp3 = iTemp2
                                        iTemp2 = iTemp1
                                        iTemp1 = iTemp3
                                    elif (iTemp1 == iTemp2):
                                        SymmFlg = 2.0
                                Syst.T[iT-1].Proc[iP].Rates[iBins, jBins, iTemp1, iTemp2]         = Syst.T[iT-1].Proc[iP].Rates[iBins, jBins, iTemp1, iTemp2]        + RatesTempAll[iProc] * SymmFlg
            
    for iProc in range(2, 4):
        for iExch in range(2, Syst.NProcTypes):
            if ( ( Syst.MolToCFDComp[Syst.Pair[iProc-1].ToMol] == Syst.MolToCFDComp[Syst.ExchtoMol[iExch-2,0]] ) and ( Syst.MolToCFDComp[Syst.Pair[iPOppVec[iProc-1]].ToMol] == Syst.MolToCFDComp[Syst.ExchtoMol[iExch-2,1]] ) ):
                print('    [Read_RatesAtT_4Atoms]: iPair ' + str(iProc) + ' corresponds to iExch ' + str(iExch-1) ) 
                Syst.T[iT-1].ProcExch[iExch-2].Rates = Syst.T[iT-1].ProcExch[iExch-2].Rates + Syst.T[iT-1].Proc[iProc].Rates

    return Syst


def Read_Rates(Syst, Temp, InputData):

    print('\n  [Read_Rates]: Uploading Rates')

    if (InputData.Kin.WindAvrgFlg) and (InputData.Kin.Write_Flg):
        print('\n  [Read_Rates]: Computing Data for Window-Averaging the Kinetics')
        Syst = Compute_WindAvrg_Matrix(Syst, InputData)
    
    for iT in Temp.iTVec:
        TTra = Temp.TranVec[iT-1]
        TInt = TTra

        PathToFile  = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        f           = h5py.File(PathToFile, 'a')
        TStr        = 'T_' + str(int(TTra)) + '_' + str(int(TInt)) + '/Rates/'
        TPresentFlg = TStr in f.keys()
        f.close()
        if (TPresentFlg): print('  [Read_Rates]: Found Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File')

        if ( (TPresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):
            print('  [Read_Rates]: Loading Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) from HDF5 File')
            
            if (Syst.NAtoms == 3):
                Syst = Load_RatesAtT_HDF5_3Atoms(Syst, TTra, TInt, iT)
            elif (Syst.NAtoms == 4):
                Syst = Load_RatesAtT_HDF5_4Atoms(Syst, TTra, TInt, iT)

        else:
            print('  [Read_Rates]: Reading Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) from set of .dat Files')

            if (Syst.NAtoms == 3):
                Syst = Read_RatesAtT_3Atoms(InputData, Syst, TTra, TInt, iT)
                if (InputData.HDF5.Save_Flg):
                    print('  [Read_Rates]: Saving Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File\n')
                    Save_RatesAtT_HDF5_3Atoms(Syst, iT, TTra, TInt)
            elif (Syst.NAtoms == 4):
                Syst = Read_RatesAtT_4Atoms(InputData, Syst, TTra, TInt, iT)
                if (InputData.HDF5.Save_Flg):
                    print('  [Read_Rates]: Saving Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File\n')
                    Save_RatesAtT_HDF5_4Atoms(Syst, iT, TTra, TInt)
            

        if ((InputData.Rates.PrefJumps_Flg) or ( (not InputData.Kin.CGQCTFlg) and (InputData.Kin.Groups.Flg[0]) ) ):
            print('  [Read_Rates]: Computing Backward Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)')
            Syst = Compute_BackwardRates(Syst, iT)


        if ( (not InputData.Kin.CGQCTFlg) and (InputData.Kin.Groups.Flg[0]) ):
            print('  [Read_Rates]: Computing Grouped Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)')
            NbLevels = np.array( [Syst.Molecule[0].NBins]            + [Syst.Molecule[Syst.ExchtoMol[iExch-2]].NBins            for iExch in range(2, Syst.NProcTypes)] ) 
            NbGroups = np.array( [Syst.Molecule[0].Grouped.NbGroups] + [Syst.Molecule[Syst.ExchtoMol[iExch-2]].Grouped.NbGroups for iExch in range(2, Syst.NProcTypes)] )
            Syst.Molecule[0].Grouped.Compute_GroupRates(Syst, iT-1, NbLevels, NbGroups)


        if (InputData.Kin.CorrFactor != 1.0):   
                print('  [Read_Rates]: Correcting Dissociation Rates (Corr Factor = ' + str(InputData.Kin.CorrFactor) + ') for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
                Syst = Compute_Correction_To_DissRates(InputData, Syst, iT)
                if ( (not InputData.Kin.CGQCTFlg) and (InputData.Kin.Groups.Flg[0]) ):
                    print('  [Read_Rates]: Correcting Grouped Dissociation Rates (Corr Factor = ' + str(InputData.Kin.CorrFactor) + ') for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
                    Syst = Compute_Correction_To_GroupedDissRates(InputData, Syst, iT)


        print('  [Read_Rates]: Computing Overall Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
        Syst = Compute_Rates_Overall(InputData, Syst, iT)


        print('  [Read_Rates]: Computing Thermal Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
        Syst = Compute_Rates_Thermal(Syst, iT)


        if (InputData.Rates.PrefJumps_Flg):
            print('  [Read_Rates]: Computing Preferred Jumps for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)')
            Syst = Compute_PrefJumps(InputData, Syst, iT)
            print('  [Read_Rates]: Writing   Preferred Jumps for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)')
            Write_PrefJumps(Syst, Temp, InputData, iT)
            print('\n')


        if (InputData.Kin.Write_Flg):
            print('  [Read_Rates]: Writing Kinetics File for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
            if (Syst.NAtoms == 3):
                Write_Kinetics_3Atoms(Syst, Temp, InputData, iT)
            elif (Syst.NAtoms == 4):
                Write_Kinetics_4Atoms(Syst, Temp, InputData, iT)

        if ( ( (not InputData.Kin.CGQCTFlg) and (InputData.Kin.Groups.Flg[0]) ) ) and (InputData.Kin.Groups.Write_Flg):
            print('  [Read_Rates]: Writing Grouped Kinetics File for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K)\n')
            NbGroups = np.array( [Syst.Molecule[0].Grouped.NbGroups] + [Syst.Molecule[Syst.ExchtoMol[iExch-2]].Grouped.NbGroups for iExch in range(2, Syst.NProcTypes)] )
            Write_GroupedKinetics(Syst, Temp, InputData, iT, NbGroups)


        if (InputData.DelRateMat_Flg):
            for iProc in range(4):
                del Syst.T[iT-1].Proc[iProc].Rates
            for iProc in range(2, Syst.NProcTypes):
                del Syst.T[iT-1].ProcExch[iProc-2].Rates

    if (Syst.NAtoms == 3):
        print('  [Read_Rates]: Saving Thermal Rates\n')
        Write_Rates_Thermal(Syst, Temp, InputData)

        print('  [Read_Rates]: Plotting Thermal Rates\n')
        Plot_Rates_Thermal(Syst, Temp, InputData)

    return Syst