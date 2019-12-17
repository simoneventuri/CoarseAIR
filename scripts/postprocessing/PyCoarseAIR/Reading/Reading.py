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

def Read_PartFuncsAndEnergies(Syst, Temp):

    for iMol in range(Syst.NMolecules):

        for iTra in range(Temp.NTran):

            PathToFile = Syst.PathToFolder + '/' + Syst.Molecule[iMol].Name + '/' + Syst.Molecule[iMol].Name + '_' + str(Syst.Molecule[iMol].NBins) + '/T' + str(int(Temp.TranVec[iTra])) + '.dat'
            print('Reading Data From File: ', PathToFile)

            Data = pandas.read_csv(PathToFile, header=None, skiprows=1, delimiter=r"\s+")
            Data = Data.apply(pandas.to_numeric, errors='coerce')

            Syst.Molecule[iMol].T[iTra].QRatio = np.array(Data.values[:,0])
            Syst.Molecule[iMol].T[iTra].Q      = np.array(Data.values[:,1])
            Syst.Molecule[iMol].LevelEeV       = np.array(Data.values[:,2])

    return Syst



def Read_qnsEnBin(Syst):

    for iMol in range(Syst.NMolecules):

        PathToFile = Syst.PathToFolder + '/' + Syst.Molecule[iMol].Name + '/' + Syst.Molecule[iMol].Name + '_' + str(Syst.Molecule[iMol].NBins) + '/qnsEnBin.dat'
        print('Reading Data From File: ', PathToFile)
        
        Data = pandas.read_csv(PathToFile, header=None, skiprows=1, delimiter=r"\s+")
        Data = Data.apply(pandas.to_numeric, errors='coerce')

        Syst.Molecule[iMol].Levelvqn   = np.array(Data.values[:,1], dtype=np.int16  )
        Syst.Molecule[iMol].Leveljqn   = np.array(Data.values[:,2], dtype=np.int16  )
        Syst.Molecule[iMol].LevelEEh   = np.array(Data.values[:,3], dtype=np.float64)
        Syst.Molecule[iMol].Levelg     = np.array(Data.values[:,4], dtype=np.float64)
        Syst.Molecule[iMol].LevelToBin = np.array(Data.values[:,5], dtype=np.int32  )

    return Syst



def Read_RatesFile_CGQCT(Syst, TTra, TInt, iLevel):

    PathToFile = Syst.PathToFolder + '/' + Syst.Molecule[0].Name + '/Rates/T_' + str(int(TTra)) + '_' + str(int(TInt)) + '/Bin' + str(iLevel+1) + '.dat'
    #print(PathToFile)
    
    Data  = pandas.read_csv(PathToFile, header=None, skiprows=5, delimiter=r"\s+")
    Data  = Data.apply(pandas.to_numeric, errors='coerce')

    ProcessesTemp = np.array(Data[1].values, dtype=np.int64)
    RatesTemp     = np.array(Data[2].values, dtype=np.float64)
    RatesSDTemp   = np.array(Data[3].values, dtype=np.float64)

    return ProcessesTemp, RatesTemp, RatesSDTemp



def Read_Rates_CGQCT(Syst, Temp):

    for iT in Temp.iTVec:
        TTra = Temp.TranVec[iT-1]
        TInt = TTra
        print('\nTemperature Nb ', iT, '; T = ', TTra, 'K')

        Syst.T[iT-1].Proc[0].Rates = np.zeros((Syst.Molecule[0].NBins, 3))
        Syst.T[iT-1].Proc[1].Rates = np.zeros((Syst.Molecule[0].NBins, Syst.Molecule[0].NBins))
        for iProc in range(2, 4):
            Syst.T[iT-1].Proc[iProc].Rates     = np.zeros((Syst.Molecule[0].NBins, Syst.Molecule[Syst.Pair[iProc-1].ToMol].NBins))
        for iProc in range(2, Syst.NProcTypes):
            Syst.T[iT-1].ProcExch[iProc-2].Rates = np.zeros((Syst.Molecule[0].NBins, Syst.Molecule[Syst.ExchtoMol[iProc-2]].NBins))

        for iBins in range(10):#range(Syst.Molecule[0].NBins):
            print('\nBins Nb ', iBins+1)
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
                #print('iP = ', iProc, '; ToMol     = ', Syst.Pair[iProc-1].ToMol+1)
                #print('ExchtoMol = ', Syst.ExchtoMol[jProc-2]+1 )
                if (Syst.Pair[iProc-1].ToMol == Syst.ExchtoMol[jProc-2] ):
                    Syst.T[iT-1].ProcExch[jProc-2].Rates = Syst.T[iT-1].ProcExch[jProc-2].Rates + Syst.T[iT-1].Proc[iProc].Rates

        Syst.T[iT-1].ProcTot[0].Rates = np.sum(Syst.T[iT-1].Proc[0].Rates, axis=1)
        Syst.T[iT-1].ProcTot[1].Rates = np.sum(Syst.T[iT-1].Proc[1].Rates, axis=1)
        for jProc in range(2, Syst.NProcTypes):
            Syst.T[iT-1].ProcTot[jProc].Rates = np.sum(Syst.T[iT-1].ProcExch[jProc-2].Rates, axis=1)

    return Syst