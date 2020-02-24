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
import scipy
import scipy.signal as scsig
from scipy.optimize import curve_fit
from os import path


from matplotlib import rc 
import matplotlib.pyplot as plt

sys.path.insert(0, '../Plotting/')

from Plotting      import Plot_DissRates_Thermal
from Writing       import Write_Arrhenius_Diss, Write_Arrhenius_Inel, Write_Arrhenius_Exch
from Saving        import Save_Arr_Diss_HDF5, Save_Arr_Inel_HDF5, Save_Arr_Exch_HDF5


def Compute_Correction_To_DissRates(InputData, Syst, iT):

    Syst.T[iT-1].Proc[0].Rates[:,0] = Syst.T[iT-1].Proc[0].Rates[:,0] * InputData.Kin.CorrFactor

    return Syst


def Compute_Correction_To_GroupedDissRates(InputData, Syst, iT):

    Syst.Molecule[0].Grouped.T[iT].Proc[0].Rates =  Syst.Molecule[0].Grouped.T[iT].Proc[0].Rates * InputData.Kin.CorrFactor

    return Syst


def Compute_BackwardRates(Syst, iT):

    Syst.T[iT-1].Proc[1].BckRates = Syst.T[iT-1].Proc[1].Rates
    for iLevel in range(Syst.Molecule[0].NBins):
        for jLevel in range(Syst.Molecule[0].NBins):
            if (Syst.Molecule[0].LevelEeV[iLevel] < Syst.Molecule[0].LevelEeV[jLevel]):
                Syst.T[iT-1].Proc[1].BckRates[iLevel, jLevel] = Syst.T[iT-1].Proc[1].Rates[jLevel, iLevel] * Syst.Molecule[0].T[iT-1].LevelQExp[jLevel] / Syst.Molecule[0].T[iT-1].LevelQExp[iLevel]


    for iProc in range(2, Syst.NProcTypes):
        jMol = Syst.ExchtoMol[iProc-2]

        Syst.T[iT-1].ProcExch[iProc-2].BckRates = Syst.T[iT-1].ProcExch[iProc-2].Rates
        for iLevel in range(Syst.Molecule[0].NBins):
            for jLevel in range(Syst.Molecule[jMol].NBins):
                if (Syst.Molecule[0].LevelEeV[iLevel] < Syst.Molecule[jMol].LevelEeV[jLevel]):
                    Syst.T[iT-1].ProcExch[iProc-2].BckRates[iLevel, jLevel] = Syst.T[iT-1].ProcExch[iProc-2].Rates[jLevel, iLevel] * Syst.Molecule[0].T[iT-1].LevelQExp[jLevel] / Syst.Molecule[jMol].T[iT-1].LevelQExp[iLevel]
                
    return Syst


def Compute_Rates_Overall(Syst, iT):

    Syst.T[iT-1].ProcTot[0].Rates = np.sum(Syst.T[iT-1].Proc[0].Rates, axis=1)
    Syst.T[iT-1].ProcTot[1].Rates = np.sum(Syst.T[iT-1].Proc[1].Rates, axis=1)
    for jProc in range(2, Syst.NProcTypes):
        Syst.T[iT-1].ProcTot[jProc].Rates = np.sum(Syst.T[iT-1].ProcExch[jProc-2].Rates, axis=1)

    return Syst


def Compute_Rates_Thermal(Syst, iT):

    for jProc in range(Syst.NProcTypes):
        Syst.RatesTh[iT-1,jProc] = sum( Syst.Molecule[0].T[iT-1].LevelQExp * Syst.T[iT-1].ProcTot[jProc].Rates )

    return Syst


def Compute_Rates_Thermal_FromOverall(Syst, Temp, InputData):

    DissFile = InputData.Kin.ReadFldr   + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name + '.csv'
    print('  [Compute_Rates_Thermal_FromOverall]: Reading Dissociation Rates From File: ' + DissFile)
    for iT in Temp.iTVec:
        LevelKDiss = np.zeros(Syst.Molecule[0].NBins)
        with open(DissFile) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',')
            next(readCSV)
            for row in readCSV:
                if (float(row[iT]) > 0.0):
                    LevelKDiss[int(row[0])-1] = float(row[iT])
        csvfile.close()

        Syst.RatesTh[iT-1,0] = sum( Syst.Molecule[0].T[iT-1].LevelQExp * LevelKDiss )
    
    print('  [Compute_Rates_Thermal_FromOverall]: Thermal Dissociation Rates = ', Syst.T[iT-1].ProcTot[0].Rates)
    Write_DissRates_Thermal(Syst, Temp, InputData)

    print('  [Compute_Rates_Thermal_FromOverall]: Plotting Dissociation Thermal Rates')
    Plot_DissRates_Thermal(Syst, Temp, InputData)

    return Syst


def Compute_QSS(Syst, ME, iT):

    KQSS_Eps = 1.e-9
    EpsT     = 1.e-4

    #yy   = scsig.savgol_filter(Syst.T[iT-1].ProcTot[0].RatesAveraged, 21, 3)
    yy   = Syst.T[iT-1].ProcTot[0].RatesAveraged
    #KDer = np.gradient(yy, ME.Time)
    # print(KDer)
    # plt.figure()
    # plt.plot(ME.Time, Syst.T[iT-1].ProcTot[0].RatesAveraged)
    # plt.xscale("log")
    # plt.yscale("log")

    # plt.figure()
    # plt.plot(ME.Time, KDer)
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.show()
    iQSS = 0
    while (yy[iQSS] <= yy[0]*1.e1):
        iQSS = iQSS + 1
    while ( np.abs(np.log10(yy[iQSS]) - np.log10(yy[iQSS+1])) / np.abs(np.log10(yy[iQSS])) > EpsT ) and ( iQSS < yy.size-3):
        iQSS = iQSS + 1
    iQSS_Start = iQSS - 1
    while ( np.abs(np.log10(yy[iQSS]) - np.log10(yy[iQSS+1])) / np.abs(np.log10(yy[iQSS])) <= EpsT ) and ( iQSS < yy.size-3):
        iQSS = iQSS + 1
    iQSS_End = iQSS - 1

    if (iQSS < yy.size-2):
        Syst.T[iT-1].QSS.iTime[0] = iQSS_Start
        Syst.T[iT-1].QSS.iTime[1] = iQSS_End
        Syst.T[iT-1].QSS.Time[0]  = ME.Time[iQSS_Start]
        Syst.T[iT-1].QSS.Time[1]  = ME.Time[iQSS_End]
        print('QSS Starting Time: t = ', Syst.T[iT-1].QSS.Time[0], 's\n')
        print('QSS Final    Time: t = ', Syst.T[iT-1].QSS.Time[1], 's\n')
        Syst.T[iT-1].QSS.Rate[0]  = ( Syst.T[iT-1].ProcTot[0].RatesAveraged[iQSS_Start] + Syst.T[iT-1].ProcTot[0].RatesAveraged[iQSS_End] ) / 2.0
        for jProc in range(2, Syst.NProcTypes):
            Syst.T[iT-1].QSS.Rate[jProc] = ( Syst.T[iT-1].ProcTot[jProc].RatesAveraged[iQSS_Start] + Syst.T[iT-1].ProcTot[jProc].RatesAveraged[iQSS_End] ) / 2.0
    else:
        print('QSS Not Found!\n')
        iQSS_Start                = KDer.size
        iQSS_End                  = KDer.size
        Syst.T[iT-1].QSS.Time[0]  = ME.Time[iQSS_Start-1]
        Syst.T[iT-1].QSS.Time[1]  = ME.Time[iQSS_End-1]
        Syst.T[iT-1].QSS.Rate[0]  = 0.0
        for jProc in range(2, Syst.NProcTypes):
            Syst.T[iT-1].QSS.Rate[jProc] = 0.0

    return Syst


def Compute_PrefJumps(InputData, Syst, iT):

    NJumps = InputData.Rates.NPrefJumps

    Syst.T[iT-1].Proc[1].PrefJumps = np.zeros((Syst.Molecule[0].NBins, NJumps), dtype=np.int32)
    for iLevel in range(Syst.Molecule[0].NBins):
        TempVec = Syst.T[iT-1].Proc[1].BckRates[iLevel, :]
        Syst.T[iT-1].Proc[1].PrefJumps[iLevel,:] = np.argsort(TempVec)[-NJumps:] + 1

    for iProc in range(2, Syst.NProcTypes):
        jMol = Syst.ExchtoMol[iProc-2]
        
        Syst.T[iT-1].ProcExch[iProc-2].PrefJumps = np.zeros((Syst.Molecule[0].NBins, NJumps), dtype=np.int32)
        for iLevel in range(Syst.Molecule[0].NBins):
            TempVec =  Syst.T[iT-1].ProcExch[iProc-2].BckRates[iLevel, :]
            Syst.T[iT-1].ProcExch[iProc-2].PrefJumps[iLevel,:] = np.argsort(TempVec)[-NJumps:] + 1

    return Syst


def Read_DissRates_For_Arrhenius(Syst, Temp, InputData):

    AllDissMat = np.zeros((Syst.Molecule[0].NBins,Temp.NTran))

    for iT in Temp.iTVec:
        TTra = Temp.TranVec[iT-1]
        TInt = TTra

        PathToFile  = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        f           = h5py.File(PathToFile, 'a')
        TStr        = 'T_' + str(int(TTra)) + '_' + str(int(TInt)) + '/Rates/'
        TPresentFlg = TStr in f.keys()
        f.close()
        if (TPresentFlg): print('  [Read_DissRates_For_Arrhenius]: Found Dissociation Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(TTra) + 'K) in the HDF5 File')

        if ( (TPresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):
            print('  [Read_DissRates_For_Arrhenius]: Loading Dissociation Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(TTra) + 'K) from HDF5 File')
        
            Load_DissRatesAtT_HDF5(Syst, TTra, TInt, iT)

        else:
            print('  [Read_DissRates_For_Arrhenius]: HDF5 File does not contain Dissociation Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(TTra) + 'K)')

        AllDissMat[:,iT-1] = Syst.T[iT-1].Proc[0].Rates[:,0] * InputData.Kin.CorrFactor
        print('  [Read_DissRates_For_Arrhenius]: Updated AllDissMat with Dissociation Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(TTra) + 'K)')

    return AllDissMat


def Read_InelRates_For_Arrhenius_CGQCT(Syst, Temp, InputData, iLevel):

    RatesInel  = np.zeros((Temp.NTran, Syst.Molecule[0].NBins))
    RatesExch  = np.zeros((Temp.NTran, Syst.Molecule[0].NBins))

    for iT in Temp.iTVec:
        TTra = Temp.TranVec[iT-1]
        TInt = TTra

        #print(       InputData.QCTOutFldr + '/' + Syst.Name + '/' + Syst.Molecule[0].Name + '/Rates/T_' + str(int(TTra)) + '_' + str(int(TInt)) + '/Bin' + str(iLevel+1) + '.dat' )
        PathToFile = InputData.QCTOutFldr + '/' + Syst.Name + '/' + Syst.Molecule[0].Name + '/Rates/T_' + str(int(TTra)) + '_' + str(int(TInt)) + '/Bin' + str(iLevel+1) + '.dat'
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
                print('      [Read_InelRates_For_Arrhenius_CGQCT]: Rate File Does Not Contain Rates: ' + PathToFile )
                ProcessesTemp = np.zeros(1, dtype=int)+1
                RatesTemp     = np.zeros(1)
                RatesSDTemp   = np.zeros(1)
        else:
            print('      [Read_InelRates_For_Arrhenius_CGQCT]: Rate File Does Not Exist: ' + PathToFile )
            ProcessesTemp = np.zeros(1, dtype=int)+1
            RatesTemp     = np.zeros(1)
            RatesSDTemp   = np.zeros(1)

        RatesTempAll                     = np.zeros(Syst.Pair[-1].NProcTot+1)
        RatesTempAll[ProcessesTemp[:]-1] = RatesTemp[:]

        RatesSplitted     = np.split( RatesTempAll, np.array([1, Syst.Pair[0].NProcTot, Syst.Pair[1].NProcTot, Syst.Pair[2].NProcTot]) )
        Temp1             = RatesSplitted[0]
        RatesInel[iT-1,:] = RatesSplitted[1]
        for iProc in range(2, 4):
            RatesExch[iT-1,:] = RatesExch[iT-1,:] + RatesSplitted[iProc]

    return np.transpose(RatesInel), np.transpose(RatesExch)


def ArrhFunc(x, a, b, c):
    return a + b*np.log(x) - c/x


def ArrhEval(x, a, b, c):
    return a * x**b * np.exp(c/x)


def Compute_ArrheniusFit(TVec, RatesVec):
    
    #Fit for the parameters a, b, c of the function ArrhFunc:
    popt, pcov = curve_fit(ArrhFunc, TVec, np.log(RatesVec), p0=[1.0, 10.0, 100.0], absolute_sigma=True)
    #perr       = np.sqrt(np.diag(pcov))
    #print(perr)
    #popt, pcov = curve_fit(ArrhFunc, TVec, RatesVec)

    # plt.plot(TVec, np.log10(RatesVec), 'b-', label='data')
    # plt.plot(TVec, np.log10(np.exp(ArrhFunc(TVec, *popt))), 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    # plt.xlabel('T [K]')
    # plt.ylabel('log10(KDiss) [cm^3/s]')
    # plt.legend()
    # plt.show()

    return popt


def Compute_ArrheniusError(TVec, RatesVec, Coeff):
    
    PredVec  = ArrhFunc(TVec, Coeff[0], Coeff[1], Coeff[2])
    ErrorVec = PredVec - np.log(np.maximum(RatesVec,1.e-20)) 
    Error    = np.amax( ErrorVec )
    return Error


def Compute_ArrheniusFit_Diss(Syst, Temp, InputData):

    AllDissMat = Read_DissRates_For_Arrhenius(Syst, Temp, InputData)

    if (InputData.Kin.WriteArr_Flg):
        DissKinetics = Write_Arrhenius_Diss(Syst, InputData, 0, 0, InputData.Kin.MaxEntOrPlato, 0)
    if (InputData.Kin.SaveArr_Flg):
        Coeffs = []
        Idxs   = []
    for iLevel in range(Syst.Molecule[0].NBins):
        RatesVec = AllDissMat[iLevel,:]
        TVec     = Temp.TranVec
        Mask     = np.ma.masked_equal(RatesVec,0.0).mask
        if (np.any(Mask)):
            TVec     = TVec[np.invert(Mask)]
            RatesVec = RatesVec[np.invert(Mask)]
        if (TVec.size >= Syst.Arr.MinNbTs):
            if (TVec.size == 1):
                TVec     = np.array([     0.5*TVec[0],     TVec[0],      2.0*TVec[0]])
                RatesVec = np.array([RatesVec[0]/10.0, RatesVec[0], RatesVec[0]*10.0])
            elif (TVec.size == 2):
                TVec     = np.array([         TVec[0],     TVec[1],      2.0*TVec[1]])
                RatesVec = np.array([     RatesVec[0], RatesVec[1], RatesVec[1]*10.0])
            CoeffsTemp    = Compute_ArrheniusFit(TVec, RatesVec)
            CoeffsTemp[0] = np.exp(CoeffsTemp[0])
            if (InputData.Kin.WriteArr_Flg):
                DissKinetics = Write_Arrhenius_Diss(Syst, InputData, iLevel, CoeffsTemp, InputData.Kin.MaxEntOrPlato, DissKinetics)
            if (InputData.Kin.SaveArr_Flg):
                Coeffs.append(CoeffsTemp)
                Idxs.append(iLevel+1)
        else:
            print("\n    [Compute_ArrheniusFit_Diss]: WARNING: Level " + str(iLevel) + " does not have at least " + str(Syst.Arr.MinNbTs) + " Dissociation Rates higher than " + str(Syst.Arr.MinRate) + " [cm^3/s]!" )
    
    if (InputData.Kin.WriteArr_Flg):
        DissKinetics.close()

    if (InputData.Kin.SaveArr_Flg):
        Syst.Arr.Proc[0].Coeffs = np.array(Coeffs               )
        Syst.Arr.Proc[0].Idxs   = np.array(Idxs,  dtype=np.int64)
        Save_Arr_Diss_HDF5(Syst, 0)

    return Syst


def Compute_ArrheniusFit_Inel(Syst, Temp, InputData):

    if (InputData.Kin.WriteArr_Flg):
        InelKinetics = Write_Arrhenius_Inel(Syst, InputData, 0, 0, InputData.Kin.MaxEntOrPlato, 0, 0)
        if (InputData.Kin.MaxEntOrPlato == 1):
            ExchKinetics = Write_Arrhenius_Inel(Syst, InputData, 0, 0, 1, 1, 0)
        else:
            ExchKinetics = Write_Arrhenius_Exch(Syst, InputData, 0, 0, 2, 1, 0)
    if (InputData.Kin.SaveArr_Flg):
        InelCoeffs = []
        InelIdxs   = []
        ExchCoeffs = []
        ExchIdxs   = []

        for iLevel in range(Syst.Molecule[0].NBins):
            print("\n    [Compute_ArrheniusFit_Inel]: Computing Arrhenius for Inelastic Rates starting from Level Nb. " + str(iLevel+1) )

            [InelRates, ExchRates] = Read_InelRates_For_Arrhenius_CGQCT(Syst, Temp, InputData, iLevel)
            AllRates               = InelRates  + ExchRates
            
            for jLevel in range(0,iLevel):
                
                if (InputData.Kin.Exch_Flg):

                    InelRatesVecOrig = InelRates[jLevel,:]
                    TVecOrig         = Temp.TranVec
                    InelRatesVec     = InelRatesVecOrig
                    TVec             = TVecOrig
                    Mask         = np.ma.masked_equal(InelRatesVec,0.0).mask
                    if (np.any(Mask)):
                        TVec     = TVec[np.invert(Mask)]
                        InelRatesVec = InelRatesVec[np.invert(Mask)]
                    if (TVec.size >= Syst.Arr.MinNbTs):
                        if (TVec.size == 1):
                            TVec     = np.array([             0.5*TVec[0],         TVec[0],          2.0*TVec[0]])
                            InelRatesVec = np.array([InelRatesVec[0]/10.0, InelRatesVec[0], InelRatesVec[0]*10.0])
                        elif (TVec.size == 2):
                            TVec     = np.array([             TVec[0],         TVec[1],          2.0*TVec[1]])
                            RatesVec = np.array([     InelRatesVec[0], InelRatesVec[1], InelRatesVec[1]*10.0])
                        InelCoeffsTemp    = Compute_ArrheniusFit(TVec, InelRatesVec)
                        ArrErr            = Compute_ArrheniusError(TVecOrig, InelRatesVecOrig, InelCoeffsTemp)
                        if (ArrErr <= np.absolute(np.log(InputData.Kin.MaxErrArr))):
                            InelCoeffsTemp[0] = np.exp(InelCoeffsTemp[0])
                            if (InputData.Kin.WriteArr_Flg):
                                InelKinetics = Write_Arrhenius_Inel(Syst, InputData, np.array([iLevel, jLevel]), InelCoeffsTemp, InputData.Kin.MaxEntOrPlato, 0, InelKinetics)
                            if (InputData.Kin.SaveArr_Flg):
                                InelCoeffs.append(InelCoeffsTemp)
                                InelIdxs.append([iLevel,jLevel])
                        #else:
                            #print("\n    [Compute_ArrheniusFit_Inel]: WARNING: Level " + str(iLevel+1) + " to Level " + str(jLevel+1) + " Resulted in too large Fitting Error for Inelastic Rates!" )
                    #else:
                        #print("\n    [Compute_ArrheniusFit_Inel]: WARNING: Level " + str(iLevel+1) + " to Level " + str(jLevel+1) + " does not have at least " + str(Syst.Arr.MinNbTs) + " Inelastic Rates higher than " + str(Syst.Arr.MinRate) + " [cm^3/s]!" )


                if (InputData.Kin.Inel_Flg):

                    if (InputData.Kin.MaxEntOrPlato == 1):
                        InelRatesVec = AllRates[jLevel,:]
                        TVec         = Temp.TranVec
                        Mask         = np.ma.masked_equal(InelRatesVec,0.0).mask
                        if (np.any(Mask)):
                            TVec     = TVec[np.invert(Mask)]
                            InelRatesVec = InelRatesVec[np.invert(Mask)]
                        if (TVec.size >= Syst.Arr.MinNbTs):
                            if (TVec.size == 1):
                                TVec     = np.array([             0.5*TVec[0],         TVec[0],          2.0*TVec[0]])
                                InelRatesVec = np.array([InelRatesVec[0]/10.0, InelRatesVec[0], InelRatesVec[0]*10.0])
                            elif (TVec.size == 2):
                                TVec     = np.array([             TVec[0],         TVec[1],          2.0*TVec[1]])
                                RatesVec = np.array([     InelRatesVec[0], InelRatesVec[1], InelRatesVec[1]*10.0])
                            InelCoeffsTemp    = Compute_ArrheniusFit(TVec, InelRatesVec)
                            ArrErr            = Compute_ArrheniusError(TVecOrig, InelRatesVecOrig, InelCoeffsTemp)
                            if (ArrErr <= np.absolute(np.log(InputData.Kin.MaxErrArr))):
                                InelCoeffsTemp[0] = np.exp(InelCoeffsTemp[0])
                                if (InputData.Kin.WriteArr_Flg):
                                    ExchKinetics = Write_Arrhenius_Inel(Syst, InputData, np.array([iLevel, jLevel]), InelCoeffsTemp, 1, 1, ExchKinetics)
                                if (InputData.Kin.SaveArr_Flg):
                                    ExchCoeffs.append(InelCoeffsTemp)
                                    ExchIdxs.append([iLevel,jLevel])
                            #else:
                                #print("\n    [Compute_ArrheniusFit_Inel]: WARNING: Level " + str(iLevel+1) + " to Level " + str(jLevel+1) + " Resulted in too large Fitting Error for Inelastic+Exhange Rates!" )
                        #else:
                            #print("\n    [Compute_ArrheniusFit_Inel]: WARNING: Level " + str(iLevel+1) + " to Level " + str(jLevel+1) + " does not have at least " + str(Syst.Arr.MinNbTs) + " Inelastic+Exhange Rates higher than " + str(Syst.Arr.MinRate) + " [cm^3/s]!" )


                    if (InputData.Kin.MaxEntOrPlato == 2):
                        InelRatesVec = ExchRates[jLevel,:]
                        TVec         = Temp.TranVec
                        Mask         = np.ma.masked_equal(InelRatesVec,0.0).mask
                        if (np.any(Mask)):
                            TVec     = TVec[np.invert(Mask)]
                            InelRatesVec = InelRatesVec[np.invert(Mask)]
                        if (TVec.size >= Syst.Arr.MinNbTs):
                            if (TVec.size == 1):
                                TVec     = np.array([             0.5*TVec[0],         TVec[0],          2.0*TVec[0]])
                                InelRatesVec = np.array([InelRatesVec[0]/10.0, InelRatesVec[0], InelRatesVec[0]*10.0])
                            elif (TVec.size == 2):
                                TVec     = np.array([             TVec[0],         TVec[1],          2.0*TVec[1]])
                                RatesVec = np.array([     InelRatesVec[0], InelRatesVec[1], InelRatesVec[1]*10.0])
                            InelCoeffsTemp    = Compute_ArrheniusFit(TVec, InelRatesVec)
                            ArrErr            = Compute_ArrheniusError(TVecOrig, InelRatesVecOrig, InelCoeffsTemp)
                            if (ArrErr <= np.absolute(np.log(InputData.Kin.MaxErrArr))):
                                InelCoeffsTemp[0] = np.exp(InelCoeffsTemp[0])
                                if (InputData.Kin.WriteArr_Flg):
                                    InelKinetics = Write_Arrhenius_Exch(Syst, InputData, np.array([iLevel, jLevel]), InelCoeffsTemp, 2, 1, InelKinetics)
                                if (InputData.Kin.SaveArr_Flg):
                                    ExchCoeffs.append(InelCoeffsTemp)
                                    ExchIdxs.append([iLevel,jLevel])
                            #else:
                                #print("\n    [Compute_ArrheniusFit_Inel]: WARNING: Level " + str(iLevel+1) + " to Level " + str(jLevel+1) + " Resulted in too large Fitting Error for Inelastic+Exhange Rates!" )
                        #else:
                            #print("\n    [Compute_ArrheniusFit_Inel]: WARNING: Level " + str(iLevel+1) + " to Level " + str(jLevel+1) + " does not have at least " + str(Syst.Arr.MinNbTs) + " Inelastic+Exhange Rates higher than " + str(Syst.Arr.MinRate) + " [cm^3/s]!" )
                    

        if (InputData.Kin.WriteArr_Flg):
            if (InputData.Kin.Inel_Flg):
                InelKinetics.close()
            if (InputData.Kin.Exch_Flg):
                ExchKinetics.close()

        if (InputData.Kin.SaveArr_Flg):
            if (InputData.Kin.Inel_Flg):
                Syst.Arr.Proc[1].Coeffs = np.array(InelCoeffs               )
                Syst.Arr.Proc[1].Idxs   = np.array(InelIdxs,  dtype=np.int64) + 1
                Save_Arr_Inel_HDF5(Syst, 0)
            if (InputData.Kin.Exch_Flg):
                Syst.Arr.Proc[2].Coeffs = np.array(ExchCoeffs               )
                Syst.Arr.Proc[2].Idxs   = np.array(ExchIdxs,  dtype=np.int64) + 1
                Save_Arr_Exch_HDF5(Syst, 0)

    return Syst



def Compute_WindAvrg_Matrix(Syst, InputData):

    print('        [Compute_WindAvrg_Matrix]: Computing Data for Kinetics Window-Averaging')
    Syst.Molecule[0].WindAvrgDet   = int( (2*InputData.Kin.WindAvrgJs+1) * (2*InputData.Kin.WindAvrgVs+1) )
    Syst.Molecule[0].WindAvrgMat   = np.zeros( (Syst.Molecule[0].NBins,  Syst.Molecule[0].WindAvrgDet), dtype=np.int64)
    Syst.Molecule[0].WindAvrgFound = np.zeros( (Syst.Molecule[0].NBins,  1), dtype=np.int64 )

    for iLevel in range(Syst.Molecule[0].NBins):
        iFound = -1
        for jLevel in range(Syst.Molecule[0].NBins):
            if (np.absolute(Syst.Molecule[0].Levelvqn[iLevel] - Syst.Molecule[0].Levelvqn[jLevel]) <= InputData.Kin.WindAvrgVs) and  (np.absolute(Syst.Molecule[0].Leveljqn[iLevel] - Syst.Molecule[0].Leveljqn[jLevel]) <= InputData.Kin.WindAvrgJs):
                iFound     = iFound  + 1
                Syst.Molecule[0].WindAvrgMat[iLevel,iFound] = jLevel
        Syst.Molecule[0].WindAvrgFound[iLevel]              = iFound 
    
    print('        [Compute_WindAvrg_Matrix]: Done Computing Data for Kinetics Window-Averaging\n')

    return Syst