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
import sys
import os.path
from os import path

from matplotlib import rc 
import matplotlib.pyplot as plt
sys.path.insert(0, '../Parameters/')
from PlotParameters import *

from Atom     import atom
from Molecule import molecule
from Pair     import pair
from CFDComp  import cfdcomp



class component(object):

    def __init__(self, Syst, iComp):

        self.MolFrac    = 0.0

        self.Idx        = iComp
        self.Name       = Syst.CFDComp[iComp].Name
        self.ToMol      = Syst.CFDComp[iComp].ToMol
        self.Mass       = Syst.CFDComp[iComp].Mass
        self.Deg        = Syst.CFDComp[iComp].Deg
        self.Color      = Syst.CFDComp[iComp].Color
        self.LineStyle  = Syst.CFDComp[iComp].LineStyle
        self.RxLxIdx    = Syst.CFDComp[iComp].RxLxIdx
        if ( self.ToMol > -1 ):
            self.Pop   = 0.0
            self.NBins = Syst.Molecule[self.ToMol].NBins



    def Read_Pop(self, InputData, NTime):

        self.Pop = np.zeros((NTime,self.NBins))
        PathToFile = InputData.ME.ReadFldr + '/pop_' + self.Name + '.dat'
        print('\n  [Read_Pop]: Reading Level Populations from the File ' + PathToFile)
        iTime = 0
        for chunk in pandas.read_csv(PathToFile, header=1, chunksize=self.NBins, comment='&', delimiter=r"\s+"):
            Data  = chunk['[eV]'].apply(pandas.to_numeric, errors='coerce')
            self.Pop[iTime,:]  = np.array(Data.values, dtype=np.float64)
            iTime = iTime + 1



    def Compute_DistFunc(self, Syst, NTime):
        print('\n    [Compute_DistFunc]: Computing Distribution Function')

        self.DistFunc = np.zeros((NTime,self.NBins))
        for iTime in range(NTime):
            self.DistFunc[iTime,:] = self.Pop[iTime,:]      * Syst.Molecule[self.ToMol].Levelg[:]
            self.DistFunc[iTime,:] = self.DistFunc[iTime,:] / np.sum( self.DistFunc[iTime,:] )



    def Plot_Pop(self, InputData, Syst, ME, Temp, iT):

        iT=0

        plt.figure()
        plt.title(r'Ro-Vibrational Level Populations', fontsize=FontSize)

        for iTime in (ME.iTimeVec):
            LabelStr = ('t={:.2e}s'.format(ME.Time[iTime]))
            plt.scatter(Syst.Molecule[self.ToMol].T[iT].LevelEeV-Syst.Molecule[self.ToMol].DissEn, self.Pop[iTime,:], label=LabelStr, s=PointSize)
            plt.yscale("log")
        
        plt.xlabel(r'$E_i$ [eV]',          fontsize=FontSize)
        plt.ylabel(r'$log_{10}(N_i/g_i)$', fontsize=FontSize)
        plt.legend(fontsize=20)
        plt.tight_layout()
        if (InputData.PlotShow_Flg):
            plt.show()
        FigSavePath = InputData.FinalFldrT + '/Pops.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_Pop]: Saved Level Population Plot in: ' + FigSavePath)



    def Compute_KAveraged(self, Syst, iT, NTime):
        print('\n    [Compute_KAveraged]: Computing Average Rates')

        for jProc in range(Syst.NProcTypes):
            Syst.T[iT-1].ProcTot[jProc].RatesAveraged = np.zeros(NTime)
            for iTime in range(NTime):
                Syst.T[iT-1].ProcTot[jProc].RatesAveraged[iTime] = np.sum( Syst.T[iT-1].ProcTot[jProc].Rates[:] * self.DistFunc[iTime,:] )

        return Syst



    def Compute_EAveraged(self, Syst, iT, NTime):
        print('\n    [Compute_EAveraged]: Computing Average Energies')

        self.ERot = np.zeros(NTime)
        self.EVib = np.zeros(NTime)
        self.EInt = np.zeros(NTime)
        for iTime in range(NTime):
            self.ERot[iTime] = np.sum( Syst.Molecule[self.ToMol].LevelERot[:] * self.DistFunc[iTime,:] )
            self.EVib[iTime] = np.sum( Syst.Molecule[self.ToMol].LevelEVib[:] * self.DistFunc[iTime,:] )
            self.EInt[iTime] = np.sum( Syst.Molecule[self.ToMol].LevelEeV[:]  * self.DistFunc[iTime,:] )



    def Plot_KAveraged(self, InputData, Syst, ME, Temp, iT):

        plt.figure()
        plt.title(r'$\bar{K}^{Proc}$ Evolution', fontsize=FontSize)

        plt.plot(ME.Time, Syst.T[iT-1].ProcTot[0].RatesAveraged, label=r'$\bar{K}^{D}$')
        plt.plot(Syst.T[iT-1].QSS.Time[0], Syst.T[iT-1].QSS.Rate[0], 'ko')
        plt.plot(Syst.T[iT-1].QSS.Time[1], Syst.T[iT-1].QSS.Rate[0], 'ko')

        for jProc in range(2, Syst.NProcTypes):
            LabelStr = r'$\bar{K}_{' + Syst.Molecule[Syst.ExchtoMol[jProc-2]].Name + '}^{E}$'
            plt.plot(ME.Time, Syst.T[iT-1].ProcTot[jProc].RatesAveraged, label=LabelStr)
            plt.plot(Syst.T[iT-1].QSS.Time[0], Syst.T[iT-1].QSS.Rate[jProc], 'ko')
            plt.plot(Syst.T[iT-1].QSS.Time[1], Syst.T[iT-1].QSS.Rate[jProc], 'ko')
            
        plt.xscale("log")
        plt.yscale("log")
            
        plt.xlabel(r'time [s]', fontsize=FontSize)
        plt.ylabel(r'$\bar{K}^{Proc}$ $[cm^3/s]$', fontsize=FontSize)
        plt.tight_layout()
        plt.legend(fontsize=20)
        if (InputData.PlotShow_Flg):
            plt.show()
        FigSavePath = InputData.FinalFldrT + '/KAveraged_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_TTran_Evolution]: Saved Averaged Rates Evolution Plot in: ' + FigSavePath)



    def Plot_EAveraged(self, InputData, Syst, ME, Temp, iT):

        plt.figure()
        plt.title(r'$\bar{K}^{Proc}$ Evolution', fontsize=FontSize)

        plt.plot(ME.Time, self.ERot, '-.r', label=r'$\bar{E}_{R}$')
        plt.plot(ME.Time, self.EVib, '-k',  label=r'$\bar{E}_{V}$')
        plt.plot(ME.Time, self.EInt, '--b', label=r'$\bar{E}_{I}$')
            
        plt.xscale("log")
            
        plt.xlabel(r'time [s]', fontsize=FontSize)
        plt.ylabel(r'$\bar{E}$ $[eV]$', fontsize=FontSize)
        plt.tight_layout()
        plt.legend(fontsize=20)
        if (InputData.PlotShow_Flg):
            plt.show()
        FigSavePath = InputData.FinalFldrT + '/EAveraged_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_EAveraged]: Saved Averaged Energy Evolution Plot in: ' + FigSavePath)



    def Compute_Taus(self, InputData, Syst, ME, iT):

        iComp    = Syst.ColPartToComp
        Pressure = ME.Component[iComp].MolFrac[-1] * ME.P[-1] / 101325.0
        print('\n    [Compute_Taus]: Partial Pressure for iComp ' + str(iComp+1) + ' = ' + str(Pressure*100.0) + "%")

        EVibLim = (self.EVib[-1] - self.EVib[0]) * 0.632 + self.EVib[0]
        iVib=0
        while (self.EVib[iVib] < EVibLim):
          iVib = iVib+1;
        DeltaT =   ME.Time[iVib] -   ME.Time[iVib-1]
        DeltaE = self.EVib[iVib] - self.EVib[iVib-1]
        SemiE  =         EVibLim - self.EVib[iVib-1]
        self.TauVib = DeltaT / DeltaE * SemiE + ME.Time[iVib-1]
        self.TauVib = Pressure * self.TauVib
        print('    [Compute_Taus]: Vibrational Tau = ' + str(self.TauVib))

        ERotLim = (self.ERot[-1] - self.ERot[0]) * 0.632 + self.ERot[0]
        iRot=0
        while (self.ERot[iRot] < ERotLim):
          iRot = iRot+1;
        DeltaT =   ME.Time[iRot] -   ME.Time[iRot-1]
        DeltaE = self.ERot[iRot] - self.ERot[iRot-1]
        SemiE  =         ERotLim - self.ERot[iRot-1]
        self.TauRot = DeltaT / DeltaE * SemiE + ME.Time[iRot-1]
        self.TauRot = Pressure * self.TauRot
        print('    [Compute_Taus]: Rotational Tau = ' + str(self.TauRot))

        EIntLim = (self.EInt[-1] - self.EInt[0]) * 0.632 + self.EInt[0]
        iInt=0
        while (self.EInt[iInt] < EIntLim):
          iInt = iInt+1;
        DeltaT =   ME.Time[iInt] -   ME.Time[iInt-1]
        DeltaE = self.EInt[iInt] - self.EInt[iInt-1]
        SemiE  =         EIntLim - self.EInt[iInt-1]
        self.TauInt = DeltaT / DeltaE * SemiE + ME.Time[iInt-1]
        self.TauInt = Pressure * self.TauInt
        print('    [Compute_Taus]: Internal Tau = ' + str(self.TauInt))


    def Write_Taus(self, Temp, InputData, iT):
 
        PathToFile = InputData.FinalFldr + '/Taus.csv'
        print('\n    [Write_Taus]: Writing Taus in File: ' + PathToFile )
        if (not path.exists(PathToFile) ):
            WriteFlg = True
        else:
            WriteFlg = False        
        with open(PathToFile, 'a') as csvTaus:
            if (WriteFlg):
                Line       = '# T,RotTau,VibTau,IntTau\n'
                csvTaus.write(Line)
            TempVec = np.array([self.TauRot, self.TauVib, self.TauInt])
            TempMat = np.transpose( np.expand_dims( np.concatenate( [np.array([Temp.TranVec[iT-1]], float), TempVec] ), axis=1 ) )
            np.savetxt(csvTaus, TempMat, delimiter=',')
        csvTaus.close()



class me_output(object):

    def __init__(self, Syst):

        self.NCFDComp  = Syst.NCFDComp
        self.Component = [component(Syst, iComp) for iComp in range(self.NCFDComp)]



    def Read_Box(self, InputData):

        PathToFile = InputData.ME.ReadFldr + 'box.dat'
        print('\n  [Read_Box]: Reading Master Equation Solution from the File ' + PathToFile)
        Data  = pandas.read_csv(PathToFile, header=None, delimiter=r"\s+")
        Data  = Data.apply(pandas.to_numeric, errors='coerce')

        self.Time  = np.array(Data[0].values, dtype=np.float64)
        self.NTime = self.Time.shape[0]
        print('  [Read_Box]: Found ' + str(self.NTime+1) + ' Time Instants')

        self.TimeVec  = InputData.ME.TimeVec
        self.iTimeVec = np.zeros((self.TimeVec.size), dtype=np.int32)
        jT = 0
        for Time in self.TimeVec:
            iTime = 0
            while ( (iTime < self.NTime-1) and (self.Time[iTime] < Time) ):
                iTime = iTime + 1
            self.iTimeVec[jT] = iTime - 1
            jT = jT + 1
        print('  [Read_Box]: Vector of Times ' + str(self.TimeVec) + ' corresponds to Vector of Time Instants ' + str(self.iTimeVec) + '; (' + str(self.Time[self.iTimeVec[:]]) + ')')

        for iComp in range(self.NCFDComp):
            self.Component[iComp].MolFrac = np.array(Data[1+iComp].values, dtype=np.float64)

        self.TTran = np.array(Data[1+self.NCFDComp].values, dtype=np.float64)
        self.Rho   = np.array(Data[2+self.NCFDComp].values, dtype=np.float64)
        self.P     = np.array(Data[3+self.NCFDComp].values, dtype=np.float64)
        self.Nd    = np.array(Data[4+self.NCFDComp].values, dtype=np.float64)
        self.EEh   = np.array(Data[5+self.NCFDComp].values, dtype=np.float64)



    def Plot_MolFracs_Evolution(self, InputData, Temp, iT):

        #rc.style.use('seaborn')
        plt.figure()
        plt.title(r'Mole Fraction Evolutions', fontsize=FontSize)

        # for iComp in range(self.NCFDComp):
        #     LabelStr = self.Component[iComp].Name
        #     TempStr  = 'C' + str(iComp+1)
        #     plt.plot(self.Time, self.Component[iComp].MolFrac, linestyle=self.Component[iComp].LineStyle, color=self.Component[iComp].Color, label=LabelStr, linewidth=LineWidth)
        #     #plt.plot(self.Time, self.Component[iComp].MolFrac, TempStr, linestyle=self.Component[iComp].LineStyle, label=LabelStr)

        iComp    = 1        
        LabelStr = self.Component[iComp].Name
        TempStr  = 'C' + str(iComp+1)
        plt.plot(self.Time, self.Component[iComp].MolFrac, linestyle=self.Component[iComp].LineStyle, color=self.Component[iComp].Color, label=LabelStr, linewidth=LineWidth)
        #plt.plot(self.Time, self.Component[iComp].MolFrac, TempStr, linestyle=self.Component[iComp].LineStyle, label=LabelStr)
        
        plt.xscale("log")
        
        plt.xlabel(r'time [s]', fontsize=FontSize)
        plt.ylabel(r'Mole Fraction', fontsize=FontSize)
        plt.legend(fontsize=20)
        plt.tight_layout()
        if (InputData.PlotShow_Flg):
            plt.show()
        FigSavePath = InputData.FinalFldrT + '/MoleFracs_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_MolFracs_Evolution]: Saved Mole Fraction Evolutions Plot in: ' + FigSavePath)



    def Plot_TTran_Evolution(self, InputData, Temp, iT):

        plt.figure()
        plt.title(r'$T_{Tran}$ Evolution', fontsize=FontSize)

        LabelStr = '$T_{Tran}$'
        plt.plot(self.Time, self.TTran, '-k', label=LabelStr)
        plt.xscale("log")
        
        plt.xlabel(r'time [s]', fontsize=FontSize)
        plt.ylabel(r'T [K]', fontsize=FontSize)
        plt.tight_layout()
        if (InputData.PlotShow_Flg):
            plt.show()
        FigSavePath = InputData.FinalFldrT + '/T_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_TTran_Evolution]: Saved Temperature Evolution Plot in: ' + FigSavePath)



    def Plot_Rho_Evolution(self, InputData, Temp, iT):

        plt.figure()
        plt.title(r'Density Evolution', fontsize=FontSize)

        LabelStr = '$\rho$'
        plt.plot(self.Time, self.Rho, '-k', label=LabelStr)
        plt.xscale("log")
        
        plt.xlabel(r'time [s]', fontsize=FontSize)
        plt.ylabel(r'\rho $[Kg/m^3]$', fontsize=FontSize)
        plt.tight_layout()
        if (InputData.PlotShow_Flg):
            plt.show()
        FigSavePath = InputData.FinalFldrT + '/Rho_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_Rho_Evolution]: Saved Density Evolution Plot in: ' + FigSavePath)
        


    def Plot_P_Evolution(self, InputData, Temp, iT):

        plt.figure()
        plt.title(r'Pressure Evolution', fontsize=FontSize)

        LabelStr = 'P'
        plt.plot(self.Time, self.P, '-k', label=LabelStr)
        plt.xscale("log")
        
        plt.xlabel(r'time [s]', fontsize=FontSize)
        plt.ylabel(r'P [Pa]', fontsize=FontSize)
        plt.tight_layout()
        if (InputData.PlotShow_Flg):
            plt.show()
        FigSavePath = InputData.FinalFldrT + '/P_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_P_Evolution]: Saved Pressure Evolution Plot in: ' + FigSavePath)



    def Plot_Nd_Evolution(self, InputData, Temp, iT):

        plt.figure()
        plt.title(r'Number Density Evolution', fontsize=FontSize)

        LabelStr = 'Nd'
        plt.plot(self.Time, self.Nd, '-k', label=LabelStr)
        plt.xscale("log")
        
        plt.xlabel(r'time [s]', fontsize=FontSize)
        plt.ylabel(r'$N_D$ $[mol/m^3]$', fontsize=FontSize)
        plt.tight_layout()
        if (InputData.PlotShow_Flg):
            plt.show()
        FigSavePath = InputData.FinalFldrT + '/Nd_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_Nd_Evolution]: Saved Number Density Evolution Plot in: ' + FigSavePath)



    def Plot_Energy_Evolution(self, InputData, Temp, iT):

        plt.figure()
        plt.title(r'Energy Evolution', fontsize=FontSize)

        LabelStr = '$E_{Int}$'
        plt.plot(self.Time, self.EEh, '-k', label=LabelStr)
        plt.xscale("log")
        
        plt.xlabel(r'time [s]', fontsize=FontSize)
        plt.ylabel(r'E [Eh]', fontsize=FontSize)
        plt.tight_layout()
        if (InputData.PlotShow_Flg):
            plt.show()
        FigSavePath = InputData.FinalFldrT + '/EEh_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_Energy_Evolution]: Saved Energy Evolution Plot in: ' + FigSavePath)
