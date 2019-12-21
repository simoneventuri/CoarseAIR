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

        for iComp in range(self.NCFDComp):
            self.Component[iComp].MolFrac = np.array(Data[1+iComp].values, dtype=np.float64)

        self.TTran = np.array(Data[1+self.NCFDComp].values, dtype=np.float64)
        self.Rho   = np.array(Data[2+self.NCFDComp].values, dtype=np.float64)
        self.P     = np.array(Data[3+self.NCFDComp].values, dtype=np.float64)
        self.Nd    = np.array(Data[4+self.NCFDComp].values, dtype=np.float64)
        self.EEh   = np.array(Data[5+self.NCFDComp].values, dtype=np.float64)



    def Plot_MolFracs_Evolution(self, InputData):

        #rc.style.use('seaborn')
        plt.figure()
        plt.title(r'Mole Fraction Evolutions', fontsize=FontSize)

        for iComp in range(self.NCFDComp):
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
        FigSavePath = InputData.FinalFldr + '/' + InputData.SystNameLong + '_MoleFracs_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_MolFracs_Evolution]: Saved Mole Fraction Evolutions Plot in: ' + FigSavePath)



    def Plot_TTran_Evolution(self, InputData):

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
        FigSavePath = InputData.FinalFldr + '/' + InputData.SystNameLong + '_T_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_TTran_Evolution]: Saved Temperature Evolution Plot in: ' + FigSavePath)



    def Plot_Rho_Evolution(self, InputData):

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
        FigSavePath = InputData.FinalFldr + '/' + InputData.SystNameLong + '_Rho_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_Rho_Evolution]: Saved Density Evolution Plot in: ' + FigSavePath)
        


    def Plot_P_Evolution(self, InputData):

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
        FigSavePath = InputData.FinalFldr + '/' + InputData.SystNameLong + '_P_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_P_Evolution]: Saved Pressure Evolution Plot in: ' + FigSavePath)



    def Plot_Nd_Evolution(self, InputData):

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
        FigSavePath = InputData.FinalFldr + '/' + InputData.SystNameLong + '_Nd_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_Nd_Evolution]: Saved Number Density Evolution Plot in: ' + FigSavePath)



    def Plot_Energy_Evolution(self, InputData):

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
        FigSavePath = InputData.FinalFldr + '/' + InputData.SystNameLong + '_EEh_Evo.png'
        plt.savefig(FigSavePath)
        print('\n    [Plot_Energy_Evolution]: Saved Energy Evolution Plot in: ' + FigSavePath)