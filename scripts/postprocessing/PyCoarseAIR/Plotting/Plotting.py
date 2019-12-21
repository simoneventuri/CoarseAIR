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

from matplotlib import rc 
import matplotlib.pyplot as plt
sys.path.insert(0, '../Parameters/')
from PlotParameters import *


def Plot_Rates_Thermal(Syst, Temp, InputData):

    plt.figure()
    plt.title(r'$K_{i}^{Th}$', fontsize=FontSize)

    LabelStr = '$K_{D}^{Th}$ for ' + Syst.Molecule[0].Name + ', Interpolated'
    plt.plot(10000/Temp.TranVec, Syst.RatesTh[:,0], '-k', label=LabelStr)
    
    LabelStr = '$K_{D}^{Th}$ for ' + Syst.Molecule[0].Name
    plt.plot(10000/Temp.TranVec, Syst.RatesTh[:,0], 'ok')#, label=LabelStr)
    
    for jProc in range(2, Syst.NProcTypes):

        LabelStr = '$K_{Ex' + str(jProc-1) + '}^{Th}$ for ' + Syst.Molecule[0].Name + ', Interpolated'
        plt.plot(10000/Temp.TranVec, Syst.RatesTh[:,jProc], '-', label=LabelStr)
        
        LabelStr = '$K_{Ex' + str(jProc-1) + '}^{Th}$ for ' + Syst.Molecule[0].Name
        plt.plot(10000/Temp.TranVec, Syst.RatesTh[:,jProc], 'ok')#, label=LabelStr)

    plt.legend(fontsize=FontSize)
    plt.yscale("log")
    plt.xlabel(r'10,000/T [1/K]',             fontsize=FontSize)
    plt.ylabel(r'$K_{Proc}^{Th}$ $[cm^3/s]$', fontsize=FontSize)
    plt.tight_layout()
    if (InputData.PlotShow_Flg):
        plt.show()
    FigSavePath = InputData.FinalFldr + '/' + InputData.SystNameLong + '_KTh.png'
    plt.savefig(FigSavePath)
    print('    [Plot_Rates_Thermal]: Saved Thermal Rates Plot in: ' + FigSavePath)



def Plot_DissRates_Thermal(Syst, Temp, InputData):

    plt.figure()
    plt.title(r'$K_{i}^{D}$', fontsize=FontSize)

    LabelStr = '$K_{i}^{D}$ for ' + Syst.Molecule[0].Name + ', Interpolated'
    plt.plot(10000/Temp.TranVec, Syst.RatesTh[:,0], '-k', label=LabelStr)
    
    LabelStr = '$K_{i}^{D}$ for ' + Syst.Molecule[0].Name
    plt.plot(10000/Temp.TranVec, Syst.RatesTh[:,0], 'ok', label=LabelStr)
    plt.yscale("log")
    
    plt.xlabel(r'10,000/T [1/K]',         fontsize=FontSize)
    plt.ylabel(r'$K_{i}^{D}$ $[cm^3/s]$', fontsize=FontSize)
    plt.tight_layout()
    if (InputData.PlotShow_Flg):
        plt.show()
    FigSavePath = InputData.FinalFldr + '/' + InputData.SystNameLong + '_KTh_Diss.png'
    plt.savefig(FigSavePath)
    print('    [Plot_DissRates_Thermal]: Saved Dissociation Thermal Rates Plot in: ' + FigSavePath)
