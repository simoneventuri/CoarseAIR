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

from Parameters import UKb, Ue

def Compute_ThermalDissRates_FromOverall(Syst, Temp, InputData):

    Syst.KDissTh = np.zeros(Temp.NTran)
    for iT in Temp.iTVec:

        LevelKDiss = np.zeros(Syst.Molecule[0].NBins)

        DissFile = InputData.ReadKinFolder   + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name + '.csv'
        with open(DissFile) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',')
            next(readCSV)
            for row in readCSV:
                if (float(row[iT]) > 0.0):
                    LevelKDiss[int(row[0])-1] = float(row[iT])



        Syst.Molecule[0].LevelQ = Syst.Molecule[0].Levelg * np.exp( - Syst.Molecule[0].LevelEeV * Ue / (Temp.TranVec[iT-1] * UKb) )
        Syst.Molecule[0].LevelQ = Syst.Molecule[0].LevelQ / np.sum(Syst.Molecule[0].LevelQ)

        Syst.KDissTh[iT-1] = sum( Syst.Molecule[0].LevelQ * LevelKDiss )
    
    print('Thermal Dissociation Rates = ', Syst.KDissTh)

    PathToFile = InputData.PostprocessingFldr + '/' + Syst.Molecule[0].Name + '_KDiss_Thermal.csv'
    with open(PathToFile, 'w') as f:
        TStr = str(Temp.TranVec[Temp.iTVec[0]-1])
        KStr = str(Syst.KDissTh[Temp.iTVec[0]-1])
        for iT in Temp.iTVec[1:-1]:
            TStr = TStr + ',' + str(Temp.TranVec[Temp.iTVec[iT]-1])
            KStr = KStr + ',' + str(Syst.KDissTh[Temp.iTVec[iT]-1])
        TStr = TStr + '\n'
        KStr = KStr + '\n'
        Line = '# Dissociation Thermal Rates (Line 2) @ Multiple Translational Temperatures (Line 1)\n'
        f.write(Line)
        f.write(TStr)
        f.write(KStr)

    plt.figure()
    plt.title(r'$K_{i}^{D}$', fontsize=20)
    LabelStr = '$K_{i}^{D}$ for ' + Syst.Molecule[0].Name + ', Interpolated'
    plt.plot(10000/Temp.TranVec, Syst.KDissTh, '-k', label=LabelStr)
    LabelStr = '$K_{i}^{D}$ for ' + Syst.Molecule[0].Name
    plt.plot(10000/Temp.TranVec, Syst.KDissTh, 'ok', label=LabelStr)
    plt.yscale("log")
    plt.xlabel(r'10,000/T [1/K]', fontsize=20)
    plt.ylabel(r'$K_{i}^{D} [cm^3/s]$', fontsize=20)
    if (InputData.PlotShowFlg):
        plt.show()
    FigSavePath = InputData.PostprocessingFldr + '/' + Syst.Molecule[0].Name + '_KDiss_Thermal.png'
    plt.savefig(FigSavePath)