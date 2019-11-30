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


def ReadPartFuncsAndEnergies(Syst, Temp):

    for iMol in range(Syst.NMolecules):

        for iTra in range(Temp.NTra):

            PathToFile = Syst.PathToFolder + '/' + Syst.Molecule[iMol].Name + '/' + Syst.Molecule[iMol].Name + '_' + str(Syst.Molecule[iMol].NBins) + '/T' + str(int(Temp.TraVec[iTra])) + '.dat'
            print('Reading Data From File: ', PathToFile)

            Data = pandas.read_csv(PathToFile, header=None, skiprows=1, delimiter=r"\s+")
            Data = Data.apply(pandas.to_numeric, errors='coerce')

            Syst.Molecule[iMol].T[iTra].QRatio = np.array([Data.values[:,0]])
            Syst.Molecule[iMol].T[iTra].Q      = np.array([Data.values[:,1]])
            Syst.Molecule[iMol].EeV            = np.array([Data.values[:,2]])


    return Syst