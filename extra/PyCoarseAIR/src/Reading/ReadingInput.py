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

from System       import system
from Temperatures import temperatures

def (PathToOutput)

    FileName = PathToOutput + '/InputForBash.inp'
    with open(FileName) as InputForBash:


        iLine    = 1
        SystName = InputForBash.readline()


        ## Reading 
        Temp = temperatures()

        iLine   += 1
        Temp.TranFlg  = InputForBash.readline()

        iLine   += 1
        Temp.NTran    = int(InputForBash.readline())
        Temp.TranVec  = np.empty(Temp.NTran, float)
        for iTran in range(Temp.NTran):
            iLine += 1
            Temp.TranVec[iTran] = float(InputForBash.readline())

        iLine   += 1
        Temp.NInt   = int(InputForBash.readline())
        Temp.IntVec = np.empty(Temp.NInt, float)
        if (Temp.NInt > 1):
            for iInt in range(Temp.NInt):
                iLine += 1
                Temp.IntVec[iInt] = float(InputForBash.readline())


        iLine   += 1
        NMolecules = int(InputForBash.readline())
        Molecule   = [molecule() for iMol in range(NMolecules)]
        for iMol in range(NMolecules):
            
            iLine += 1
            Molecule[iMol].Name    = InputForBash.readline()
            
            iLine += 1
            Molecule[iMol].KinMthd = InputForBash.readline()

            if (Molecule[iMol].KinMthd == 'Sta'):
               Molecule[iMol].KinMthd = 'STS'

               iLine += 1
               Molecule[iMol].StartBin = int(InputForBash.readline())
               
               iLine += 1
               Molecule[iMol].FinalBin = int(InputForBash.readline())
               
               iLine += 1

            elif (Molecule[iMol].KinMthd == 'Vib'):
               Molecule[iMol].KinMthd  = 'VIB'

               iLine += 1
               Molecule[iMol].StartBin = int(InputForBash.readline())

               iLine += 1
               Molecule[iMol].FinalBin = int(InputForBash.readline())

               iLine += 1

            elif ( (Molecule[iMol].KinMthd == 'RoV') or (Molecule[iMol].KinMthd == 'Fro') ):
               Molecule[iMol].KinMthd = 'CGM'

               iLine += 1
               Molecule[iMol].StartBin = int(InputForBash.readline())

               iLine += 1
               Molecule[iMol].FinalBin = int(InputForBash.readline())

               iLine += 1
               Molecule[iMol].NBins    = int(InputForBash.readline())


    return SystName, Temp, NMolecules, Molecule