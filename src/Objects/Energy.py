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
import os.path
from os import path

class energies(object):

    def __init__(self, NTime):

        self.Vib = np.zeros(NTime)
        self.Rot = np.zeros(NTime)
        self.Int = np.zeros(NTime)



    def Compute(self, Syst, ME, iT):

        KQSS_Eps   = 1.e-7

        KDer = np.gradient(Syst.T[iT-1].ProcTot[0].RatesAveraged, ME.Time)
        iQSS = 0
        while (KDer[iQSS] > KQSS_Eps):
            iQSS = iQSS + 1
        iQSS_Start = iQSS - 1
        while (KDer[iQSS] < KQSS_Eps):
            iQSS = iQSS + 1
        iQSS_End = iQSS - 1
        self.iTime[0] = iQSS_Start
        self.iTime[1] = iQSS_End
        self.Time[0]  = ME.Time[iQSS_Start]
        self.Time[1]  = ME.Time[iQSS_End]
        
        self.Rate[0]  = ( Syst.T[iT-1].ProcTot[0].RatesAveraged[iQSS_Start] + Syst.T[iT-1].ProcTot[0].RatesAveraged[iQSS_End] ) / 2.0
        for jProc in range(2, Syst.NProcTypes):
            self.Rate[jProc] = ( Syst.T[iT-1].ProcTot[jProc].RatesAveraged[iQSS_Start] + Syst.T[iT-1].ProcTot[jProc].RatesAveraged[iQSS_End] ) / 2.0



    def Write(self, Syst, Temp, InputData, iT):
 
        PathToFile = InputData.FinalFldr + '/' + InputData.SystNameLong + '_KQSS.csv'
        print('    [Write]: Writing QSS Rates in File: ' + PathToFile )
        with open(PathToFile, 'a') as csvQSS:
            if (not path.exists(PathToFile) ):
                Line       = '# T,KDiss,KInel' 
                for iExch in range(2, Syst.NProcTypes):
                    Line = Line + ',KExch' + str(iExch-1)
                Line = Line + '\n'
                csvQSS.write(Line)
            TempMat = np.transpose( np.expand_dims( np.concatenate( [np.array([Temp.TranVec[iT-1]], float), self.Rate] ), axis=1 ) )
            np.savetxt(csvQSS, TempMat, delimiter=',')
        csvQSS.close()