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
import os, errno
import os.path
from os import path
import shutil

import numpy as np
import pandas
import h5py
import csv
from   matplotlib import rc 
import matplotlib.pyplot as plt

sys.path.insert(0, '../src/Parameters/')
from Parameters     import *
from PlotParameters import *

from Atom      import atom
from Molecule  import molecule
from Pair      import pair
from CFDComp   import cfdcomp
from Processes import processes
#from QSS       import qss

def mkdirs(newdir, mode=0o777):
    os.makedirs(newdir, mode, exist_ok=True)
    




# ===================================================================================================================== CLASS ===
class t_properties(object):

    def __init__(self, NProcTypes):

        self.Proc     = [processes() for iProc in range(4)]
        self.ProcExch = [processes() for iProc in range(NProcTypes-2)]
        self.ProcTot  = [processes() for iProc in range(NProcTypes)]
        #self.QSS      = qss(NProcTypes)



    # ***************************************************************************************************************************
    def Load_RatesAtT_HDF5( self, Syst ):

        if (Syst.NAtoms == 3):
            print('    [System.py - Load_RatesAtT_HDF5]: Loading Rates for a System of 3 Atoms' )
        
            PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
            HDF5Exist_Flg = path.exists(PathToFile)
            if (HDF5Exist_Flg):
                f = h5py.File(PathToFile, 'a')
            else:
                f = {'key': 'value'}

            TStr = 'T_' + str(int(self.TTra)) + '_' + str(int(self.TInt)) + '/Rates/'
            grp  = f[TStr]

            Data               = grp["Diss"]
            self.Proc[0].Rates = Data[...]
            Data               = grp["Inel"]
            self.Proc[1].Rates = Data[...]

            for iProc in range(2, Syst.NProcTypes):
                ExchStr                      = "Exch_" + str(iProc-1)
                Data                         = grp[ExchStr]
                self.ProcExch[iProc-2].Rates = Data[...]

            f.close()
            
        elif (Syst.NAtoms == 4):
            print('    [System.py - Load_RatesAtT_HDF5]: Loading Rates for a System of 4 Atoms' )

            PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
            f          = h5py.File(PathToFile, "r")

            TStr = 'T_' + str(int(self.TTra)) + '_' + str(int(self.TInt)) + '/Rates/'
            grp  = f[TStr]

            Data                       = grp["Diss"]
            self.DissRates     = Data[...]
            Data                       = grp["DissInel"]
            self.Proc[0].Rates = Data[...]
            Data                       = grp["Inel"]
            self.Proc[1].Rates = Data[...]

            for iProc in range(2, Syst.NProcTypes):
                ExchStr                              = "Exch_" + str(iProc-1)
                Data                                 = grp[ExchStr]
                self.ProcExch[iProc-2].Rates = Data[...]

            f.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Read_RatesAtT_N3( self, InputData, Syst ):       

        if (Syst.NAtoms == 3):
            print('    [System.py - Read_RatesAtT_N3]: Reading Rates for a System of 3 Atoms' )

            NStates0 = Syst.EqNStatesIn[0]
            print('    [System.py - Read_RatesAtT_N3]: Nb of Levels Initial Molecules = ' + str(NStates0) )

            NProcTot = 1
            print('    [System.py - Read_RatesAtT_N3]: Nb of Levels Per Final Pair: ')
            for iP in range(3):
                iMol                   = Syst.Pair[iP].ToMol
                Syst.Pair[iP].NStates  = Syst.EqNStatesIn[iMol]  
                NProcTot               = NProcTot + Syst.Pair[iP].NStates 
                Syst.Pair[iP].NProcTot = NProcTot
                print('    [System.py - Read_RatesAtT_N3]:   Pair ' + str(iP) + ' = ' + str(Syst.Pair[iP].NStates) )


            self.Proc[0].Rates = np.zeros((NStates0, 3))
            self.Proc[1].Rates = np.zeros((NStates0, Syst.Pair[0].NStates))
            for iProc in range(2, 4):
                self.Proc[iProc].Rates       = np.zeros((NStates0, Syst.Pair[iP-1].NStates))


            kT = 1
            for TTemp in Syst.N3Data_TVec:
                if (TTemp == self.TTra):
                    jT = kT
                kT = kT + 1 
            iTemp = jT - 1
            jTemp = jT + 1


            print('    [System.py - Read_RatesAtT_N3]: Reading the ' + str(iTemp+1) + '-th column in the Dissociation Rates File. It corresponds to ' + str(self.TTra) + ' K' )
            self.Proc[0].Rates[Syst.N3Mapping[:],0] = Syst.N3Data_Diss[:,iTemp]


            print('    [System.py - Read_RatesAtT_N3]: Reading the ' + str(jTemp+1) + '-th column in the Inelastic Rates File. It corresponds to ' + str(self.TTra) + ' K' )
            iLevelTemp = np.array(Syst.N3Data_Inel[0].values,     dtype=np.int64) - 1
            jLevelTemp = np.array(Syst.N3Data_Inel[1].values,     dtype=np.int64) - 1
            RatesTemp  = np.array(Syst.N3Data_Inel[jTemp].values, dtype=np.float64)
            for iProc in range( np.size(iLevelTemp) ):
                iiLevel = iLevelTemp[iProc]
                jjLevel = jLevelTemp[iProc]
                iLevel  = iiLevel #Syst.N3Mapping[iiLevel]
                jLevel  = jjLevel #Syst.N3Mapping[jjLevel]
                #print(iiLevel,jjLevel,iLevel+1,jLevel+1)
                self.Proc[1].Rates[iLevel,jLevel] = RatesTemp[iProc]


        if (Syst.NAtoms == 4):
            print('    [System.py - Read_RatesAtT_N3]: ERROR! Read_RatesAtT_N3 Only Inplemented for a System of 3 Atoms' )
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Read_RatesFile_OldVersion( self, Syst, iLevel ):
    #sed -i 's/D-/E-/g' *

        PathToFile = Syst.PathToReadFldr + '/' + Syst.Molecule[0].Name + '/Rates/T_' + str(int(self.TTra)) + '_' + str(int(self.TInt)) + '/Bin' + str(iLevel+1) + '.dat'
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
                print('      [System.py - Read_RatesFile_OldVersion]: WARNING! Rate File Does Not Contain Rates: ' + PathToFile )
                ProcessesTemp = np.zeros(1, dtype=int)+1
                RatesTemp     = np.zeros(1)
                RatesSDTemp   = np.zeros(1)
        else:
            print('      [System.py - Read_RatesFile_OldVersion]: WARNING! Rate File Does Not Exist: ' + PathToFile )
            ProcessesTemp = np.zeros(1, dtype=int)+1
            RatesTemp     = np.zeros(1)
            RatesSDTemp   = np.zeros(1)

        return ProcessesTemp, RatesTemp, RatesSDTemp
    # ...........................................................................................................................


    
    # ***************************************************************************************************************************
    def Read_RatesAtT_OldVersion( self, InputData, Syst ):

        NProcTot = 1
        for iP in range(3):
            Syst.Pair[iP].NLevels  = Syst.Molecule[Syst.Pair[iP].ToMol].NLevels
            NProcTot               = NProcTot + Syst.Pair[iP].NLevels
            Syst.Pair[iP].NProcTot = NProcTot
            print('    [System.py - Read_RatesAtT_OldVersion]: Pair ' + str(iP) + ' has ' + str(Syst.Pair[iP].NProcTot) + ' Nb of Processes' )
            print('    [System.py - Read_RatesAtT_OldVersion]: Total Nb of Processes Up to Pair Nb ' + str(iP) + ' = ' + str(Syst.Pair[iP].NLevels) )
        print('    [System.py - Read_RatesAtT_OldVersion]:   Total Nb of Processes = ' + str(NProcTot) )

        self.Proc[0].Rates = np.zeros((Syst.Molecule[0].NLevels, Syst.NProcTypes))
        self.Proc[1].Rates = np.zeros((Syst.Molecule[0].NLevels, Syst.Molecule[0].NLevels))
        for iProc in range(2, 4):
            self.Proc[iProc].Rates       = np.zeros((Syst.Molecule[0].NLevels, Syst.Molecule[Syst.Pair[iProc-1].ToMol].NLevels))
        for iExch in range(2, Syst.NProcTypes):
            self.ProcExch[iExch-2].Rates = np.zeros((Syst.Molecule[0].NLevels, Syst.Molecule[Syst.ExchtoMol[iExch-2]].NLevels))
        print('    [System.py - Read_RatesAtT_OldVersion]:     Allocated Zero-Tensors of Processes')

        i  = 0
        ii = 10
        for iStates in range(Syst.Molecule[0].NLevels):
            if (int(iStates/Syst.Molecule[0].NLevels*100) == int(i*ii)):
                print('    [System.py - Read_RatesAtT_OldVersion]:       Read ' + str(i*ii) + '% of the Rate Files')
                i = i+1 
            if ( (iStates >= InputData.Kin.MinStateIn[0] - 1) and (iStates <= InputData.Kin.MaxStateIn[0] - 1) ):
                RatesTempAll                            = np.zeros(Syst.Pair[-1].NProcTot+1)
                [ProcessesTemp, RatesTemp, RatesSDTemp] = self.Read_RatesFile_OldVersion( Syst, iStates )
                RatesTempAll[ProcessesTemp[:]-1]        = RatesTemp[:]

                RatesSplitted                 = np.split( RatesTempAll, np.array([1, Syst.Pair[0].NProcTot, Syst.Pair[1].NProcTot, Syst.Pair[2].NProcTot]) )
                self.Proc[0].Rates[iStates,0] = RatesSplitted[0]
                self.Proc[1].Rates[iStates,:] = RatesSplitted[1]
                for iProc in range(2, 4):
                    self.Proc[iProc].Rates[iStates,:] = RatesSplitted[iProc]
            
        for iProc in range(2, 4):
            for iExch in range(2, Syst.NProcTypes):
                if (Syst.Pair[iProc-1].ToMol == Syst.ExchtoMol[iExch-2] ):
                    self.ProcExch[iExch-2].Rates = self.ProcExch[iExch-2].Rates + self.Proc[iProc].Rates
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Read_RatesFile( self, Syst, iProc ):
        #sed -i 's/D-/E-/g' *

        PathToFile = Syst.PathToReadFldr + '/Rates/T_' + str(int(self.TTra)) + '_' + str(int(self.TInt)) + '/Proc' + str(iProc+1) + '.csv'
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

                ProcessesTemp = np.maximum(ProcessesTemp, 1)
            else:
                print('      [System.py - Read_RatesFile]: Rate File Does Not Contain Rates: ' + PathToFile )
                ProcessesTemp = np.zeros(1, dtype=int)+1
                RatesTemp     = np.zeros(1)
                RatesSDTemp   = np.zeros(1)
        else:
            print('      [System.py - Read_RatesFile]: Rate File Does Not Exist: ' + PathToFile )
            ProcessesTemp = np.zeros(1, dtype=int)+1
            RatesTemp     = np.zeros(1)
            RatesSDTemp   = np.zeros(1)

        return ProcessesTemp, RatesTemp, RatesSDTemp
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Read_RatesAtT( self, InputData, Syst ):       

        if (Syst.NAtoms == 3):
            print('    [System.py - Read_RatesAtT]: Reading Rates for a System of 3 Atoms' )

            NStates0 = Syst.EqNStatesIn[0]
            print('    [System.py - Read_RatesAtT]: Nb of Levels Initial Molecules = ' + str(NStates0) )

            NProcTot = 1
            print('    [System.py - Read_RatesAtT]: Nb of Levels Per Final Pair: ')
            for iP in range(3):
                iMol                   = Syst.Pair[iP].ToMol
                Syst.Pair[iP].NStates  = Syst.EqNStatesIn[iMol]  
                NProcTot               = NProcTot + Syst.Pair[iP].NStates 
                Syst.Pair[iP].NProcTot = NProcTot
                print('    [System.py - Read_RatesAtT]:   Pair ' + str(iP) + ' = ' + str(Syst.Pair[iP].NStates) )


            self.Proc[0].Rates = np.zeros((NStates0, Syst.NProcTypes))
            self.Proc[1].Rates = np.zeros((NStates0, Syst.Pair[0].NStates))
            for iProc in range(2, 4):
                self.Proc[iProc].Rates       = np.zeros((NStates0, Syst.Pair[iP-1].NStates))
            for iExch in range(2, Syst.NProcTypes):
                self.ProcExch[iExch-2].Rates = np.zeros((NStates0, Syst.EqNStatesIn[Syst.ExchtoMol[iExch-2]]))

            ProcType = np.array([1, 2, 3]) 
            if (Syst.Pair[2].ToMol == Syst.Pair[1].ToMol ):
                ProcType[2] = ProcType[1]


            i  = 0
            ii = 10
            for iStates in range(Syst.Molecule[0].NStates):
                if (int(iStates/Syst.Molecule[0].NStates*100) == int(i*ii)):
                    print('    [System.py - Read_RatesAtT]: Read ' + str(i*ii) + '% of the Rate Files')
                    i = i+1 
                if ( (iStates >= InputData.Kin.MinStateIn[0] - 1) and (iStates <= InputData.Kin.MaxStateIn[0] - 1) ):

                    RatesTempAll                            = np.zeros(Syst.Pair[-1].NProcTot+1)
                    [ProcessesTemp, RatesTemp, RatesSDTemp] = self.Read_RatesFile( Syst, iStates )
                    RatesTempAll[ProcessesTemp[:]-1]        = RatesTemp[:]

                    RatesSplitted                       = np.split( RatesTempAll, np.array([1, Syst.Pair[0].NProcTot, Syst.Pair[1].NProcTot, Syst.Pair[2].NProcTot]) )
                    #self.Proc[0].Rates[iStates,0] = RatesSplitted[0]
                    self.Proc[0].Rates[iStates,0] = RatesSplitted[1][0:0]
                    self.Proc[1].Rates[iStates,:] = RatesSplitted[1][1:]
                    for iProc in range(2, 4):
                        self.Proc[0].Rates[iStates,ProcType[iProc-2]] = RatesSplitted[iProc][0:0]
                        self.Proc[0].Rates[iStates,0]                 = self.Proc[0].Rates[iStates,0] + RatesSplitted[iProc][0:0]
                        self.Proc[iProc].Rates[iStates,:]             = RatesSplitted[iProc][1:]
                
            for iProc in range(2, 4):   
                for iExch in range(2, Syst.NProcTypes):
                    if (Syst.Pair[iProc-1].ToMol == Syst.ExchtoMol[iExch-2] ):
                        self.ProcExch[iExch-2].Rates = self.ProcExch[iExch-2].Rates + self.Proc[iProc].Rates


        elif (Syst.NAtoms == 4):
            print('    [System.py - Read_RatesAtT]: Reading Rates for a System of 4 Atoms' )


            iPVec    = np.array([0, 1, 2])
            iPOppVec = np.array([5, 4, 3])

            NStates0_1 = Syst.EqNStatesIn[Syst.Pair[0].ToMol]
            NStates0_2 = Syst.EqNStatesIn[Syst.Pair[5].ToMol]
            print('    [System.py - Read_RatesAtT]: Nb of Levels in Initial Molecule 1 = ' + str(NStates0_1) )
            print('    [System.py - Read_RatesAtT]: Nb of Levels in Initial Molecule 2 = ' + str(NStates0_2) )

            maxNStates = 0
            print('    [System.py - Read_RatesAtT]: Nb of Levels Per Final Pair: ')
            for iP in range(6):
                iMol                   = Syst.Pair[iP].ToMol
                maxNStates             = np.maximum(Syst.EqNStatesIn[iMol], maxNStates)
                Syst.Pair[iP].NStates  = Syst.EqNStatesIn[iMol]  
                print('    [System.py - Read_RatesAtT]:   Pair ' + str(iP) + ' = ' + str(Syst.Pair[iP].NStates) )

            NProcTot = 1
            for iP in range(3):
                NProcTot               = NProcTot + (Syst.Pair[iP].NStates + 1) * (Syst.Pair[iPOppVec[iP]].NStates + 1) 
                Syst.Pair[iP].NProcTot = NProcTot


            self.DissRates     = np.zeros((NStates0_1, NStates0_2, Syst.NProcTypes))
            self.Proc[0].Rates = np.zeros((NStates0_1, NStates0_2, maxNStates, Syst.NDistMolecules+6))
            for iP in range(1, 4):
                print('    [System.py - Read_RatesAtT]: Pair ' + str(iP) + '; Rate Matrix shape = (' + str(NStates0_1) + '; ' + str(NStates0_2) + '; ' + str(Syst.Pair[iP-1].NStates) + '; ' + str(Syst.Pair[iPOppVec[iP-1]].NStates) + ')')
                self.Proc[iP].Rates          = np.zeros((NStates0_1, NStates0_2, Syst.Pair[iP-1].NStates, Syst.Pair[iPOppVec[iP-1]].NStates))
            for iExch in range(2, Syst.NProcTypes):
                self.ProcExch[iExch-2].Rates = np.zeros((NStates0_1, NStates0_2, Syst.EqNStatesIn[Syst.ExchtoMol[iExch-2,0]], Syst.EqNStatesIn[Syst.ExchtoMol[iExch-2,1]]))

            ProcType = np.array([1, 2, 3]) 
            if ( Syst.MolToCFDComp[Syst.Pair[2].ToMol] == Syst.MolToCFDComp[Syst.Pair[1].ToMol] ):
                ProcType[2] = ProcType[1]


            i  = 0
            ii = 10
            for iStates in range(NStates0_1):
                jStatesStart = 0 
                if (Syst.SymmFlg):
                    jStatesStart = iStates
                for jStates in range(NStates0_2):
                    i = i+1 
                    if (int((iStates*jStates)/(NStates0_1*NStates0_2)*100) == int(i*ii)):
                        print('    [System.py - Read_RatesAtT]: Read ' + str(i*ii) + '% of the Rate Files')

                    if ( ( (iStates >= InputData.Kin.MinStateIn[0] - 1) and (iStates <= InputData.Kin.MaxStateIn[0] - 1) ) and ( (jStates >= InputData.Kin.MinStateIn[1]) and (jStates <= InputData.Kin.MaxStateIn[1]) ) ):

                        if (jStates >= jStatesStart):
                            RatesTempAll                            = np.zeros(Syst.Pair[2].NProcTot)
                            [ProcessesTemp, RatesTemp, RatesSDTemp] = self.Read_RatesFile( Syst, i-1 )
                            RatesTempAll[ProcessesTemp[:]-1]        = RatesTemp[:]

                            iProc = -1

                            iProc = iProc + 1
                            self.DissRates[iStates, jStates, 0] = self.DissRates[iStates, jStates, 0] + RatesTempAll[iProc]

                            for iP in range(1, 4):
                                for kStates in range(Syst.Pair[iP-1].NStates+1):
                                    for lStates in range(Syst.Pair[iPOppVec[iP-1]].NStates+1):
                                        iProc = iProc + 1

                                        if ( (kStates == 0) and (lStates == 0) ): 
                                            self.DissRates[iStates, jStates, ProcType[iP-1]]           = self.DissRates[iStates, jStates, ProcType[iP-1]]           + RatesTempAll[iProc]
                                            self.DissRates[iStates, jStates,  0]                       = self.DissRates[iStates, jStates,  0]                       + RatesTempAll[iProc]
                                        elif (kStates == 0):
                                            iPair                                                      = iP-1
                                            iPairTemp                                                  = Syst.NDistMolecules + iPair 
                                            ToMol                                                      = Syst.CFDComp[ Syst.MolToCFDComp[ Syst.Pair[iPair].ToMol ] ].ToMol
                                            self.Proc[0].Rates[iStates, jStates, lStates-1, ToMol]     = self.Proc[0].Rates[iStates, jStates, lStates-1, ToMol]     + RatesTempAll[iProc]
                                            self.Proc[0].Rates[iStates, jStates, lStates-1, iPairTemp] = self.Proc[0].Rates[iStates, jStates, lStates-1, iPairTemp] + RatesTempAll[iProc]
                                        elif (lStates == 0):
                                            iPair                                                      = iPOppVec[iP-1]
                                            iPairTemp                                                  = Syst.NDistMolecules + iPair 
                                            ToMol                                                      = Syst.CFDComp[ Syst.MolToCFDComp[ Syst.Pair[iPair].ToMol ] ].ToMol
                                            self.Proc[0].Rates[iStates, jStates, kStates-1, ToMol]     = self.Proc[0].Rates[iStates, jStates, kStates-1, ToMol]     + RatesTempAll[iProc]
                                            self.Proc[0].Rates[iStates, jStates, kStates-1, iPairTemp] = self.Proc[0].Rates[iStates, jStates, kStates-1, iPairTemp] + RatesTempAll[iProc]
                                        else:
                                            iTemp1  = kStates-1
                                            iTemp2  = lStates-1
                                            SymmFlg = 1.0
                                            if (Syst.SymmFlg):
                                                if (iTemp1 > iTemp2):
                                                    iTemp3 = iTemp2
                                                    iTemp2 = iTemp1
                                                    iTemp1 = iTemp3
                                                elif (iTemp1 == iTemp2):
                                                    SymmFlg = 2.0
                                            self.Proc[iP].Rates[iStates, jStates, iTemp1, iTemp2]      = self.Proc[iP].Rates[iStates, jStates, iTemp1, iTemp2]      + RatesTempAll[iProc] * SymmFlg
                        
            for iProc in range(2, 4):
                for iExch in range(2, Syst.NProcTypes):
                    if ( ( Syst.MolToCFDComp[Syst.Pair[iProc-1].ToMol] == Syst.MolToCFDComp[Syst.ExchtoMol[iExch-2,0]] ) and ( Syst.MolToCFDComp[Syst.Pair[iPOppVec[iProc-1]].ToMol] == Syst.MolToCFDComp[Syst.ExchtoMol[iExch-2,1]] ) ):
                        print('    [System.py - Read_RatesAtT]: iPair ' + str(iProc) + ' corresponds to iExch ' + str(iExch-1) ) 
                        self.ProcExch[iExch-2].Rates = self.ProcExch[iExch-2].Rates + self.Proc[iProc].Rates
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Save_RatesAtT_HDF5( self, Syst ):

        if (Syst.NAtoms == 3):

            PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
            HDF5Exist_Flg = path.exists(PathToFile)
            #if (HDF5Exist_Flg):
            f = h5py.File(PathToFile, 'a')
            #else:
            #    f = {'key': 'value'}

            TStr = 'T_' + str(int(self.TTra)) + '_' + str(int(self.TInt)) + '/Rates/'       
            if TStr in f.keys():
                grp               = f[TStr]

                Data      = grp["Diss"]
                Data[...] = self.Proc[0].Rates
                Data      = grp["Inel"]
                Data[...] = self.Proc[1].Rates

                for iProc in range(2, Syst.NProcTypes):
                    ExchStr   = "Exch_" + str(iProc-1)
                    Data      = grp[ExchStr]
                    Data[...] = self.ProcExch[iProc-2].Rates

            else:
                grp           = f.create_group(TStr)

                Proc0         = grp.create_dataset("Diss", data=self.Proc[0].Rates, compression="gzip", compression_opts=9)
                Proc1         = grp.create_dataset("Inel", data=self.Proc[1].Rates, compression="gzip", compression_opts=9)

                for iProc in range(2, Syst.NProcTypes):
                    ExchStr   = "Exch_" + str(iProc-1)
                    ProcExchi = grp.create_dataset(ExchStr, data=self.ProcExch[iProc-2].Rates, compression="gzip", compression_opts=9)

            f.close()

        else:
            
            PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
            f          = h5py.File(PathToFile, 'a')

            TStr = 'T_' + str(int(self.TTra)) + '_' + str(int(self.TInt)) + '/Rates/'       
            if TStr in f.keys():
                grp               = f[TStr]

                Data      = grp["Diss"]
                Data[...] = self.DissRates
                Data      = grp["DissInel"]
                Data[...] = self.Proc[0].Rates
                Data      = grp["Inel"]
                Data[...] = self.Proc[1].Rates

                for iProc in range(2, Syst.NProcTypes):
                    ExchStr   = "Exch_" + str(iProc-1)
                    Data      = grp[ExchStr]
                    Data[...] = self.ProcExch[iProc-2].Rates

            else:
                grp           = f.create_group(TStr)

                Proc00        = grp.create_dataset("Diss",     data=self.DissRates,     compression="gzip", compression_opts=9)
                Proc0         = grp.create_dataset("DissInel", data=self.Proc[0].Rates, compression="gzip", compression_opts=9)
                Proc1         = grp.create_dataset("Inel",     data=self.Proc[1].Rates, compression="gzip", compression_opts=9)

                for iProc in range(2, Syst.NProcTypes):
                    ExchStr   = "Exch_" + str(iProc-1)
                    ProcExchi = grp.create_dataset(ExchStr,    data=self.ProcExch[iProc-2].Rates, compression="gzip", compression_opts=9)

            f.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Transform_ProcToDiss( self, Syst ):
        if (Syst.NAtoms == 3):
            print('    [System.py - Transform_ProcToDiss]: Transforming Processes in Dissociation for 3 Atoms System')

            for iLevel in range(Syst.Molecule[0].NLevels):
                if (Syst.Molecule[0].LevelWrite_Flg[iLevel]):
                    
                    for jLevel in range(Syst.Molecule[0].NLevels):
                        if (not Syst.Molecule[0].LevelWrite_Flg[jLevel]):
                            TempRate                      = self.Proc[1].Rates[iLevel, jLevel]
                            self.Proc[0].Rates[iLevel, 0] = self.Proc[0].Rates[iLevel, 0] + float(TempRate)
                            self.Proc[0].Rates[iLevel, 1] = self.Proc[0].Rates[iLevel, 1] + float(TempRate)

                    for iProc in range(2, Syst.NProcTypes):
                        jMol = Syst.ExchtoMol[iProc-2]
                        for jLevel in range(Syst.Molecule[jMol].NLevels):
                            if (not Syst.Molecule[0].LevelWrite_Flg[jLevel]):
                                TempRate                          = self.ProcExch[iProc-2].Rates[iLevel, jLevel]
                                self.Proc[0].Rates[iLevel, 0]     = self.Proc[0].Rates[iLevel, 0]     + float(TempRate)
                                self.Proc[0].Rates[iLevel, iProc] = self.Proc[0].Rates[iLevel, iProc] + float(TempRate)

        else:
            print('    [System.py - Transform_ProcToDiss]: ERROR! Transforming Processes in Dissociation for 4 Atoms System NOT IMPLEMENTED yet!')
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Compute_BackwardRates( self, Syst ):

        if (Syst.NAtoms == 3):
            print('    [System.py - Compute_BackwardRates]: Computing Backweard Rates for 3 Atoms System')

            self.Proc[1].BckRates = self.Proc[1].Rates
            for iLevel in range(Syst.Molecule[0].NLevels):
                for jLevel in range(Syst.Molecule[0].NLevels):
                    if (Syst.Molecule[0].LevelEeV[iLevel] < Syst.Molecule[0].LevelEeV[jLevel]):
                        self.Proc[1].BckRates[iLevel, jLevel] = self.Proc[1].Rates[jLevel, iLevel] * Syst.Molecule[0].T[self.iT-1].LevelQ[jLevel] / Syst.Molecule[0].T[self.iT-1].LevelQ[iLevel]


            for iProc in range(2, Syst.NProcTypes):
                jMol = Syst.ExchtoMol[iProc-2]

                self.ProcExch[iProc-2].BckRates = self.ProcExch[iProc-2].Rates
                for iLevel in range(Syst.Molecule[0].NLevels):
                    for jLevel in range(Syst.Molecule[jMol].NLevels):
                        if (Syst.Molecule[0].LevelEeV[iLevel] < Syst.Molecule[jMol].LevelEeV[jLevel]):
                            self.ProcExch[iProc-2].BckRates[iLevel, jLevel] = self.ProcExch[iProc-2].Rates[jLevel, iLevel] * Syst.Molecule[0].T[self.iT-1].LevelQ[jLevel] / Syst.Molecule[jMol].T[self.iT-1].LevelQ[iLevel]
        
        else:
            print('    [System.py - Compute_BackwardRates]: ERROR! Computation of Backweard Rates for 4 Atoms System NOT IMPLEMENTED yet!')
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Compute_GroupRates(self, Syst ):

        if (Syst.NAtoms == 3):
            print('    [System.py - Compute_GroupRates]: Computing Grouped Rates for 3 Atoms System')
            
            print('      [System.py - Compute_GroupRates]:   Computing Grouped Dissociation Rates for Temperature Nb ' + str(self.iT+1) )
            self.GroupedProc[0].Rates = np.zeros((Syst.EqNStatesOut[0],4))
            
            for iLevel in range(Syst.NLevels[0]):
                iGroup = Syst.Molecule[0].Mapping[iLevel]
                for iP in range(4):
                    self.GroupedProc[0].Rates[iGroup,iP] = self.GroupedProc[0].Rates[iGroup,iP] + self.Proc[0].Rates[iLevel,iP] * Syst.Molecule[0].GroupsOut.T[self.iT+1].Expvec[iLevel] / Syst.Molecule[0].GroupsOut.T[self.iT+1].Q[iGroup]

            print('      [System.py - Compute_GroupRates]:   Computing Grouped Inelastic    Rates for Temperature Nb ' + str(self.iT+1) )
            self.GroupedProc[1].Rates = np.zeros((Syst.EqNStatesOut[0],Syst.EqNStatesOut[0]))
            
            for iLevel in range(Syst.NLevels[0]):
                iGroup = self.Mapping[iLevel]
                for jLevel in range(Syst.NLevels[0]):
                    jGroup = self.Mapping[jLevel]
                    self.GroupedProc[1].Rates[iGroup,jGroup] = self.GroupedProc[1].Rates[iGroup,jGroup] + self.Proc[1].BckRates[iLevel,jLevel] * Syst.Molecule[0].GroupsOut.T[self.iT+1].Expvec[iLevel] / Syst.Molecule[0].GroupsOut.T[self.iT+1].Q[iGroup]

            print('      [System.py - Compute_GroupRates]:   Computing Grouped Exchange     Rates for Temperature Nb ' + str(self.iT+1) )
            for iProc in range(2, Syst.NProcTypes):
                self.ProcExch[iProc-2].Rates = np.zeros((Syst.EqNStatesOut[0],Syst.EqNStatesOut[iProc-1]))

                for iLevel in range(Syst.NLevels[0]):
                    iGroup = self.Mapping[iLevel]
                    for jLevel in range(Syst.NLevels[iProc-1]):
                        jGroup = self.Mapping[jLevel]
                        self.GroupedProcExch[iProc-2].Rates[iGroup,jGroup] = self.GroupedProcExch[iProc-2].Rates[iGroup,jGroup] + self.ProcExch[iProc-2].BckRates[iLevel,jLevel] * Syst.Molecule[0].GroupsOut.T[self.iT+1].Expvec[iLevel] / Syst.Molecule[0].GroupsOut.T[self.iT+1].Q[iGroup]
        
        if (Syst.NAtoms == 4):
            print('    [System.py - Compute_GroupRates]: ERROR! Computing Grouped Rates for 4 Atoms System NOT IMPLEMENTED Yet!')
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Correcting_DissRates( self, InputData ):
        print('    [System.py - Correcting_DissRates]: Correcting Dissociation Rates by Factor ' + str(InputData.Kin.CorrFactor) )

        self.Proc[0].Rates = self.Proc[0].Rates * InputData.Kin.CorrFactor
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Correcting_GroupedDissRates( self, InputData ):
        print('    [System.py - Correcting_GroupedDissRates]: Correcting Grouped Dissociation Rates by Factor ' + str(InputData.Kin.CorrFactor) )

        self.GroupedProc[0].Rates =  self.GroupedProc[0].Rates * InputData.Kin.CorrFactor
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Compute_Rates_Overall( self, InputData, Syst ):

        if (Syst.NAtoms == 3):
            print('    [System.py - Compute_Rates_Overall]: Computing Overall Rates for 3 Atoms System')

            self.ProcTot[0].Rates = self.Proc[0].Rates[:,0]
            self.ProcTot[1].Rates = np.sum(self.Proc[1].Rates, axis=1)
            for jProc in range(2, Syst.NProcTypes):
                self.ProcTot[jProc].Rates = np.sum(self.ProcExch[jProc-2].Rates, axis=1)
        
        elif (Syst.NAtoms == 4):
            print('    [System.py - Compute_Rates_Overall]: Computing Overall Rates for 4 Atoms System')

            iPVec    = np.array([0, 1, 2])
            iPOppVec = np.array([5, 4, 3])

            NStates0_1 = Syst.EqNStatesIn[Syst.Pair[0].ToMol]
            NStates0_2 = Syst.EqNStatesIn[Syst.Pair[5].ToMol]
            print('    [System.py - Compute_Rates_Overall]:   Nb of States ( = Levels / Groups) in Initial Molecule 1 = ' + str(NStates0_1) )
            print('    [System.py - Compute_Rates_Overall]:   Nb of States                      in Initial Molecule 2 = ' + str(NStates0_2) )

            maxNStates = 0
            for iP in range(6):
                iMol                   = Syst.Pair[iP].ToMol
                maxNStates             = np.maximum(Syst.EqNStatesIn[iMol], maxNStates)
                Syst.Pair[iP].NLevels  = Syst.EqNStatesIn[iMol]
                print('    [System.py - Compute_Rates_Overall]:     Nb of States in Pair ' + str(iP) + ' = ' + str(Syst.Pair[iP].NLevels) )
            print('    [System.py - Compute_Rates_Overall]:    Max Nb of States in Pair = ' + str(maxNStates) )

            self.DissRatesTot = self.DissRates[:,:,0]

            if (Syst.SymmFlg):
                iTemp = 1
            else:
                iTemp = 2
            
            self.ProcTot[0].Rates = np.zeros((NStates0_1, NStates0_2, iTemp))
            for iMol in range(iTemp):
                for iP in range(6):
                    if ( Syst.MolToCFDComp[Syst.Pair[iP].ToMol] == Syst.MolToCFDComp[iMol] ):
                        self.ProcTot[0].Rates[:,:,iMol] = self.ProcTot[0].Rates[:,:,iMol] + np.squeeze( np.sum(self.Proc[0].Rates[:,:,:,iP], axis=2) )
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Compute_PrefJumps( self, InputData, Syst ):

        NJumps = InputData.Kin.RatesNPrefJumps

        self.Proc[1].PrefJumps = np.zeros((Syst.Molecule[0].NLevels, NJumps), dtype=np.int32)
        for iLevel in range(Syst.Molecule[0].NLevels):
            TempVec = self.Proc[1].BckRates[iLevel, :]
            self.Proc[1].PrefJumps[iLevel,:] = np.argsort(TempVec)[-NJumps:] + 1

        for iProc in range(2, Syst.NProcTypes):
            jMol = Syst.ExchtoMol[iProc-2]
            
            self.ProcExch[iProc-2].PrefJumps = np.zeros((Syst.Molecule[0].NLevels, NJumps), dtype=np.int32)
            for iLevel in range(Syst.Molecule[0].NLevels):
                TempVec =  self.ProcExch[iProc-2].BckRates[iLevel, :]
                self.ProcExch[iProc-2].PrefJumps[iLevel,:] = np.argsort(TempVec)[-NJumps:] + 1
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Write_PrefJumps( self, InputData, Syst ):

        TempFldr   = Syst.PathToReadFldr + '/' + Syst.Molecule[0].Name + '/Rates/T_' + str(int(self.TTra)) + '_' + str(int(self.TInt))

        PathToFile = TempFldr + '/PrefJumps_Inel.csv'
        print('    [System.py - Write_PrefJumps]: Writing Jumps in File: ' + PathToFile ) 
        with open(PathToFile, 'w') as csvJumps:
            Line    = '# Top 5, Less Probable -> More Probable\n' 
            csvJumps.write(Line)
            TempMat = self.Proc[1].PrefJumps
            np.savetxt(csvJumps, TempMat.astype(int), delimiter=',', fmt='%d')
        csvJumps.close()

        for iProc in range(2, Syst.NProcTypes):
            PathToFile = TempFldr + '/PrefJumps_Exch_Type' + str(iProc-1) + '.csv'
            print('    [System.py - Write_PrefJumps]: Writing Jumps in File: ' + PathToFile ) 
            with open(PathToFile, 'w') as csvJumps:
                Line    = '# Top 5, Less Probable -> More Probable\n' 
                csvJumps.write(Line)
                TempMat = self.ProcExch[iProc-2].PrefJumps
                np.savetxt(csvJumps, TempMat.astype(int), delimiter=',', fmt='%d')
            csvJumps.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Write_Kinetics( self, InputData, Syst, Temp ):

        if (InputData.Kin.WriteExoth_Flg):
            TempCoeff =  1.0
        else:
            TempCoeff = -1.0


        if (InputData.Kin.WriteFormat == 'PLATO'):
            print('    [System.py - Write_Kinetics]: Writing Kinetics in PLATO Format' )
        elif (InputData.Kin.WriteFormat == 'csv'):
            print('    [System.py - Write_Kinetics]: Writing Kinetics in .csv Format' )
        elif (InputData.Kin.WriteFormat == 'custom'):
            print('    [System.py - Write_Kinetics]: Writing Kinetics in Custom Format' )
        else:
            raise NameError("    [System.py - Write_Kinetics]: InputData.Kin.WriteFormat must be specified. Select 'PLATO', 'csv', or 'custom'")


        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' ) 
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.GroupsOutSuffix )
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.GroupsOutSuffix + '/T' + str(int(Temp.TranVec[self.iT-1])) + 'K/' )    
        TempFldr = InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.GroupsOutSuffix + '/T' + str(int(Temp.TranVec[self.iT-1])) + 'K/'

        if (Syst.NAtoms == 3):
            print('    [System.py - Write_Kinetics]: Writing Kinetics for 3 Atoms System' )

            if (InputData.Kin.WriteDiss_Flg):
                if (InputData.Kin.CorrFactor != 1.0):
                    if (InputData.Kin.PackUnpackDiss_Flg):
                        DissKinetics = TempFldr + '/Diss' + InputData.Kin.PackUnpackSuffix + '.dat' 
                    else:
                        DissKinetics = TempFldr + '/Diss_Corrected.dat'
                    print('    [System.py - Write_Kinetics]: Writing Corrected Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
                else:
                    if (InputData.Kin.PackUnpackDiss_Flg):
                        DissKinetics = TempFldr + '/Diss' + InputData.Kin.PackUnpackSuffix + '.dat' 
                    else:
                        DissKinetics = TempFldr + '/Diss.dat' 
                    print('    [System.py - Write_Kinetics]: Writing Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
                csvkinetics  = open(DissKinetics, 'w')
                
                if (InputData.Kin.WriteFormat == 'csv'):
                    csvkinetics.write('#i,$k_i^{D}~[cm^3/s]$\n')
                elif (InputData.Kin.WriteFormat == 'custom'):
                    print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )
                    csvkinetics.write('# HEADER')


                for iLevel in range(Syst.Molecule[0].NLevels):
                    TempRate = self.Proc[0].Rates[iLevel,0]
                    if ( (TempRate > 0.0) and ( (iLevel >= InputData.Kin.MinStateOut[0] - 1) and (iLevel <= InputData.Kin.MaxStateOut[0] - 1) ) and (Syst.Molecule[0].LevelWrite_Flg[iLevel]) ):
                        iiLevel  = Syst.Molecule[0].LevelNewMapping[iLevel]

                        if (InputData.Kin.WriteFormat == 'PLATO'):
                            ProcName = Syst.Molecule[0].Name + '(' + str(iiLevel+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name
                            Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,2\n' % float(TempRate)
                        
                        elif (InputData.Kin.WriteFormat == 'csv'):
                            Line     = '%d,%e\n' % (iiLevel+1, float(TempRate))

                        elif (InputData.Kin.WriteFormat == 'custom'):
                            print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )

                        csvkinetics.write(Line)
                    
                csvkinetics.close()


            if (InputData.Kin.WriteInel_Flg):

                if (InputData.Kin.WindAvrg_Flg):
                    InelKinetics = TempFldr + '/Inel_WindAvrg.dat' 
                    print('    [System.py - Write_Kinetics]: Writing Window-Averaged Inelastic: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name )
                else:
                    InelKinetics = TempFldr + '/Inel.dat' 
                    print('    [System.py - Write_Kinetics]: Writing Inelastic: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name )
                csvkinetics  = open(InelKinetics, 'w')

                if (InputData.Kin.WriteFormat == 'csv'):
                    csvkinetics.write('#i,j,$k_{i,j}^{Inel}~[cm^3/s]$\n')
                elif (InputData.Kin.WriteFormat == 'custom'):
                    print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )
                    csvkinetics.write('# HEADER')

                for iLevel in range(Syst.Molecule[0].NLevels):
                    if ( ( (iLevel >= InputData.Kin.MinStateOut[0] - 1) and (iLevel <= InputData.Kin.MaxStateOut[0] - 1) ) and (Syst.Molecule[0].LevelWrite_Flg[iLevel]) ):
                        TempRates = self.Proc[1].Rates[iLevel,:]
                        if (InputData.Kin.WindAvrg_Flg):
                            TempRates = Syst.Compute_WindAvrg_Rates( TempRates )                
                        iiLevel  = Syst.Molecule[0].LevelNewMapping[iLevel]

                        for jLevel in range(Syst.Molecule[0].NLevels):
                            if ( ( (jLevel >= InputData.Kin.MinStateOut[1] - 1) and (jLevel <= InputData.Kin.MaxStateOut[1] - 1) ) and (Syst.Molecule[0].LevelWrite_Flg[jLevel]) ):
                                if ( (TempRates[jLevel] > 0.0) and ( TempCoeff * Syst.Molecule[0].LevelEEh[iLevel] > TempCoeff * Syst.Molecule[0].LevelEEh[jLevel]) ):
                                    jjLevel  = Syst.Molecule[0].LevelNewMapping[jLevel]

                                    if (InputData.Kin.WriteFormat == 'PLATO'):
                                        ProcName = Syst.Molecule[0].Name + '(' + str(iiLevel+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '(' + str(jjLevel+1) + ')+' + Syst.Atom[2].Name
                                        Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,5\n' % TempRates[jLevel]
                                    
                                    elif (InputData.Kin.WriteFormat == 'csv'):
                                        Line     = '%d,%d,%e\n' % (iiLevel+1, jjLevel+1, float(TempRate))

                                    elif (InputData.Kin.WriteFormat == 'custom'):
                                        print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )

                                    csvkinetics.write(Line)
                    
                csvkinetics.close()


            if (InputData.Kin.WriteExch_Flg):

                for iExch in range (2, Syst.NProcTypes):
                    
                    jToMol  = Syst.ExchtoMol[iExch-2]
                    jToAtom = Syst.ExchtoAtom[iExch-2] 
                    print('    [System.py - Write_Kinetics]: iExch =  ' + str(iExch-1) + ', Corresponding to Final Molecule Nb ' + str(jToMol) + ' and Atom Nb ' + str(jToAtom) )

                    if (InputData.Kin.WindAvrg_Flg):
                        ExchKinetics = TempFldr + 'Exch_Type' + str(iExch-1) + '_WindAvrg.dat' 
                        print('    [System.py - Write_Kinetics]: Writing Window-Averaged Exchange Nb. '+ str(iExch-1) + ': ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[jToMol].Name  + '+' + Syst.Atom[jToAtom].Name  )
                    else:
                        ExchKinetics = TempFldr + 'Exch_Type' + str(iExch-1) + '.dat' 
                        print('    [System.py - Write_Kinetics]: Writing Exchange Nb. '+ str(iExch-1) + ': ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[jToMol].Name  + '+' + Syst.Atom[jToAtom].Name  )
                    csvkinetics  = open(ExchKinetics, 'w')

                    if (InputData.Kin.WriteFormat == 'csv'):
                        csvkinetics.write('#i,j,$k_{i,j}^{Exch}~[cm^3/s]$\n')
                    elif (InputData.Kin.WriteFormat == 'custom'):
                        print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )
                        csvkinetics.write('# HEADER')

                    for iLevel in range(Syst.Molecule[0].NLevels):
                        if ( ( (iLevel >= InputData.Kin.MinStateOut[0] - 1) and (iLevel <= InputData.Kin.MaxStateOut[0] - 1) ) and (Syst.Molecule[0].LevelWrite_Flg[iLevel]) ):
                            TempRates = self.ProcExch[iExch-2].Rates[iLevel,:]
                            if (InputData.Kin.WindAvrg_Flg):
                                TempRates = Syst.Compute_WindAvrg_Rates( TempRates )                    
                            iiLevel  = Syst.Molecule[0].LevelNewMapping[iLevel]

                            for jLevel in range(Syst.Molecule[jToMol].NLevels):
                                if ( ( (jLevel >= InputData.Kin.MinStateOut[1] - 1) and (jLevel <= InputData.Kin.MaxStateOut[1] - 1) ) and (Syst.Molecule[jToMol].LevelWrite_Flg[jLevel]) ):
                                    if ((TempRates[jLevel] > 0.0) and ( TempCoeff * Syst.Molecule[0].LevelEEh[iLevel] > TempCoeff * Syst.Molecule[jToMol].LevelEEh[jLevel]) ):
                                        jjLevel  = Syst.Molecule[jToMol].LevelNewMapping[jLevel]
                                        
                                        if (InputData.Kin.WriteFormat == 'PLATO'):
                                            ProcName = Syst.Molecule[0].Name + '(' + str(iiLevel+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[jToMol].Name + '(' + str(jjLevel+1) + ')+' + Syst.Atom[jToAtom].Name
                                            Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,6\n' % TempRates[jLevel]
                                        
                                        elif (InputData.Kin.WriteFormat == 'csv'):
                                            Line     = '%d,%d,%e\n' % (iiLevel+1, jjLevel+1, float(TempRate))

                                        elif (InputData.Kin.WriteFormat == 'custom'):
                                            print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )

                                        csvkinetics.write(Line)

                    csvkinetics.close()


        elif (Syst.NAtoms == 4):
            print('    [System.py - Write_Kinetics]: Writing Kinetics for 4 Atoms System' )

            iPVec    = np.array([0, 1, 2])
            iPOppVec = np.array([5, 4, 3])

            NStates0_1 = Syst.EqNStatesIn[Syst.Pair[0].ToMol]
            NStates0_2 = Syst.EqNStatesIn[Syst.Pair[5].ToMol]
            print('    [System.py - Write_Kinetics]: Nb of States in Initial Molecule 1 = ' + str(NStates0_1) )
            print('    [System.py - Write_Kinetics]: Nb of States in Initial Molecule 2 = ' + str(NStates0_2) )


            if (InputData.Kin.WriteDiss_Flg):

                if (InputData.Kin.CorrFactor != 1.0):
                    DissKinetics = TempFldr + '/Diss_Corrected.dat' 
                    print('    [System.py - Write_Kinetics]: Writing Corrected Dissociation: ' + Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[0].ToMol]].Name + '+' + Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[5].ToMol]].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name + '+' + Syst.Atom[3].Name )
                else:
                    DissKinetics = TempFldr + '/Diss.dat' 
                    print('    [System.py - Write_Kinetics]: Writing Dissociation: '           + Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[0].ToMol]].Name + '+' + Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[5].ToMol]].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name + '+' + Syst.Atom[3].Name )
                csvkinetics  = open(DissKinetics, 'w')

                if (InputData.Kin.WriteFormat == 'csv'):
                    csvkinetics.write('#i,j,$k_{i,j}^{D}~[cm^3/s]$\n')
                elif (InputData.Kin.WriteFormat == 'custom'):
                    print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )
                    csvkinetics.write('# HEADER')

                for iStates in range(NStates0_1):
                    iiStates = Syst.Molecule[Syst.Pair[0].ToMol].LevelNewMapping[iStates]
                    jStatesStart = 0 
                    if (Syst.SymmFlg):
                        jStatesStart = iStates
                    for jStates in range(jStatesStart, NStates0_2):
                        jjStates = Syst.Molecule[Syst.Pair[5].ToMol].LevelNewMapping[jStates]
                        if ( ( (iStates >= InputData.Kin.MinStateOut[0] - 1) and (iStates <= InputData.Kin.MaxStateOut[0] - 1) ) and (Syst.Molecule[Syst.Pair[0].ToMol].LevelWrite_Flg[iStates]) and ( (jStates >= InputData.Kin.MinStateOut[1] - 1) and (jStates <= InputData.Kin.MaxStateOut[1] - 1) and (Syst.Molecule[Syst.Pair[5].ToMol].LevelWrite_Flg[jStates]) ) ):

                            TempRate = self.DissRates[iStates, jStates, 0]
                            if (TempRate > 0.0):

                                if (InputData.Kin.WriteFormat == 'PLATO'):
                                    Mol1_Str = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[0].ToMol]].Name + '(' + str(iiStates+1) + ')'
                                    Mol2_Str = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[5].ToMol]].Name + '(' + str(jjStates+1) + ')'
                                    LHS_Str  = Mol1_Str + '+' + Mol2_Str
                                    RHS_Str  = Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name + '+' + Syst.Atom[3].Name
                                    ProcName = LHS_Str  + '=' + RHS_Str
                                    Line     = ProcName + ':+%.4e,+0.0000E+00,+0.0000E+00,2\n' % float(TempRate)
                                
                                elif (InputData.Kin.WriteFormat == 'csv'):
                                    Line     = '%d,%d,%e\n' % (iiStates+1, jjStates+1, float(TempRate))

                                elif (InputData.Kin.WriteFormat == 'custom'):
                                    print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )

                                csvkinetics.write(Line)
                    
                csvkinetics.close()


            if (InputData.Kin.WriteDissInel_Flg):

                Mol1_Str     = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[0].ToMol]].Name 
                Mol2_Str     = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[5].ToMol]].Name
                LHS_Str      = Mol1_Str + '+' + Mol2_Str
                print('    [System.py - Write_Kinetics]: Writing Inelastic + Dissociation: ' + LHS_Str + '= Molecule + Atom + Atom' )
              
                InelKinetics = TempFldr + '/DissInel.dat' 
                csvkinetics  = open(InelKinetics, 'w')

                if (InputData.Kin.WriteFormat == 'csv'):
                    csvkinetics.write('#i,j,k,$k_{i,j,k}^{D\-Inel}~[cm^3/s]$\n')
                elif (InputData.Kin.WriteFormat == 'custom'):
                    print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )
                    csvkinetics.write('# HEADER')

                for iStates in range(NStates0_1):
                    iiStates = Syst.Molecule[Syst.Pair[0].ToMol].LevelNewMapping[iStates]
                    jStatesStart = 0 
                    SymmFct    = 1.0
                    if (Syst.SymmFlg):
                        jStatesStart = iStates
                        #SymmFct    = 2.0
                    for jStates in range(jStatesStart, NStates0_2):
                        jjStates = Syst.Molecule[Syst.Pair[5].ToMol].LevelNewMapping[jStates]
                        if ( ( (iStates >= InputData.Kin.MinStateOut[0] - 1) and (iStates <= InputData.Kin.MaxStateOut[0] - 1) ) and (Syst.Molecule[Syst.Pair[0].ToMol].LevelWrite_Flg[iStates]) and ( (jStates >= InputData.Kin.MinStateOut[1] - 1) and (jStates <= InputData.Kin.MaxStateOut[1] - 1) and (Syst.Molecule[Syst.Pair[5].ToMol].LevelWrite_Flg[jStates]) ) ):

                            for iComp in range(Syst.NCFDComp):
                                if (Syst.CFDComp[iComp].ToMol >= 0):
                                    iMol    = Syst.CFDComp[iComp].ToMol
                                    NBins_1 = Syst.EqNStatesIn[iMol]
                                    for kStates in range(NBins_1):
                                        kkStates = Syst.Molecule[iMol].LevelNewMapping[kStates]
                                        if ( ( (kStates >= InputData.Kin.MinStateOut[2] - 1) and (kStates <= InputData.Kin.MaxStateOut[2] - 1) ) and (Syst.Molecule[iMol].LevelWrite_Flg[kStates]) ):

                                            TempRate = self.Proc[0].Rates[iStates, jStates, kStates, iMol]
                                            if (TempRate > 0.0):

                                                if (InputData.Kin.WriteFormat == 'PLATO'):
                                                    Mol1_Str = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[0].ToMol]].Name + '(' + str(iiStates+1) + ')'
                                                    Mol2_Str = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[5].ToMol]].Name + '(' + str(jjStates+1) + ')'
                                                    LHS_Str  = Mol1_Str + '+' + Mol2_Str
                                                    Mol1_Str = Syst.CFDComp[iComp].Name                                 + '(' + str(kkStates+1) + ')'
                                                    Atms_Str = Syst.CFDComp[Syst.CFDComp[iComp].ToOppAtoms[0]].Name     + '+' + Syst.CFDComp[Syst.CFDComp[iComp].ToOppAtoms[1]].Name
                                                    RHS_Str  = Mol1_Str + '+' + Atms_Str
                                                    ProcName = LHS_Str  + '=' + RHS_Str
                                                    Line     = ProcName + ':+%.4e,+0.0000E+00,+0.0000E+00,2\n' % float(TempRate)
                                                
                                                elif (InputData.Kin.WriteFormat == 'csv'):
                                                    Line     = '%d,%d,%d,%e\n' % (iiStates+1, jjStates+1, kkStates+1, float(TempRate))

                                                elif (InputData.Kin.WriteFormat == 'custom'):
                                                    print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )

                                                csvkinetics.write(Line)

                csvkinetics.close()


            if (InputData.Kin.WriteInel_Flg):
              
                Mol1_Str     = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[0].ToMol]].Name 
                Mol2_Str     = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[5].ToMol]].Name
                LHS_Str      = Mol1_Str + '+' + Mol2_Str
                print('    [System.py - Write_Kinetics]: Writing Inelastic: ' + LHS_Str + '=' + LHS_Str )
                
                InelKinetics = TempFldr + '/Inel.dat' 
                csvkinetics  = open(InelKinetics, 'w')

                if (InputData.Kin.WriteFormat == 'csv'):
                    csvkinetics.write('#i,j,k,l,$k_{i,j,k,l}^{Inel}~[cm^3/s]$\n')
                elif (InputData.Kin.WriteFormat == 'custom'):
                    print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )
                    csvkinetics.write('# HEADER')

                for iStates in range(NStates0_1):
                    iiStates     = Syst.Molecule[Syst.Pair[0].ToMol].LevelNewMapping[iStates]
                    jStatesStart = 0 
                    if (Syst.SymmFlg):
                        jStatesStart = iStates
                    for jStates in range(jStatesStart, NStates0_2):
                        jjStates = Syst.Molecule[Syst.Pair[5].ToMol].LevelNewMapping[jStates]
                        if ( ( (iStates >= InputData.Kin.MinStateOut[0] - 1) and (iStates <= InputData.Kin.MaxStateOut[0] - 1) ) and (Syst.Molecule[Syst.Pair[0].ToMol].LevelWrite_Flg[iStates]) and ( (jStates >= InputData.Kin.MinStateOut[1] - 1) and (jStates <= InputData.Kin.MaxStateOut[1] - 1) and (Syst.Molecule[Syst.Pair[5].ToMol].LevelWrite_Flg[jStates]) ) ):

                            for kStates in range(NStates0_1):
                                kkStates = Syst.Molecule[Syst.Pair[0].ToMol].LevelNewMapping[kStates]
                                lStatesStart = 0
                                SymmFct    = 1.0
                                if (Syst.SymmFlg):
                                    lStatesStart = kStates
                                    #SymmFct    = 2.0
                                for lStates in range(lStatesStart, NStates0_2):
                                    llStates = Syst.Molecule[Syst.Pair[5].ToMol].LevelNewMapping[lStates]
                                    if ( ( (kStates >= InputData.Kin.MinStateOut[2] - 1) and (kStates <= InputData.Kin.MaxStateOut[2] - 1) ) and (Syst.Molecule[Syst.Pair[0].ToMol].LevelWrite_Flg[kStates]) and ( (lStates >= InputData.Kin.MinStateOut[3] - 1) and (lStates <= InputData.Kin.MaxStateOut[3] - 1) and (Syst.Molecule[Syst.Pair[5].ToMol].LevelWrite_Flg[lStates]) ) ):

                                        TempRate = self.Proc[1].Rates[iStates, jStates, kStates, lStates]           * SymmFct
                                        InEEh    = Syst.Molecule[Syst.Pair[0].ToMol].T[self.iT-1].EqEeV0In[iStates] + Syst.Molecule[Syst.Pair[5].ToMol].T[self.iT-1].EqEeV0In[jStates]
                                        FinEEh   = Syst.Molecule[Syst.Pair[0].ToMol].T[self.iT-1].EqEeV0In[kStates] + Syst.Molecule[Syst.Pair[5].ToMol].T[self.iT-1].EqEeV0In[lStates]
                                        
                                        if ((TempRate > 0.0) and ( TempCoeff * InEEh >= TempCoeff * FinEEh )):

                                            if (InputData.Kin.WriteFormat == 'PLATO'):
                                                Mol1_Str     = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[0].ToMol]].Name + '(' + str(iiStates+1) + ')'
                                                Mol2_Str     = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[5].ToMol]].Name + '(' + str(jjStates+1) + ')'
                                                LHS_Str      = Mol1_Str + '+' + Mol2_Str
                                                Mol1_Str = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[0].ToMol]].Name + '(' + str(kkStates+1) + ')'
                                                Mol2_Str = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[5].ToMol]].Name + '(' + str(llStates+1) + ')'
                                                RHS_Str  = Mol1_Str + '+' + Mol2_Str
                                                ProcName = LHS_Str  + '=' + RHS_Str
                                                Line     = ProcName + ':+%.4e,+0.0000E+00,+0.0000E+00,5\n' % float(TempRate)

                                            elif (InputData.Kin.WriteFormat == 'csv'):
                                                Line     = '%d,%d,%d,%d,%e\n' % (iiStates+1, jjStates+1, kkStates+1, llStates+1, float(TempRate))

                                            elif (InputData.Kin.WriteFormat == 'custom'):
                                                print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )

                                            csvkinetics.write(Line)

                csvkinetics.close()


            if (InputData.Kin.WriteExch_Flg):

                for iExch in range (2, Syst.NProcTypes):

                    Mol1_Str     = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[0].ToMol]].Name 
                    Mol2_Str     = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[5].ToMol]].Name
                    LHS_Str      = Mol1_Str + '+' + Mol2_Str
                    Mol1_Str     = Syst.CFDComp[Syst.MolToCFDComp[Syst.ExchtoMol[iExch-2,0]]].Name
                    Mol2_Str     = Syst.CFDComp[Syst.MolToCFDComp[Syst.ExchtoMol[iExch-2,1]]].Name
                    RHS_Str      = Mol1_Str + '+' + Mol2_Str
                    print('    [System.py - Write_Kinetics]: Writing Exchange Nb. '+ str(iExch-1) + ': ' + LHS_Str + '=' + RHS_Str )

                    ExchKinetics = TempFldr + 'Exch_Type' + str(iExch-1) + '.dat' 
                    csvkinetics  = open(ExchKinetics, 'w')

                    if (InputData.Kin.WriteFormat == 'csv'):
                        csvkinetics.write('#i,j,k,l,$k_{i,j,k,l}^{Exch}~[cm^3/s]$\n')
                    elif (InputData.Kin.WriteFormat == 'custom'):
                        print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )
                        csvkinetics.write('# HEADER')

                    NBins_1 = Syst.EqNStatesIn[Syst.ExchtoMol[iExch-2,0]]
                    NBins_2 = Syst.EqNStatesIn[Syst.ExchtoMol[iExch-2,1]]

                    for iStates in range(NStates0_1):
                        iiStates = Syst.Molecule[Syst.Pair[0].ToMol].LevelNewMapping[iStates]
                        jStatesStart = 0 
                        if (Syst.SymmFlg):
                            jStatesStart = iStates
                        for jStates in range(jStatesStart, NStates0_2):
                            jjStates = Syst.Molecule[Syst.Pair[5].ToMol].LevelNewMapping[jStates]
                            if ( ( (iStates >= InputData.Kin.MinStateOut[0] - 1) and (iStates <= InputData.Kin.MaxStateOut[0] - 1) ) and (Syst.Molecule[Syst.Pair[0].ToMol].LevelWrite_Flg[iStates]) and ( (jStates >= InputData.Kin.MinStateOut[1] - 1) and (jStates <= InputData.Kin.MaxStateOut[1] - 1) and (Syst.Molecule[Syst.Pair[5].ToMol].LevelWrite_Flg[jStates]) ) ):

                                for kStates in range(NBins_1):
                                    kMol     = Syst.ExchtoMol[iExch-2,0]
                                    kkStates = Syst.Molecule[kMol].LevelNewMapping[kStates]
                                    lStatesStart = 0
                                    SymmFct      = 1.0
                                    if (Syst.SymmFlg):
                                        lStatesStart = kStates
                                        #SymmFct    = 2.0
                                    for lStates in range(lStatesStart, NBins_2):
                                        lMol     = Syst.ExchtoMol[iExch-2,1]
                                        llStates = Syst.Molecule[lMol].LevelNewMapping[lStates]
                                        if ( ( (kStates >= InputData.Kin.MinStateOut[2] - 1) and (kStates <= InputData.Kin.MaxStateOut[2] - 1) ) and (Syst.Molecule[kMol].LevelWrite_Flg[kStates]) and ( (lStates >= InputData.Kin.MinStateOut[3] - 1) and (lStates <= InputData.Kin.MaxStateOut[3] - 1) and (Syst.Molecule[lMol].LevelWrite_Flg[lStates]) ) ):
 
                                            TempRate = self.ProcExch[iExch-2].Rates[iStates, jStates, kStates, lStates]        * SymmFct
                                            InEEh    = Syst.Molecule[Syst.Pair[0].ToMol].T[self.iT-1].EqEeV0In[iStates]        + Syst.Molecule[Syst.Pair[5].ToMol].T[self.iT-1].EqEeV0In[jStates]
                                            FinEEh   = Syst.Molecule[Syst.ExchtoMol[iExch-2,0]].T[self.iT-1].EqEeV0In[kStates] + Syst.Molecule[Syst.ExchtoMol[iExch-2,1]].T[self.iT-1].EqEeV0In[lStates]
                                            
                                            if ((TempRate > 0.0) and ( TempCoeff * InEEh >= TempCoeff * FinEEh )):

                                                if (InputData.Kin.WriteFormat == 'PLATO'):
                                                    Mol1_Str = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[0].ToMol]].Name + '(' + str(iiStates+1) + ')'
                                                    Mol2_Str = Syst.CFDComp[Syst.MolToCFDComp[Syst.Pair[5].ToMol]].Name + '(' + str(jjStates+1) + ')'
                                                    LHS_Str  = Mol1_Str + '+' + Mol2_Str
                                                    Mol1_Str = Syst.CFDComp[Syst.MolToCFDComp[kMol]].Name + '(' + str(kkStates+1) + ')'
                                                    Mol2_Str = Syst.CFDComp[Syst.MolToCFDComp[lMol]].Name + '(' + str(llStates+1) + ')'
                                                    RHS_Str  = Mol1_Str + '+' + Mol2_Str
                                                    ProcName = LHS_Str  + '=' + RHS_Str
                                                    Line     = ProcName + ':+%.4e,+0.0000E+00,+0.0000E+00,6\n' % float(TempRate)
                                                
                                                elif (InputData.Kin.WriteFormat == 'csv'):
                                                    Line     = '%d,%d,%d,%d,%e\n' % (iiStates+1, jjStates+1, kkStates+1, llStates+1, float(TempRate))

                                                elif (InputData.Kin.WriteFormat == 'custom'):
                                                    print('    [System.py - Write_Kinetics]: Customize the File by Specifying the Desired Format' )

                                                csvkinetics.write(Line)

                    csvkinetics.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Write_GroupedKinetics( self, InputData, Syst, Temp ):

        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' ) 
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.GroupsOutSuffix )
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.GroupsOutSuffix + '/T' + str(int(Temp.TranVec[self.iT-1])) + 'K/' )    
        TempFldr = InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.GroupsOutSuffix + '/T' + str(int(Temp.TranVec[self.iT-1])) + 'K/'

        if (InputData.Kin.WriteDiss_Flg):
            if (InputData.Kin.CorrFactor != 1.0):
                DissKinetics = TempFldr + '/Diss_Corrected.dat' 
                print('    [System.py - Write_GroupedKinetics]: Writing Corrected Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
            else:
                DissKinetics = TempFldr + '/Diss.dat' 
                print('    [System.py - Write_GroupedKinetics]: Writing Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
            csvkinetics  = open(DissKinetics, 'w')

            for iGroup in range(self.EqNStatesOut[0]):

                if (self.GroupedProc[0].Rates[iGroup,0] > 0.0):
                    ProcName = Syst.Molecule[0].Name + '(' + str(iGroup+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name
                    Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,2\n' % float(self.GroupedProc[0].Rates[iGroup,0])
                    csvkinetics.write(Line)
            
            csvkinetics.close()


        if (InputData.Kin.WriteInel_Flg):
  
            InelKinetics = TempFldr + '/Inel.dat' 
            print('    [System.py - Write_GroupedKinetics]: Writing Inelastic: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name )
            csvkinetics  = open(InelKinetics, 'w')
            #InelFile     = InputData.Kin.ReadFldr  + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '.csv'


            for iGroup in range(self.EqNStatesOut[0]):
                TempRates = self.GroupedProc[1].Rates[iGroup,:]
                if (InputData.Kin.WindAvrg_Flg):
                    TempRates = Syst.Compute_WindAvrg_Rates( TempRates )                

                for jGroup in range(self.EqNStatesOut[0]):
                    if ((TempRates[jGroup] > 0.0) and (Syst.Molecule[0].GroupsOut.T[iT].EeV[iGroup] > Syst.Molecule[0].GroupsOut.T[iT].EeV[jGroup]) ):
                        ProcName = Syst.Molecule[0].Name + '(' + str(iGroup+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '(' + str(jGroup+1) + ')+' + Syst.Atom[2].Name
                        Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,5\n' % TempRates[jGroup]
                        csvkinetics.write(Line)
                        
            csvkinetics.close()


        if (InputData.Kin.WriteExch_Flg):

            for iExch in range (2, Syst.NProcTypes):
                print('    [System.py - Write_GroupedKinetics]: iExch =  ' + str(iExch-1) )

                ExchKinetics = TempFldr + 'Exch_Type' + str(iExch-1) + '.dat' 
                print('    [System.py - Write_GroupedKinetics]: Writing Exchange Nb. '+ str(iExch-1) + ': ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name  + '+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name  )
                csvkinetics  = open(ExchKinetics, 'w')
                #InelFile     = InputData.Kin.ReadFldr                               + '/'  + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name  + '+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name + '.csv'

                for iGroup in range(self.EqNStatesOut[0]):
                    TempRates = self.GroupedProcExch[iExch-2].Rates[iGroup,:]                  

                    for jGroup in range(self.EqNStatesOut[iExch-1]):
                        if ((TempRates[jGroup] > 0.0) and (Syst.Molecule[0].GroupsOut.T[iT].EeV[iGroup] > Syst.Molecule[Syst.ExchtoMol[iExch-2]].GroupsOut.T[iT].EeV[jGroup]) ):
                            ProcName = Syst.Molecule[0].Name + '(' + str(iGroup+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name + '(' + str(jGroup+1) + ')+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name
                            Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,6\n' % TempRates[jGroup]
                            csvkinetics.write(Line)

                csvkinetics.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def PackUnpackDiss( self, InputData, Syst ):

        if (Syst.NAtoms == 3):
            print('    [System.py - PackUnpackDiss]: Packing and Unpacking Dissociation StS Rates for 3 Atoms System for Temperature Nb ' + str(self.iT-1) + ' (T = ' + str(int(self.TTra)) + 'K)')

            NGroups = Syst.Molecule[0].PackUnpack.NGroups
            NLevels = Syst.Molecule[0].NLevels
            print('    [System.py - PackUnpackDiss]:   Nb of Groups = ' + str(NGroups) + '; Nb of Levels = ' + str(NLevels) )

            GroupsDissRates = np.zeros((NGroups)) 
            for iLevel in range(NLevels):
                if (Syst.Molecule[0].LevelWrite_Flg[iLevel]):
                    iGroup                  = Syst.Molecule[0].PackUnpack.Mapping[iLevel]
                    GroupsDissRates[iGroup] = GroupsDissRates[iGroup] + self.Proc[0].Rates[iLevel, 0] * Syst.Molecule[0].T[self.iT-1].LevelQ[iLevel] / Syst.Molecule[0].PackUnpack.T[self.iT].Q[iGroup]
            
            self.Proc[0].Rates = np.zeros((NLevels,4)) 
            for iLevel in range(NLevels):
                if (Syst.Molecule[0].LevelWrite_Flg[iLevel]):
                    iGroup                        = Syst.Molecule[0].PackUnpack.Mapping[iLevel]
                    self.Proc[0].Rates[iLevel, 0] = float(GroupsDissRates[iGroup])

        if (Syst.NAtoms == 4):
            print('    [System.py - PackUnpackDiss]: ERROR! Packing and Unpacking Dissociation StS Rates for 4 Atoms System NOT IMPLEMENTED Yet!')
    # ...........................................................................................................................



# --------------------------------------------------------------------------------------------------------------------- CLASS ---





# ===================================================================================================================== CLASS ===
class arrhenius(object):

    def __init__(self, NProcTypes):

        self.Proc    = [processes() for iProc in range(NProcTypes)]
        self.MinRate = 0.0
        self.MinNbTs = 3

# --------------------------------------------------------------------------------------------------------------------- CLASS ---





# ===================================================================================================================== CLASS ===
class system(object):

    def __init__(self, SystNameLong, SystName, NAtoms, NMolecules, NDistMolecules, NPairs, NCFDComp, NTTran, NProcTypes):

        self.NameLong     = SystNameLong
        self.Name         = SystName

        self.NAtoms       = NAtoms
        self.Atom         = [atom() for iA in range(NAtoms)]

        self.NMolecules     = NMolecules
        self.NDistMolecules = NDistMolecules
        self.Molecule       = [molecule(NTTran) for iMol in range(NMolecules)]

        self.NPairs       = NPairs
        self.Pair         = [pair() for iP in range(NPairs)]

        self.NCFDComp     = NCFDComp
        self.CFDComp      = [cfdcomp() for iComp in range(NCFDComp)]
        self.MolToCFDComp = 0

        self.PathToFolder = ''
        self.PathToHDF5   = ''

        self.NProcTypes   = NProcTypes

        self.T            = [t_properties(NProcTypes) for iT in range(NTTran)]

        self.Arr          = arrhenius(NProcTypes)

        self.Proc         = [processes() for iProc in range(NProcTypes)]

        self.RatesTh      = np.zeros((NTTran, NProcTypes))
        #self.RatesQSS     = np.zeros((NTTran, NProcTypes))

        self.MolToCFDComp = np.zeros((NMolecules), dtype=np.int8)



    # ***************************************************************************************************************************
    def Read_Rates( self, InputData, Temp ):
        print('  [System.py - Read_Rates]: Uploading Rates')


        if ( (InputData.Kin.WindAvrg_Flg) and (InputData.Kin.Write_Flg) ):
            print('  [System.py - Read_Rates]: Initializing Objects for Window-Averaging the Kinetics')
            
            if (self.NAtoms == 3):
                self.Molecule[0].Initialize_WindAvrg_Objects( InputData, self )
            else:
                print('  [System.py - Read_Rates]:   ERROR! Window-Averaging not implemented for Nb Atoms > 3!')


        if ( (InputData.OldVersion_IntFlg == 2) and (not self.HDF5Exist_Flg) ):
            print("  [System.py - Read_Rates]: Asked To Read N3 Rates in An Old Format (Panesi's RVC 2013)" )
            self.Read_RatesFile_N3( InputData )


        for iT in Temp.iTVec:
            self.T[iT-1].iT   = iT
            self.T[iT-1].TTra = Temp.TranVec[iT-1]
            self.T[iT-1].TInt = self.T[iT-1].TTra

            PathToFile  = self.PathToHDF5 + '/' + self.NameLong + '.hdf5'
            f           = h5py.File(PathToFile, 'a')
            TStr        = 'T_' + str(int(self.T[iT-1].TTra)) + '_' + str(int(self.T[iT-1].TInt)) + '/Rates/'
            TPresent_Flg = TStr in f.keys()
            f.close()
            if (TPresent_Flg): print('  [System.py - Read_Rates]: Found Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K) in the HDF5 File')

            if ( (TPresent_Flg) and (not InputData.HDF5.ForceReadDat_Flg) ):
                print('  [System.py - Read_Rates]: Loading Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K) from HDF5 File')
            
                self.T[iT-1].Load_RatesAtT_HDF5( self )

            else:
                print('  [System.py - Read_Rates]: Reading Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K) from set of .dat Files')

                if (InputData.OldVersion_IntFlg == 1):
                    self.T[iT-1].Read_RatesAtT_OldVersion( InputData, self )
                elif (InputData.OldVersion_IntFlg == 0):
                    self.T[iT-1].Read_RatesAtT( InputData, self )
                elif (InputData.OldVersion_IntFlg == 2):
                    self.T[iT-1].Read_RatesAtT_N3( InputData, self )

                if (InputData.HDF5.Save_Flg):
                    print('  [System.py - Read_Rates]: Saving Rates Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K) in the HDF5 File')
                    
                    self.T[iT-1].Save_RatesAtT_HDF5( self )
            

            if (InputData.Kin.WriteQB_IntFlg < 2):
                print('  [System.py - Read_Rates]: Considering Dissociation the Inelastic and Exchange Processes to Excluded Levels')
                self.T[iT-1].Transform_ProcToDiss( self )


            if ( (InputData.Kin.RatesPrefJumps_Flg) or (InputData.Kin.GroupsOut_Flg)):
                print('  [System.py - Read_Rates]: Computing Backweard Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K)')
                self.T[iT-1].Compute_BackwardRates( self )


            if (InputData.Kin.GroupsOut_Flg):
                print('  [System.py - Read_Rates]: Computing Grouped Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K)')
                self.T[iT-1].Compute_GroupRates( self )


            if (InputData.Kin.CorrFactor != 1.0):   
                    print('  [System.py - Read_Rates]: Correcting Dissociation Rates (Corr Factor = ' + str(InputData.Kin.CorrFactor) + ') for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K)')
                    
                    self.T[iT-1].Correcting_DissRates( InputData )
                    
                    if (InputData.Kin.GroupsOut_Flg):
                        print('  [System.py - Read_Rates]: Correcting Grouped Dissociation Rates (Corr Factor = ' + str(InputData.Kin.CorrFactor) + ') for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K)')
                    
                        self.T[iT-1].Correcting_GroupedDissRates( InputData )


            print('  [System.py - Read_Rates]: Computing Overall Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K)')
            self.T[iT-1].Compute_Rates_Overall( InputData, self )


            # print('  [System.py - Read_Rates]: Computing Thermal Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K)')
            # self.Compute_ThermalRates( iT )


            if (InputData.Kin.RatesPrefJumps_Flg):
                
                print('  [System.py - Read_Rates]: Computing Preferred Jumps for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K)')
                self.T[iT-1].Compute_PrefJumps( InputData, self )
                
                print('  [System.py - Read_Rates]: Writing   Preferred Jumps for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K)')
                self.T[iT-1].Write_PrefJumps( InputData, self )
                

            if (InputData.Kin.PackUnpackDiss_Flg):
                print('  [System.py - Read_Rates]: Packing and Unpacking Dissociation Rates for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K)')
                self.T[iT-1].PackUnpackDiss( InputData, self )

            if (InputData.Kin.Write_Flg):
                print('  [System.py - Read_Rates]: Writing Kinetics File for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K)')
                self.T[iT-1].Write_Kinetics( InputData, self, Temp )


            if  (InputData.Kin.GroupsOut_Flg) and (InputData.Kin.GroupsOutWrite_Flg):
                print('  [System.py - Read_Rates]: Writing Grouped Kinetics File for Temperature Nb ' + str(iT) + ' (T = ' + str(int(self.T[iT-1].TTra)) + 'K)')                
                self.T[iT-1].Write_GroupedKinetics( InputData, self, Temp )


            if (InputData.DelRateMat_Flg):
                if (self.NAtoms == 3):
                    for iProc in range(4):
                        del self.T[iT-1].Proc[iProc].Rates
                    for iProc in range(2, self.NProcTypes):
                        del self.T[iT-1].ProcExch[iProc-2].Rates
                else:
                    for iProc in range(4):
                        del self.T[iT-1].Proc[iProc].Rates
                    for iProc in range(2, self.NProcTypes):
                        del self.T[iT-1].ProcExch[iProc-2].Rates

        # print('  [System.py - Read_Rates]: Saving Thermal Rates\n')
        # self.Write_ThermalRates( InputData, Temp )

        # print('  [System.py - Read_Rates]: Plotting Thermal Rates\n')
        # self.Plot_ThermalRates( InputData, Temp )
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Compute_WindAvrg_Rates( self, TempTempRates ):
        
        TempRates = np.zeros( (self.Molecule[0].NLevels) )
        for jLevel in range(self.Molecule[0].NLevels):
            TempRates[jLevel] = np.sum(TempTempRates[self.Molecule[0].WindAvrgMat[jLevel,0:self.Molecule[0].WindAvrgFound[jLevel,0]+1]]) / (self.Molecule[0].WindAvrgFound[jLevel,0]+1)

        return TempRates 
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Read_RatesFile_N3( self, InputData ):
    #sed -i 's/D-/E-/g' *
        print("    [System.py - Read_RatesFile_N3]: Reading Rates for N3 System in Old Format (Panesi's RVC 2013)" )

        self.N3NLevels   = 9390

        
        PathToFile      = InputData.PathToN3 + 'Mapping_OLD_TO_NEW.csv'
        Data            = pandas.read_csv(PathToFile, header=None, skiprows=1)
        Data            = Data.apply(pandas.to_numeric, errors='coerce')
        self.N3Mapping  = np.array(Data[1].values, dtype=np.int64)
        print('    [System.py - Read_RatesFile_N3]: Old to New Levels List Mapping: ',  self.N3Mapping )


        self.N3Data_TVec = np.array( [1000.0, 2500.0, 3500.0, 7500.0, 10000.0, 12500.0, 15000.0, 20000.0, 25000.0, 30000.0, 40000.0, 50000.0] )


        self.N3Data_Diss = np.zeros( (self.N3NLevels, np.size(self.N3Data_TVec)) )
        PathToFile       = InputData.PathToN3 + '/n3.ratecoef3_dissoc'
        with open(PathToFile) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=" ", skipinitialspace=True)
            for iLine in range(9):
                next(readCSV)
            iLine = -1
            jLine = 0
            for row in readCSV:
                iLine = iLine + 1
                if ( (iLine % 2 == 1) and (jLine < self.N3NLevels) ):
                    self.N3Data_Diss[jLine,:] = row[0::2]
                    jLine                     = jLine + 1
        print("    [System.py - Read_RatesFile_N3]:   Matrix of Dissociation Rates: ", self.N3Data_Diss )


        PathToFile       = InputData.PathToN3 + '/FullSet.ascii.dat'
        Data             = pandas.read_csv(PathToFile, header=None, delimiter=r"\s+")
        self.N3Data_Inel = Data.apply(pandas.to_numeric, errors='coerce')
        print("    [System.py - Read_RatesFile_N3]:   Matrix of Inelastic + Exchange Rates: ", self.N3Data_Inel )
    # ...........................................................................................................................



    # # ***************************************************************************************************************************
    # def Compute_ThermalRates( self, iT ):

    #     if (self.NAtoms == 3):
    #         for jProc in range(self.NProcTypes):
    #             self.RatesTh[iT-1,jProc] = sum( self.Molecule[0].T[iT-1].LevelQ * self.T[iT-1].ProcTot[jProc].Rates )
    #     else:
    #         print('    [System.py - Compute_Rates_Thermal]:   ERROR! Thermal Rates Computation not implemented for Nb Atoms > 3!')
    # # ...........................................................................................................................



    # # ***************************************************************************************************************************
    # def Compute_ThermalRates_FromOverall( self, InputData, Temp ):

    #     if (self.NAtoms == 3):
            
    #         #DissFile = InputData.Kin.ReadFldr   + '/' + self.Molecule[0].Name + '+' + self.Atom[2].Name + '_' + self.Atom[0].Name + '+' + self.Atom[1].Name + '+' + self.Atom[2].Name + '.csv'
    #         print('  [Compute_Rates_Thermal_FromOverall]: Reading Dissociation Rates From File: ' + DissFile)
    #         for iT in Temp.iTVec:
    #             LevelKDiss = np.zeros(self.Molecule[0].NLevels)
    #             with open(DissFile) as csvfile:
    #                 readCSV = csv.reader(csvfile, delimiter=',')
    #                 next(readCSV)
    #                 for row in readCSV:
    #                     if (float(row[iT]) > 0.0):
    #                         LevelKDiss[int(row[0])-1] = float(row[iT])
    #             csvfile.close()

    #             self.RatesTh[iT-1,0] = sum( self.Molecule[0].T[iT-1].LevelQ * LevelKDiss )
            

    #         print('    [System.py - Compute_Rates_Thermal_FromOverall]: Thermal Dissociation Rates = ', self.T[iT-1].ProcTot[0].Rates)
    #         self.Write_ThermalRates_Diss( InputData, Temp )

    #         print('    [System.py - Compute_Rates_Thermal_FromOverall]: Plotting Dissociation Thermal Rates')
    #         self.Plot_ThermalRates_Diss( InputData, Temp )
        
    #     else:
    #         print('    [System.py - Compute_Rates_Thermal_FromOverall]:   ERROR! Thermal Rates Computation From Overall not implemented for Nb Atoms > 3!')
    # # ...........................................................................................................................



    # # ***************************************************************************************************************************
    # def Write_ThermalRates_Diss( self, InputData, Temp ):

    #     mkdirs( InputData.OutputWriteFldr )    
    #     PathToFile = InputData.OutputWriteFldr + '/KTh_Diss.csv'
    #     print('      [System.py - Write_ThermalRates_Diss]: Writing Dissociation Thermal Rates in File: ' + PathToFile )
    #     if (not path.exists(PathToFile) ):
    #         WriteFlg = True
    #     else:
    #         WriteFlg = False        
    #     with open(PathToFile, 'a') as csvTermo:
    #         if (WriteFlg):
    #             Line    = '# T,KDiss\n' 
    #             csvTermo.write(Line)
    #         TempMat = np.concatenate( (np.expand_dims(np.array(Temp.TranVec, dtype='float'), axis=1), self.RatesTh[:,0]), axis=1 )
    #         np.savetxt(csvTermo, TempMat, delimiter=',')
    #     csvTermo.close()
    # # ...........................................................................................................................



    # # ***************************************************************************************************************************
    # def Plot_ThermalRates_Diss( self, InputData, Temp ):

    #     plt.figure()
    #     plt.title(r'$K_{i}^{D}$', fontsize=FontSize)

    #     LabelStr = '$K_{i}^{D}$ for ' + self.Molecule[0].Name + ', Interpolated'
    #     plt.plot(10000/Temp.TranVec, self.RatesTh[:,0], '-k', label=LabelStr)
        
    #     LabelStr = '$K_{i}^{D}$ for ' + self.Molecule[0].Name
    #     plt.plot(10000/Temp.TranVec, self.RatesTh[:,0], 'ok', label=LabelStr)
    #     plt.yscale("log")
        
    #     plt.xlabel(r'10,000/T [1/K]',         fontsize=FontSize)
    #     plt.ylabel(r'$K_{i}^{D}$ $[cm^3/s]$\n', fontsize=FontSize)
    #     plt.tight_layout()
    #     if (InputData.PlotShow_Flg):
    #         plt.show()
    #     FigSavePath = InputData.OutputWriteFldr + '/KTh_Diss.png'
    #     plt.savefig(FigSavePath)
    #     print('      [System.py - Plot_ThermalRates_Diss]: Saved Dissociation Thermal Rates Plot in: ' + FigSavePath)
    # # ...........................................................................................................................



    # # ***************************************************************************************************************************
    # def Write_Rates_Thermal( self, InputData, Temp ):

    #     mkdirs( InputData.OutputWriteFldr )    
    #     PathToFile = InputData.OutputWriteFldr + '/KTh.csv'
    #     print('    [System.py - Write_Rates_Thermal]: Writing Thermal Rates in File: ' + PathToFile )
    #     if (not path.exists(PathToFile) ):
    #         WriteFlg = True
    #     else:
    #         WriteFlg = False        
    #     with open(PathToFile, 'a') as csvTermo:
    #         if (WriteFlg):
    #             Line       = '# T,KDiss,KInel' 
    #             for iExch in range(2, self.NProcTypes):
    #                 Line = Line + ',KExch' + str(iExch-1)
    #             Line = Line + '\n'
    #             csvTermo.write(Line)
    #         TempMat = np.concatenate( (np.expand_dims(np.array(Temp.TranVec, dtype='float'), axis=1), self.RatesTh), axis=1 )
    #         np.savetxt(csvTermo, TempMat, delimiter=',')
    #     csvTermo.close()
    # # ...........................................................................................................................



    # # ***************************************************************************************************************************
    # def Write_Kinetics_FromOverall( self, InputData, Temp ):
        
    #     mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' ) 
    #     mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + self.NameLong + InputData.Kin.GroupsOutSuffix )

    #     for iT in Temp.iTVec:
    #         print('\nTemperature Nb ', iT, '; T = ', Temp.TranVec[iT-1], 'K')

    #         mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + self.NameLong + InputData.Kin.GroupsOutSuffix + '/T' + str(int(Temp.TranVec[iT-1])) + 'K/' )    
    #         TempFldr = InputData.Kin.WriteFldr + '/kinetics/' + self.NameLong + InputData.Kin.GroupsOutSuffix + '/T' + str(int(Temp.TranVec[iT-1])) + 'K/'


    #         if (InputData.Kin.WriteInel_Flg == True):
    #             InelKinetics = TempFldr + '/Inel.dat' 
    #             csvkinetics  = open(InelKinetics, 'w')

    #             print('  Writing Inelastic: ' + self.Molecule[0].Name + '+' + self.Atom[2].Name + '=' + self.Molecule[0].Name + '+' + self.Atom[2].Name )
    #             InelFile     = InputData.Kin.ReadFldr  + '/' + self.Molecule[0].Name + '+' + self.Atom[2].Name + '_' + self.Molecule[0].Name + '+' + self.Atom[2].Name + '.csv'
    #             with open(InelFile) as csvfile:
    #                 readCSV = csv.reader(csvfile, delimiter=',')
    #                 next(readCSV)
    #                 for row in readCSV:

    #                     if (float(row[iT+1]) > 0.0):
    #                         ProcName = self.Molecule[0].Name + '(' + str(row[0]) + ')+' + self.Atom[2].Name + '=' + self.Molecule[0].Name + '(' + str(row[1]) + ')+' + self.Atom[2].Name
    #                         Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,5\n' % float(row[iT+1])
    #                         csvkinetics.write(Line)
                    
    #                 csvfile.close()
    #                 csvkinetics.close()


    #         if (InputData.Kin.WriteExch_Flg == True):

    #             for iExch in range (2, self.NProcTypes):
    #                 print('  Writing Exchange: ' + self.Molecule[0].Name + '+' + self.Atom[2].Name + '=' + self.Molecule[self.ExchtoMol[iExch-2]].Name + '+' + self.Atom[self.ExchtoAtom[iExch-2]].Name )
    #                 ExchKinetics = TempFldr + '/Exch_Type' + str(iExch-1) + '.dat' 
    #                 csvkinetics  = open(ExchKinetics, 'w')

    #                 ExchFile     = InputData.Kin.ReadFldr + '/' + self.Molecule[0].Name + '+' + self.Atom[2].Name + '_' + self.Molecule[self.ExchtoMol[iExch-2]].Name + '+' + self.Atom[self.ExchtoAtom[iExch-2]].Name + '_Exch.csv'
    #                 with open(ExchFile) as csvfile:
    #                     readCSV = csv.reader(csvfile, delimiter=',')
    #                     next(readCSV)
    #                     for row in readCSV:

    #                         if (float(row[iT+1]) > 0.0):
    #                             ProcName = self.Molecule[0].Name + '(' + str(row[0]) + ')+' + self.Atom[2].Name + '=' + self.Molecule[self.ExchtoMol[iExch-2]].Name + '(' + str(row[1]) + ')+' + self.Atom[self.ExchtoAtom[iExch-2]].Name
    #                             Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,6\n' % float(row[iT+1])
    #                             csvkinetics.write(Line)
                    
    #                     csvfile.close()
    #                     csvkinetics.close()


    #         if (InputData.Kin.WriteDiss_Flg == True):
    #             if (InputData.Kin.CorrFactor != 1.0):
    #                 print('  Writing Corrected Dissociation: ' + self.Molecule[0].Name + '+' + self.Atom[2].Name + '=' + self.Atom[0].Name + '+' + self.Atom[1].Name + '+' + self.Atom[2].Name )
    #                 DissKinetics = TempFldr + '/Diss_Corrected.dat' 
    #             else:
    #                 print('  Writing Dissociation: ' + self.Molecule[0].Name + '+' + self.Atom[2].Name + '=' + self.Atom[0].Name + '+' + self.Atom[1].Name + '+' + self.Atom[2].Name )
    #                 DissKinetics = TempFldr + '/Diss.dat' 
    #             csvkinetics  = open(DissKinetics, 'w')

    #             DissFile = InputData.Kin.ReadFldr   + '/' + self.Molecule[0].Name + '+' + self.Atom[2].Name + '_' + self.Atom[0].Name + '+' + self.Atom[1].Name + '+' + self.Atom[2].Name + '.csv'
    #             with open(DissFile) as csvfile:
    #                 readCSV = csv.reader(csvfile, delimiter=',')
    #                 next(readCSV)
    #                 for row in readCSV:
    #                     #print(row[:])

    #                     if (float(row[iT]) > 0.0):
    #                         ProcName = self.Molecule[0].Name + '(' + str(row[0]) + ')+' + self.Atom[2].Name + '=' + self.Atom[0].Name + '+' + self.Atom[1].Name + '+' + self.Atom[2].Name
    #                         TempRate = float(row[iT])
    #                         Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,2\n' % (TempRate * InputData.Kin.CorrFactor)
    #                         csvkinetics.write(Line)
    #                     else:
    #                         ProcName = self.Molecule[0].Name + '(' + str(row[0]) + ')+' + self.Atom[2].Name + '=' + self.Atom[0].Name + '+' + self.Atom[1].Name + '+' + self.Atom[2].Name
    #                         TempRate = float(row[iT])
    #                         Line     = ProcName + ':+1.0000E-20,+0.0000E+00,+0.0000E+00,2\n'
    #                         csvkinetics.write(Line)
                    
    #                 csvfile.close()
    #                 csvkinetics.close()
    # # ...........................................................................................................................



# --------------------------------------------------------------------------------------------------------------------- CLASS ---