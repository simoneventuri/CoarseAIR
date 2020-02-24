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
import shutil
import os, errno
import os.path
import pandas
from os import path
    
from Processes  import processes
from Parameters import *

def mkdirs(newdir, mode=0o777):
    os.makedirs(newdir, mode, exist_ok=True)


class grouped_t_properties(object):

    def __init__(self, NProcTypes, T):

        self.Value    = T

        self.EeV = 0.0
        self.EeV0       = 0.0
        self.QRatio   = 0.0
        self.Q        = 0.0
        self.QThermo  = 0.0

        self.Proc     = [processes() for iProc in range(4)]
        self.ProcExch = [processes() for iProc in range(NProcTypes-2)]
        self.ProcTot  = [processes() for iProc in range(NProcTypes)]


class groupedmolecule(object):

    def __init__(self, Type, PathToMapping, T0, NProcTypes, Temp, Name):

        self.Name          = Name
        self.Type          = Type
        self.PathToMapping = PathToMapping
        self.NbTs          = Temp.NTran+1
        self.TVec          = np.concatenate( (np.array([T0]), Temp.TranVec), axis=0 )
        print(self.NbTs)
        print(self.TVec)
        self.T             = [grouped_t_properties(NProcTypes, self.TVec[iT]) for iT in range(self.NbTs)]


    def Get_Mapping(self, Levelvqn):

        if (self.Type == 'VibSpecific'):
            self.NbGroups = np.max(Levelvqn)+1
            self.Mapping  = Levelvqn
        else:
            Data = pandas.read_csv(self.PathToMapping, header=None, skiprows=1)
            Data = Data.apply(pandas.to_numeric, errors='coerce')
            Idx  = np.array(Data.values[:,0], dtype=np.int64  ) - 1
            Bin  = np.array(Data.values[:,1], dtype=np.int64  ) - 1
            self.Mapping         = np.zeros((np.max(Idx)+1), dtype=np.int64) 
            self.Mapping[Idx[:]] = Bin[:]
            self.NbGroups        = np.max(Bin)+1


    def Compute_GroupProps(self, NLevels, Levelg, LevelEeV):

        Flg       = np.zeros((self.NbGroups), dtype=np.int8)
        self.EeV0 = np.zeros((self.NbGroups))
        for iLevel in range(NLevels):
            if (Flg[self.Mapping[iLevel]] == 0):
                Flg[self.Mapping[iLevel]]     = 1
                self.EeV0[self.Mapping[iLevel]] = LevelEeV[iLevel]

        for iT in range(self.NbTs):

            self.T[iT].Expvec       = Levelg * np.exp( -                 LevelEeV * Ue / (self.T[iT].Value * UKb) )
            self.T[iT].ExpvecThermo = Levelg * np.exp( - (LevelEeV - LevelEeV[0]) * Ue / (self.T[iT].Value * UKb) )
            self.T[iT].Q            = np.zeros((self.NbGroups))
            self.T[iT].QThermo      = np.zeros((self.NbGroups))
            self.T[iT].EeV          = np.zeros((self.NbGroups))
            self.T[iT].Ratio        = np.zeros((self.NbGroups))

            for iLevel in range(NLevels):
                self.T[iT].Q[self.Mapping[iLevel]]       = self.T[iT].Q[self.Mapping[iLevel]]        +                    self.T[iT].Expvec[iLevel]
                self.T[iT].QThermo[self.Mapping[iLevel]] = self.T[iT].QThermo[self.Mapping[iLevel]]  +                    self.T[iT].ExpvecThermo[iLevel]
                self.T[iT].EeV[self.Mapping[iLevel]]     = self.T[iT].EeV[self.Mapping[iLevel]] + LevelEeV[iLevel] * self.T[iT].Expvec[iLevel]
            self.T[iT].EeV   = self.T[iT].EeV / self.T[iT].Q
            self.T[iT].Ratio = self.T[iT].Q        / np.sum(self.T[iT].Q)


    def Compute_GroupRates(self, Syst, iT, NbLevels, NbGroups):
        
        print('    [Compute_GroupRates]: Computing Grouped Dissociation Rates for Temperature Nb ' + str(iT+1) )
        self.T[iT+1].Proc[0].Rates = np.zeros((NbGroups[0],3))
        
        for iLevel in range(NbLevels[0]):
            iGroup = self.Mapping[iLevel]
            for iP in range(3):
                self.T[iT+1].Proc[0].Rates[iGroup,iP] = self.T[iT+1].Proc[0].Rates[iGroup,iP] + Syst.T[iT].Proc[0].Rates[iLevel,iP] * self.T[iT+1].Expvec[iLevel] / self.T[iT+1].Q[iGroup]


        print('    [Compute_GroupRates]: Computing Grouped Inelastic    Rates for Temperature Nb ' + str(iT+1) )
        self.T[iT+1].Proc[1].Rates = np.zeros((NbGroups[0],NbGroups[0]))
        
        for iLevel in range(NbLevels[0]):
            iGroup = self.Mapping[iLevel]
            for jLevel in range(NbLevels[0]):
                jGroup = self.Mapping[jLevel]
                self.T[iT+1].Proc[1].Rates[iGroup,jGroup] = self.T[iT+1].Proc[1].Rates[iGroup,jGroup] + Syst.T[iT].Proc[1].BckRates[iLevel,jLevel] * self.T[iT+1].Expvec[iLevel] / self.T[iT+1].Q[iGroup]


        print('    [Compute_GroupRates]: Computing Grouped Exchange     Rates for Temperature Nb ' + str(iT+1) )
        for iProc in range(2, Syst.NProcTypes):
            self.T[iT+1].ProcExch[iProc-2].Rates = np.zeros((NbGroups[0],NbGroups[iProc-1]))

            for iLevel in range(NbLevels[0]):
                iGroup = self.Mapping[iLevel]
                for jLevel in range(NbLevels[iProc-1]):
                    jGroup = self.Mapping[jLevel]
                    self.T[iT+1].ProcExch[iProc-2].Rates[iGroup,jGroup] = self.T[iT+1].ProcExch[iProc-2].Rates[iGroup,jGroup] + Syst.T[iT].ProcExch[iProc-2].BckRates[iLevel,jLevel] * self.T[iT+1].Expvec[iLevel] / self.T[iT+1].Q[iGroup]

        # for iProc in range(0, Syst.NProcTypes):
        #     self.T[iT+1].ProcTot[iProc].Rates = np.zeros((self.NbGroups))
        
        print(' ')


    def Write_Groups_PartFuncsAndEnergies(self, Syst, InputData, Temp):

        mkdirs(    InputData.Kin.WriteFldr + '/thermo/' ) 
        mkdirs(    InputData.Kin.WriteFldr + '/thermo/' + Syst.NameLong + InputData.Kin.Groups.FldrName )
        TempFldr = InputData.Kin.WriteFldr + '/thermo/' + Syst.NameLong + InputData.Kin.Groups.FldrName

        MoleFracsFile = TempFldr + '/' + self.Name + '_InitialMoleFracs_T' + str(int(self.T[0].Value)) + 'K.dat' 
        print('      [Write_Groups_PartFuncsAndEnergies]: Writing Initial Mole Fractions for Molecule: ' + self.Name )
        csvmole       = open(MoleFracsFile, 'w')
        Line          = '# Percentage of ' + self.Name + ' Contained in Each Group\n'
        csvmole.write(Line)

        for iGroup in range(self.NbGroups):
            Line     = '%.10e\n' % float(self.T[0].Ratio[iGroup])
            csvmole.write(Line)

        csvmole.close()


        for iT in range(1, self.NbTs):
            print('      [Write_Groups_PartFuncsAndEnergies]: Writing Thermo File for Molecule: ' + self.Name + ' at T = ' + str(int(self.TVec[iT])) + ' K' )
            PathToFileOrig = TempFldr + '/../' + self.Name + '_Format'
            PathToFile     = TempFldr + '/'    + self.Name + '_T' + str(int(self.TVec[iT])) + 'K'
            DestTemp       = shutil.copyfile(PathToFileOrig, PathToFile)

            with open(PathToFile, 'a') as f:
                Line = 'NB_ENERGY_LEVELS = ' + str(self.NbGroups) + '\n'
                f.write(Line)
                np.savetxt(f, np.transpose(np.array([self.T[iT].QThermo, self.T[iT].QThermo*0.0])), fmt='%.8e    %.8e')
            f.close()


class t_properties(object):

    def __init__(self):

        self.LevelEeV = 0.0
        self.QRatio   = 0.0
        self.Q        = 0.0
        self.QExp     = 0.0
      

class molecule(object):

    def __init__(self, NTTran):

        self.Name              = ""
        self.DissEn            = np.array([1, 1, 1])
        self.DegeneracyFactor  = 0
        self.Mu                = 0

        self.NBins             = 0
        self.KinMthd           = ""

        self.StartBin          = 0
        self.FinalBin          = 0

        self.Levelvqn          = 0
        self.Leveljqn          = 0
        self.LevelEEh          = 0.0
        self.LevelEeV          = 0.0
        self.LevelERot         = 0.0
        self.LevelEVib         = 0.0
        self.LevelEGam         = 0.0
        self.LevelrMin         = 0.0
        self.LevelrMax         = 0.0
        self.LevelVMin         = 0.0
        self.LevelVMax         = 0.0
        self.LevelTau          = 0.0
        self.LevelrIn          = 0.0
        self.LevelrOut         = 0.0

        self.Levelg            = 0.0
        self.LevelQ            = 0.0
        self.LevelToBin        = 0

        self.NTTran            = NTTran
        self.T                 = [t_properties() for iTTran in range(NTTran)]


    def Comput_QRatio0(self, T0):

        self.QRatio0 = self.Levelg * np.exp( - self.LevelEeV * Ue / (T0 * UKb) )
        self.QRatio0 = self.QRatio0  / np.sum(self.QRatio0)


    def Compute_ERot(self):

        self.vEeV              = np.zeros(np.max(self.Levelvqn)+1)
        self.LevelEVib         = np.zeros(self.NBins)
        self.LevelERot         = np.zeros(self.NBins)
        for iLevel in range(self.NBins):
            if (self.Leveljqn[iLevel] == 0):
                self.vEeV[int(self.Levelvqn[iLevel])] = self.LevelEeV[iLevel]
            self.LevelEVib[iLevel] = self.vEeV[int(self.Levelvqn[iLevel])]
            self.LevelERot[iLevel] = self.LevelEeV[iLevel] - self.LevelEVib[iLevel]