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

from Processes  import processes
from Parameters import *


class grouped_t_properties(object):

    def __init__(self, NProcTypes, T):

        self.Value    = T

        self.GroupEeV = 0.0
        self.QRatio   = 0.0
        self.Q        = 0.0

        self.Proc     = [processes() for iProc in range(4)]
        self.ProcExch = [processes() for iProc in range(NProcTypes-2)]
        self.ProcTot  = [processes() for iProc in range(NProcTypes)]


class groupedmolecule(object):

    def __init__(self, InputData, NProcTypes, Temp):

        self.Type          = InputData.Kin.Group.Type
        self.PathToMapping = InputData.Kin.Group.PathToMapping
        self.NbTs          = Temp.NTran+1
        self.TVec          = np.array([InputData.Kin.Group.T0, Temp.TranVec])
        self.T             = [grouped_t_properties(NProcTypes, self.TVec[iT]) for iT in range(self.NbTs)]
s

    def Compute_GroupProps(self, Mol):

        for iT in range(self.NbTs):

            self.T[iT].Expvec   = Mol.Levelg * np.exp( - Mol.LevelEeV * Ue / (self.T[iT].Value * UKb) )
            self.T[iT].Q        = np.zeros((self.NbGroups))
            self.T[iT].GroupEeV = np.zeros((self.NbGroups))

            for iLevel in range(Mol.NBins):
                self.T[iT].Q[Mapping[iLevel]]        = self.T[iT].Q[Mapping[iLevel]]        +                        self.T[iT].Expvec[iLevel]
                self.T[iT].GroupEeV[Mapping[iLevel]] = self.T[iT].GroupEeV[Mapping[iLevel]] + Mol.LevelEeV[iLevel] * self.T[iT].Expvec[iLevel]


    def Compute_GroupRates(self, Syst, iT):
        

        self.T[iT+1].Proc[0].Rates = np.zeros((self.NbGroups,3))
        
        for iLevel in range(Syst.Molecule[0].NBins):
            for iP in range(3):
                self.T[iT+1].Proc[0].Rates[self.Mapping[iLevel],iP] = self.T[iT+1].Proc[0].Rates[self.Mapping[iLevel],iP] + Syst.T[iT].Proc[0].Rates[iLevel,iP] * self.T[iT+1].Expvec[iLevel]



        for iProc in range(2, 4):
            self.T[iT+1].Proc[iProc].Rates = np.zeros((self.NbGroups,self.NbGroups))

            for iLevel in range(Syst.Molecule[0].NBins):
                for jLevel in range(Syst.Molecule[0].NBins):
                    self.T[iT+1].Proc[iProc].Rates[self.Mapping[iLevel],self.Mapping[jLevel]] = self.T[iT+1].Proc[iProc].Rates[self.Mapping[iLevel],self.Mapping[jLevel]] + Syst.T[iT].Proc[iProc].Rates[iLevel,jLevel] * self.T[iT+1].Expvec[iLevel]


        
        for iProc in range(2, Syst.NProcTypes):
            self.T[iT+1].ProcExch[iProc-2].Rates = np.zeros((self.NbGroups,self.NbGroups))

            for iLevel in range(Syst.Molecule[0].NBins):
                for jLevel in range(Syst.Molecule[0].NBins):
                    self.T[iT+1].ProcExch[iProc-2].Rates[self.Mapping[iLevel],self.Mapping[jLevel]] = self.T[iT+1].ProcExch[iProc-2].Rates[self.Mapping[iLevel],self.Mapping[jLevel]] + Syst.T[iT].ProcExch[iProc-2].Rates[iLevel,jLevel] * self.T[iT+1].Expvec[iLevel]


        
        # for iProc in range(0, Syst.NProcTypes):
        #     self.T[iT+1].ProcTot[iProc].Rates = np.zeros((self.NbGroups))