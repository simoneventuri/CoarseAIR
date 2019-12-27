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



    def Compute_ERot(self):

        self.vEeV              = np.zeros(np.max(self.Levelvqn)+1)
        self.LevelEVib         = np.zeros(self.NBins)
        self.LevelERot         = np.zeros(self.NBins)
        for iLevel in range(self.NBins):
            if (self.Leveljqn[iLevel] == 0):
                self.vEeV[int(self.Levelvqn[iLevel])] = self.LevelEeV[iLevel]
            self.LevelEVib[iLevel] = self.vEeV[int(self.Levelvqn[iLevel])]
            self.LevelERot[iLevel] = self.LevelEeV[iLevel] - self.LevelEVib[iLevel]

