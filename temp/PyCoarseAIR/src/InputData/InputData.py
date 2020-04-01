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


class rates(object):

    def __init__(self):

        self.PrefJumps_Flg    = False
        self.NPrefJumps       = 5



class kinetics(object):

    def __init__( self, WORKSPACE_PATH ):

        ## Resolution of the Kinetics Data in Input? Array of 'StS' / 'VSM' / 'CGM' of size Syst.NMolecules
        self.MolResolutionIn            = ['VSM', 'VSM', 'VSM', 'VSM', 'VSM', 'VSM']
        self.NGroupsIn                  = []
        self.GroupsInPathsToMapping     = ['', '', '', '', '', '']
        self.GroupsInSuffix             = '_VS'

        ## Resolution of the Kinetics Data in Output? Array of 'StS' / 'VSM' / 'CGM' of size Syst.NMolecules
        self.MolResolutionOut           = ['VSM', 'VSM', 'VSM', 'VSM', 'VSM', 'VSM']
        self.NGroupsOut                 = []
        self.GroupsOutPathsToMapping    = ['', '', '', '', '', '']
        self.GroupsOut_Flg              = False
        self.GroupsOutWrite_Flg         = False
        self.GroupsOutSuffix            = '_VS'

        ## Reading / Writing Kinetics Data
        self.Read_Flg                   = True
        self.ReadFldr                   = WORKSPACE_PATH + '/Mars_Database/Run_0D/database/'

        ## Writing Kinetics Data
        self.Write_Flg                  = True
        self.WriteFldr                  = WORKSPACE_PATH + '/Mars_Database/Run_0D/database/'
        self.WriteDiss_Flg              = True 
        self.WriteDissInel_Flg          = True 
        self.CorrFactor                 = 1.0
        self.WriteInel_Flg              = True
        self.WriteExch_Flg              = True

        ## Correcting Kinetics Based on Window-Averaging
        self.WindAvrg_Flg               = False
        self.WindAvrgJs                 = 3
        self.WindAvrgVs                 = 2

        ## Writing Arrhenius Files
        self.MaxEntOrPlato              = 1
        self.MinRate                    = 1.e-15
        self.MinNbTs                    = 4
        self.MaxErrArr                  = 1.e-7        



class hdf5(object):

    def __init__(self, WORKSPACE_PATH):

        self.ReadFldr                   = WORKSPACE_PATH + '/Mars_Database/HDF5_Database/'
        self.WriteFolder                = ''
        self.ForceReadDat_Flg           = False
        self.Save_Flg                   = True



class ME(object):

    def __init__(self, WORKSPACE_PATH):

        self.Read_Flg                   = False
        self.ReadFldr                   = WORKSPACE_PATH + '/Mars_Database/Run_0D/'
        self.WriteFolder                = ''
        self.ProcCode                   = '0_1_1_1'
        self.TimeVec                    = np.array([1.e-10, 1.e-8, 1.e-6, 1.e-4])



class inputdata(object):

    def __init__( self, WORKSPACE_PATH, CoarseAIRFldr, PyCoarseAIRFldr ):
        
        self.WORKSPACE_PATH            = WORKSPACE_PATH
        self.CoarseAIRFldr             = CoarseAIRFldr
        self.PyCoarseAIRFldr           = PyCoarseAIRFldr

        self.Rates                     = rates()
        self.Kin                       = kinetics( self.WORKSPACE_PATH )
        self.HDF5                      = hdf5(     self.WORKSPACE_PATH )
        self.ME                        = ME(       self.WORKSPACE_PATH )

        self.OldVersion_Flg            = False
        self.SystNameLong              = 'NaNbNcNd_NASA'

        self.TranVec                   = np.array([20000.0])
        self.T0                        = 300.0
        self.NTran                     = np.size(   self.TranVec )
        self.iTVec                     = np.arange( self.NTran   ) + 1

        self.QCTOutFldr                = self.WORKSPACE_PATH + '/CoarseAIR/N4_VS/Test/'
        self.FinalFldr                 = self.WORKSPACE_PATH + '/Mars_Database/Results/'

        self.PlotShow_Flg              = False

        self.DelRateMat_Flg            = True