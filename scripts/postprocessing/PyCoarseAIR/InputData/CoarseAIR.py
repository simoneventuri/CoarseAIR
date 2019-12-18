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

class InputData(object):

    TranVec            = np.array([5000, 10000, 12500, 20000])
    iTVec              = np.arange(4) + 1

    OutputFldr         = '/home/venturi/WORKSPACE/CG-QCT/run_CO2_ALL/Test/'

    SystNameLong       = 'CO2_NASA'

    PostprocessingFldr = '/home/venturi/WORKSPACE/MarsAIR/Results/'

    ReadKinFolder      = '/home/libo/MarsRates/O3_UMN/kinetics/'
    WriteKinFolder     = '/home/venturi/WORKSPACE/CG-QCT/run_O3_ALL/Test/RunKonig/database/'
    WriteInelKin_Flg   = False
    WriteExchKin_Flg   = False
    WriteDissKin_Flg   = True 

    PlotShowFlg        = False

    PathToHDF5         = '/home/libo/MarsRates/HDF5_Database/'
    ForceReadDatFlg    = False
    SaveHDF5Flg        = True


    def __init__(self):

        self.TraVec             = TranVec
        self.iTVec              = iTVec
 
        self.OutputFldr         = OutputFldr

        self.SystNameLong       = SystNameLong

        self.KinMthd            = KinMthd
        self.NBins              = NBins

        self.PostprocessingFldr = PostprocessingFldr

        self.ReadKinFolder      = ReadKinFolder
        self.WriteKinFolder     = WriteKinFolder
        self.WriteInelKin_Flg   = WriteInelKin_Flg
        self.WriteExchKin_Flg   = WriteExchKin_Flg
        self.WriteDissKin_Flg   = WriteDissKin_Flg    

        self.PlotShowFlg        = PlotShowFlg

        self.ForceReadDatFlg    = ForceReadDatFlg
        self.SaveHDF5Flg        = SaveHDF5Flg