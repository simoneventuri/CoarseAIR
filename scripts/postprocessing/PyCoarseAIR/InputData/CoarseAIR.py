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

    # TranVec            = np.array([1500, 3000, 6000, 8000, 10000, 12000, 14000])
    # iTVec              = np.arange(7) + 1

    # OutputFldr         = '/home/venturi/WORKSPACE/CG-QCT/run_O3_ALL/Test/'

    # SystName           = 'O3'

    # KinMthd            = ['StS']
    # NBins              = np.array([6115])

    # PostprocessingFldr = '/home/venturi/WORKSPACE/MarsAIR/Results/'

    # ReadKinFolder      = '/home/libo/MarsRates/O3_UMN/kinetics/'
    # WriteKinFolder     = '/home/venturi/WORKSPACE/CG-QCT/run_O3_ALL/Test/RunKonig/database/'
    # WriteInelKin_Flg   = False
    # WriteExchKin_Flg   = False
    # WriteDissKin_Flg   = True 

    # PlotShowFlg        = False


    TranVec            = np.array([10000])
    iTVec              = np.arange(1) + 1

    OutputFldr         = '/Users/sventuri/WORKSPACE/CG-QCT/run_O2C/Test/'

    SystName           = 'O2C'

    KinMthd            = ['StS', 'StS']
    NBins              = np.array([6078, 13521])

    PostprocessingFldr = '/Users/sventuri/WORKSPACE/MarsAIR/Results/'

    ReadKinFolder      = '/home/libo/MarsRates/O3_UMN/kinetics/'
    WriteKinFolder     = '/Users/sventuri/WORKSPACE/CG-QCT/run_O2C/Test/RunKonig/database/'
    WriteInelKin_Flg   = False
    WriteExchKin_Flg   = False
    WriteDissKin_Flg   = True 

    PlotShowFlg        = False


    def __init__(self):

        self.TraVec             = TranVec
        self.iTVec              = iTVec
 
        self.OutputFldr         = OutputFldr

        self.SystName           = SystName

        self.KinMthd            = KinMthd
        self.NBins              = NBins

        self.PostprocessingFldr = PostprocessingFldr

        self.ReadKinFolder      = ReadKinFolder
        self.WriteKinFolder     = WriteKinFolder
        self.WriteInelKin_Flg   = WriteInelKin_Flg
        self.WriteExchKin_Flg   = WriteExchKin_Flg
        self.WriteDissKin_Flg   = WriteDissKin_Flg    

        self.PlotShowFlg        = PlotShowFlg