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
import os
import sys
import numpy as np

## Paths from where to Import Python Modules
#WORKSPACE_PATH = os.environ['WORKSPACE_PATH']
WORKSPACE_PATH = '/home/venturi/WORKSPACE/'
#WORKSPACE_PATH = '/Users/sventuri/WORKSPACE/'
#CoarseAIRFldr  = os.environ['COARSEAIR_SOURCE_DIR'] 
CoarseAIRFldr    = WORKSPACE_PATH + '/CoarseAIR/coarseair/'
sys.path.insert(0, CoarseAIRFldr  + '/scripts/postprocessing/PyCoarseAIR/Objects/')
sys.path.insert(0, CoarseAIRFldr  + '/scripts/postprocessing/PyCoarseAIR/ChemicalSystems/')
sys.path.insert(0, CoarseAIRFldr  + '/scripts/postprocessing/PyCoarseAIR/Reading/')
sys.path.insert(0, CoarseAIRFldr  + '/scripts/postprocessing/PyCoarseAIR/Writing/')
sys.path.insert(0, CoarseAIRFldr  + '/scripts/postprocessing/PyCoarseAIR/Computing/')
sys.path.insert(0, CoarseAIRFldr  + '/scripts/postprocessing/PyCoarseAIR/Initializing/')

if (len(sys.argv) > 1):
    InputFile = sys.argv[1]
    print("\n[PyCoarseAIR]: Calling PyCoarseAIR with Input File = ", InputFile)
    sys.path.insert(0, InputFile)
else:
    print("[PyCoarseAIR]: Calling PyCoarseAIR without Input File")
    sys.path.insert(0, CoarseAIRFldr + '/scripts/postprocessing/PyCoarseAIR/InputData/')

from Reading           import Read_Levels, Read_PartFuncsAndEnergies, Read_qnsEnBin, Read_Rates_CGQCT
from Writing           import Write_PartFuncsAndEnergies, Write_Kinetics_FromOverall, Write_QSS
from Computing         import Compute_Rates_Thermal_FromOverall, Compute_QSS
from Initializing      import Initialize_Data
from InputData         import inputdata
from ME_Output         import me_output
##--------------------------------------------------------------------------------------------------------------



##==============================================================================================================
InputData                       = inputdata()
InputData.SystNameLong          = 'O3_UMN'

InputData.TranVec               = np.array([20000])
NTran                           = np.size(InputData.TranVec)
InputData.iTVec                 = np.arange(NTran) + 1

InputData.QCTOutFldr            = WORKSPACE_PATH + '/CG-QCT/run_O3_ALL/Test/'
InputData.FinalFldr             = WORKSPACE_PATH + '/Mars_Database/Results/'

InputData.Kin.Read_Flg          = False
InputData.Kin.Write_Flg         = False
InputData.Kin.ReadArr_Flg       = False
InputData.Kin.WriteArr_Flg      = False
InputData.Kin.SaveArr_Flg       = False
InputData.Kin.MinRate           = 1.e-15
InputData.Kin.MinNbTs           = 4
InputData.Kin.ReadFldr          = WORKSPACE_PATH + '/Mars_Database/Run_0D/database/'
InputData.Kin.WriteFldr         = WORKSPACE_PATH + '/Mars_Database/Run_0D/database/'
InputData.Kin.Diss_Flg          = False    
InputData.Kin.CorrFactor        = 5.33333333333
InputData.Kin.Inel_Flg          = False
InputData.Kin.Exch_Flg          = False
InputData.Kin.MaxEntOrPlato     = 1
InputData.Kin.MaxErrArr         = 1.e-7
InputData.Kin.WindAvrgFlg           = False
InputData.Kin.WindAvrgJs            = 3
InputData.Kin.WindAvrgVs            = 2

InputData.HDF5.ReadFldr         = WORKSPACE_PATH + '/Mars_Database/HDF5_Database/'
InputData.HDF5.ForceReadDat_Flg = False
InputData.HDF5.Save_Flg         = True

InputData.ME.Read_Flg           = True
InputData.ME.ProcCode           = '0_1_1_1'
InputData.ME.ReadFldr           = WORKSPACE_PATH + '/Mars_Database/Run_0D/'
InputData.ME.TimeVec            = np.array([1.e-10, 1.e-8, 1.e-6, 1.e-4])

InputData.PlotShow_Flg          = False

InputData.DelRateMat_Flg        = True
##--------------------------------------------------------------------------------------------------------------



print("\n[PyCoarseAIR]: Initializing Data")

[Syst, Temp] = Initialize_Data(InputData)

print("\n[PyCoarseAIR]: Uploading Data")

Syst = Read_Levels(Syst, InputData)
Syst = Read_qnsEnBin(Syst, InputData)
Syst = Read_PartFuncsAndEnergies(Syst, Temp, InputData)
Syst = Read_Rates_CGQCT(Syst, Temp, InputData)


if (InputData.Kin.Write_Flg):
    Write_PartFuncsAndEnergies(Syst, Temp, InputData)
    #Write_Kinetics_FromOverall(Syst, Temp, InputData)


if (InputData.ME.Read_Flg):
    InputData.FinalFldrT = InputData.FinalFldr + '/T' + str(Temp.TranVec[1-1]) + 'K/'

if not os.path.exists(InputData.FinalFldrT):
    os.makedirs(InputData.FinalFldrT)

InputData.ME.ReadFldr = InputData.ME.ReadFldr + 'output_' + InputData.SystNameLong + '_T' + str(InputData.TranVec[1-1]) + 'K_' + InputData.ME.ProcCode + '/'
ME = me_output(Syst)
ME.Read_Box(InputData)
ME.Plot_MolFracs_Evolution(InputData, Temp, 1)
# ME.Plot_TTran_Evolution(InputData, Temp, 1)
# ME.Plot_Rho_Evolution(InputData, Temp, 1)
# ME.Plot_P_Evolution(InputData, Temp, 1)
# ME.Plot_Nd_Evolution(InputData, Temp, 1)
# ME.Plot_Energy_Evolution(InputData, Temp, 1)

for iComp in range(ME.NCFDComp):
    if ( ME.Component[iComp].ToMol > -1 ):
        ME.Component[iComp].Read_Pop(InputData, ME.NTime)
        ME.Component[iComp].Plot_Pop(InputData, Syst, ME, Temp, 1)
        ME.Component[iComp].Compute_DistFunc( Syst, ME.NTime )

        # ME.Component[iComp].Compute_KAveraged( Syst, 1, ME.NTime )
        # Syst = Compute_QSS( Syst, ME, 1 )
        # Write_QSS( Syst, Temp, InputData, 1 )
        # ME.Component[iComp].Plot_KAveraged( InputData, Syst, ME, Temp, 1 )

        ME.Component[iComp].Compute_EAveraged( Syst, 1, ME.NTime )
        ME.Component[iComp].Plot_EAveraged( InputData, Syst, ME, Temp, 1 )
        ME.Component[iComp].Compute_Taus( InputData, Syst, ME, 1)
        ME.Component[iComp].Write_Taus( Temp, InputData, 1)


#Syst = Compute_Rates_Thermal_FromOverall(Syst, Temp, InputData)
