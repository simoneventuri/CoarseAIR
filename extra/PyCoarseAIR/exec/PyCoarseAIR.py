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
##==============================================================================================================
print("\n[PyCoarseAIR]: Defining Paths")

WORKSPACE_PATH  = '/home/venturi/WORKSPACE/'                                                                    # It is NOT REQUIRED to change this path
#WORKSPACE_PATH  = os.environ['WORKSPACE_PATH']   
CoarseAIRFldr   = WORKSPACE_PATH + '/CoarseAIR/coarseair/'										         	    # It is NOT REQUIRED to change this path
#CoarseAIRFldr   = os.environ['COARSEAIR_SOURCE_DIR']
PyCoarseAIRFldr = CoarseAIRFldr  + '/extra/PyCoarseAIR/'		                							    # <--- Please, CHANGE this path to the folder of the PyCoarseAIR Code

DtbHDF5Fldr     = WORKSPACE_PATH + '/Mars_Database/HDF5_Database/'                                              # <--- Please, CHANGE this path to the folder containing the HDF5 File

DtbWriteFldr    = WORKSPACE_PATH + '/Mars_Database/Run_0D/database/'                                            # <--- Please, CHANGE this path to the folder where you want PYCoareseAIR writing the Postprocessed Kinetics Database
OutputWriteFldr = WORKSPACE_PATH + '/Mars_Database/Results/'                                                    # <--- Please, CHANGE this path to the folder where you want PYCoareseAIR writing its Output Files and Plots
##--------------------------------------------------------------------------------------------------------------



##==============================================================================================================
print("\n[PyCoarseAIR]: Loading Files")

sys.path.insert(0, PyCoarseAIRFldr  + '/src/Objects/')
sys.path.insert(0, PyCoarseAIRFldr  + '/src/ChemicalSystems/')
sys.path.insert(0, PyCoarseAIRFldr  + '/src/MolecularProperties/')
sys.path.insert(0, PyCoarseAIRFldr  + '/src/Reading/')
sys.path.insert(0, PyCoarseAIRFldr  + '/src/Writing/')
sys.path.insert(0, PyCoarseAIRFldr  + '/src/Computing/')
sys.path.insert(0, PyCoarseAIRFldr  + '/src/Initializing/')

if (len(sys.argv) > 1):
    InputFile = sys.argv[1]
    print("[PyCoarseAIR]:   Calling PyCoarseAIR with Input File = ", InputFile)
    sys.path.insert(0, InputFile)
else:
    InputFile = PyCoarseAIRFldr + '/src/InputData/'
    print("[PyCoarseAIR]:   Calling PyCoarseAIR with the PRESET Input File Located in " + InputFile )
    sys.path.insert(0, InputFile)
##--------------------------------------------------------------------------------------------------------------



##==============================================================================================================
print("\n[PyCoarseAIR]: Loading Functions")

from Initializing      import Initialize_Data
from InputData         import inputdata
##--------------------------------------------------------------------------------------------------------------



##==============================================================================================================
print("\n[PyCoarseAIR]: Initializing Input and System Data")

InputData    = inputdata(WORKSPACE_PATH, CoarseAIRFldr, PyCoarseAIRFldr, DtbHDF5Fldr, DtbWriteFldr, OutputWriteFldr)

[Syst, Temp] = Initialize_Data(InputData)

print("[PyCoarseAIR]: Done Initializing")
##--------------------------------------------------------------------------------------------------------------



##==============================================================================================================
print("\n[PyCoarseAIR]: Uploading and Processing Rates Data")

Syst.Read_Rates( InputData, Temp )

#Syst.Write_Kinetics_FromOverall( InputData, Temp )

print("[PyCoarseAIR]: Done Initializing")
##--------------------------------------------------------------------------------------------------------------



# if (InputData.ME.Read_Flg):

#     InputData.FinalFldrT = InputData.FinalFldr + '/T' + str(Temp.TranVec[1-1]) + 'K/'
#     if not os.path.exists(InputData.FinalFldrT):
#         os.makedirs(InputData.FinalFldrT)
    
#     InputData.ME.ReadFldr = InputData.ME.ReadFldr + 'output_' + InputData.SystNameLong + '_T' + str(InputData.TranVec[1-1]) + 'K_' + InputData.ME.ProcCode + '/'

#     ME = me_output(Syst)
#     ME.Read_Box(InputData)
#     ME.Plot_MolFracs_Evolution(InputData, Temp, 1)
#     # ME.Plot_TTran_Evolution(InputData, Temp, 1)
#     # ME.Plot_Rho_Evolution(InputData, Temp, 1)
#     # ME.Plot_P_Evolution(InputData, Temp, 1)
#     # ME.Plot_Nd_Evolution(InputData, Temp, 1)
#     # ME.Plot_Energy_Evolution(InputData, Temp, 1)

#     for iComp in range(ME.NCFDComp):
#         if ( ME.Component[iComp].ToMol > -1 ):
#             ME.Component[iComp].Read_Pop(InputData, ME.NTime)
#             ME.Component[iComp].Plot_Pop(InputData, Syst, ME, Temp, 1)
#             ME.Component[iComp].Compute_DistFunc( Syst, ME.NTime )
            
#             #ME.Component[iComp].Compute_KAveraged( Syst, 1, ME.NTime )
#             #Syst = Compute_QSS( Syst, ME, 1 )
#             #Write_QSS( Syst, Temp, InputData, 1 )
#             #ME.Component[iComp].Plot_KAveraged( InputData, Syst, ME, Temp, 1 )

#             ME.Component[iComp].Compute_EAveraged( Syst, 1, ME.NTime )
#             ME.Component[iComp].Plot_EAveraged( InputData, Syst, ME, Temp, 1 )
#             ME.Component[iComp].Compute_Taus( InputData, Syst, ME, 1)
#             ME.Component[iComp].Write_Taus( Temp, InputData, 1)


# #Syst = Compute_Rates_Thermal_FromOverall(Syst, Temp, InputData)
