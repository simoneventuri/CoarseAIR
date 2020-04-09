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
import numpy as np
import pandas
import h5py
from os import path

sys.path.insert(0, '../Saving/')
sys.path.insert(0, '../Loading/')
sys.path.insert(0, '../Computing/')
sys.path.insert(0, '../Parameters/')
sys.path.insert(0, '../Writing/')
sys.path.insert(0, '../Plotting/')

from Saving        import Save_Levels_HDF5, Save_qnsEnBin_HDF5, Save_PartFuncsAndEnergiesAtT_HDF5, Save_RatesAtT_HDF5_3Atoms, Save_RatesAtT_HDF5_4Atoms
from Loading       import Load_Levels_HDF5, Load_qnsEnBin_HDF5, Load_PartFuncsAndEnergiesAtT_HDF5, Load_RatesAtT_HDF5_3Atoms, Load_RatesAtT_HDF5_4Atoms, Load_DissRatesAtT_HDF5
from Computing     import Compute_Correction_To_DissRates, Compute_Rates_Overall, Compute_Rates_Thermal, Compute_BackwardRates, Compute_PrefJumps, Compute_WindAvrg_Matrix, Compute_Correction_To_GroupedDissRates
from Parameters    import *
from Writing       import Write_Rates_Thermal, Write_Kinetics_3Atoms, Write_Kinetics_4Atoms, Write_GroupedKinetics, Write_PrefJumps
from Plotting      import Plot_Rates_Thermal




#########################################################################################################################################################
#########################################################################################################################################################