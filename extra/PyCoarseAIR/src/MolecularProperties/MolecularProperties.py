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

#from Molecule import molecule


################################################################################################################
### N2
def GetProperties_N2( Molecule ):   

    Molecule.ElementNames = 'N'
    Molecule.ElementQnts  = '2'
    Molecule.MolarMass    = 28.0134e-3
    Molecule.FormationE   = 0.0
    Molecule.LinFactor    = 1
    Molecule.SymmFactor   = 2.0
    Molecule.ThetaRot     = 2.88
    Molecule.ThetaVib     = 3392.7

    return Molecule


################################################################################################################
### O2
def GetProperties_O2( Molecule ):   

    Molecule.ElementNames = 'O'
    Molecule.ElementQnts  = '2'
    Molecule.MolarMass    = 31.9988e-3
    Molecule.FormationE   = 0.0
    Molecule.LinFactor    = 1
    Molecule.SymmFactor   = 2.0
    Molecule.ThetaRot     = 2.08
    Molecule.ThetaVib     = 2273.6

    return Molecule


################################################################################################################
### CO
def GetProperties_CO( Molecule ):   

    Molecule.ElementNames = 'C O'
    Molecule.ElementQnts  = '1 1'
    Molecule.MolarMass    = 28.0104e-3
    Molecule.FormationE   = -113811.0e+00
    Molecule.LinFactor    = 1
    Molecule.SymmFactor   = 1.0
    Molecule.ThetaRot     = 2.78
    Molecule.ThetaVib     = 3121.9

    return Molecule

    