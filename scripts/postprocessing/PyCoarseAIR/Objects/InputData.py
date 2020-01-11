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


class kinetics(object):

    def __init__(self):

        self.Read_Flg         = False
        self.Write_Flg        = False
        self.ReadFldr         = ''
        self.WriteFldr        = ''
        self.WriteInel_Flg    = False
        self.WriteExch_Flg    = False
        self.WriteDiss_Flg    = False
        self.CorrFactor       = 1.0 


class hdf5(object):

    def __init__(self):

        self.ReadFolder       = ''
        self.WriteFolder      = ''
        self.ForceReadDat_Flg = False
        self.Save_Flg         = False



class ME(object):

    def __init__(self):

        self.Read_Flg         = False
        self.ReadFolder       = ''
        self.WriteFolder      = ''



class inputdata(object):

    def __init__(self):
        self.Kin              = kinetics()
        self.HDF5             = hdf5()
        self.ME               = ME()
