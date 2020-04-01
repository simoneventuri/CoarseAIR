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
import pandas
import h5py



def Save_Arr_Diss_HDF5(Syst, iMol):

    PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    f          = h5py.File(PathToFile, 'a')

    TempStr1   = Syst.Molecule[iMol].Name
    
    TempStr2   = TempStr1 + '/Arr/Diss/'

    if TempStr2 in f.keys():

        #del f[TempStr2]
        grp       = f[TempStr2]
        Data      = grp["Idxs"]
        Data[...] = Syst.Arr.Proc[0].Idxs
        Data      = grp["Coeffs"]
        Data[...] = Syst.Arr.Proc[0].Coeffs

    else:
        if (not (TempStr2 in f.keys())):
            grp       = f.create_group(TempStr2)
        else:
            grp       = f[TempStr2]
        Idxs   = grp.create_dataset("Idxs",   data=Syst.Arr.Proc[0].Idxs,   compression="gzip", compression_opts=9)
        Coeffs = grp.create_dataset("Coeffs", data=Syst.Arr.Proc[0].Coeffs, compression="gzip", compression_opts=9)

    f.close()



def Save_Arr_Inel_HDF5(Syst, iMol):

    PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    f          = h5py.File(PathToFile, 'a')

    TempStr1   = Syst.Molecule[iMol].Name
    
    TempStr2   = TempStr1 + '/Arr/Inel/'

    if TempStr2 in f.keys():

        del f[TempStr2]
        grp       = f.create_group(TempStr2)
        Idxs   = grp.create_dataset("Idxs",   data=Syst.Arr.Proc[1].Idxs,   compression="gzip", compression_opts=9)
        Coeffs = grp.create_dataset("Coeffs", data=Syst.Arr.Proc[1].Coeffs, compression="gzip", compression_opts=9)

    else:
        if (not (TempStr2 in f.keys())):
            grp       = f.create_group(TempStr2)
        else:
            grp       = f[TempStr2]
        Idxs   = grp.create_dataset("Idxs",   data=Syst.Arr.Proc[1].Idxs,   compression="gzip", compression_opts=9)
        Coeffs = grp.create_dataset("Coeffs", data=Syst.Arr.Proc[1].Coeffs, compression="gzip", compression_opts=9)

    f.close()



def Save_Arr_Exch_HDF5(Syst, iMol):

    PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
    f          = h5py.File(PathToFile, 'a')

    TempStr1   = Syst.Molecule[iMol].Name
    
    TempStr2   = TempStr1 + '/Arr/Exch/'

    if TempStr2 in f.keys():

        del f[TempStr2]
        grp       = f.create_group(TempStr2)
        Idxs   = grp.create_dataset("Idxs",   data=Syst.Arr.Proc[2].Idxs,   compression="gzip", compression_opts=9)
        Coeffs = grp.create_dataset("Coeffs", data=Syst.Arr.Proc[2].Coeffs, compression="gzip", compression_opts=9)

    else:
        if (not (TempStr2 in f.keys())):
            grp       = f.create_group(TempStr2)
        else:
            grp       = f[TempStr2]
        Idxs   = grp.create_dataset("Idxs",   data=Syst.Arr.Proc[2].Idxs,   compression="gzip", compression_opts=9)
        Coeffs = grp.create_dataset("Coeffs", data=Syst.Arr.Proc[2].Coeffs, compression="gzip", compression_opts=9)

    f.close()