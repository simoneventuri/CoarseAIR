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
import os, errno
import os.path
from os import path
import shutil

import numpy as np
import pandas
import h5py
import csv

sys.path.insert(0, '../src/Parameters/')
from Parameters import *
  
from Processes  import processes
import MolecularProperties

def mkdirs(newdir, mode=0o777):
    os.makedirs(newdir, mode, exist_ok=True)
    




# ===================================================================================================================== CLASS ===
class grouped_t_properties(object):

    def __init__(self, NProcTypes, T):

        self.Value    = T

        self.EeV      = 0.0
        self.EeV0     = 0.0
        self.QRatio   = 0.0
        self.Q        = 0.0
        self.Q0  = 0.0

        self.Proc     = [processes() for iProc in range(4)]
        self.ProcExch = [processes() for iProc in range(NProcTypes-2)]
        self.ProcTot  = [processes() for iProc in range(NProcTypes)]
# --------------------------------------------------------------------------------------------------------------------- CLASS ---





# ===================================================================================================================== CLASS ===
class groupedmolecule(object):

    def __init__(self, Temp, PathToMapping, T0, NProcTypes, Name, CFDCompName, Type, ToMol ):

        self.Name          = Name
        self.ToMol         = ToMol
        self.CFDCompName   = CFDCompName
        self.Type          = Type
        self.PathToMapping = PathToMapping
        self.NTs           = Temp.NTran+1
        self.TVec          = np.concatenate( (np.array([T0]), Temp.TranVec), axis=0 )
        self.T             = [grouped_t_properties(NProcTypes, self.TVec[iT]) for iT in range(self.NTs)]



    # ***************************************************************************************************************************
    def Initialize( self, InputData, Syst, Temp, NLevels, Levelvqn, LevelEeV, Levelg, LevelWrite_Flg, In_Flg ):
        print('      [Molecule.py - Initialize]: Initializing Grouped Molecule ' + self.Name )

        self.Get_Mapping( Levelvqn )

        self.Compute_GroupProps( NLevels, Levelg, LevelEeV, LevelWrite_Flg )

        if (In_Flg >= 0):

            for iT in Temp.iTVec:
                TTra = Temp.TranVec[iT-1]
                TInt = TTra
                
                self.Save_PartFuncsAndEnergiesAtT_HDF5( InputData, Syst, TTra, TInt, iT, In_Flg )
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Get_Mapping( self, Levelvqn ):
        print('      [Molecule.py - Get_Mapping]: Computing Nb of Groups and Mapping Levels -> Groups for Grouped Molecule ' + self.Name )

        if (self.Type == 'VSM'):
            self.NGroups = np.max(Levelvqn)+1
            self.Mapping = Levelvqn
        
        elif (self.Type == 'CGM'):
            Data   = pandas.read_csv(self.PathToMapping, header=None, skiprows=1)
            Data   = Data.apply(pandas.to_numeric, errors='coerce')
            Idx    = np.array(Data.values[:,0], dtype=np.int64  ) - 1
            Groups = np.array(Data.values[:,1], dtype=np.int64  ) - 1
            self.Mapping         = np.zeros((np.max(Idx)+1), dtype=np.int64) 
            self.Mapping[Idx[:]] = Groups[:]
            self.NGroups         = np.max(Groups)+1
        
        else:
            print('      [Molecule.py - Get_Mapping]:   ERROR! Grouping Method "' + self.Type + '" NOT IMPLEMENTED!' )

        print('      [Molecule.py - Get_Mapping]:   Nb of Groups = ' + str(self.NGroups) )
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Compute_GroupProps( self, NLevels, Levelg, LevelEeV, LevelWrite_Flg ):
        print('      [Molecule.py - Compute_GroupProps]: Computing Group Properties for Grouped Molecule ' + self.Name )

        Flg       = np.zeros((self.NGroups), dtype=np.int8)
        self.EeV0 = np.zeros((self.NGroups))
        for iLevel in range(NLevels):
            if (Flg[self.Mapping[iLevel]] == 0):
                Flg[self.Mapping[iLevel]] = 1
                self.EeV0[self.Mapping[iLevel]] = LevelEeV[iLevel]

        for iT in range(self.NTs):

            self.T[iT].Expvec       = Levelg * np.exp( -                       LevelEeV * Ue / (self.T[iT].Value * UKb) )
            self.T[iT].ExpvecThermo = Levelg * np.exp( - (LevelEeV - np.amin(LevelEeV)) * Ue / (self.T[iT].Value * UKb) )
            self.T[iT].Q            = np.zeros((self.NGroups))
            self.T[iT].Q0           = np.zeros((self.NGroups))
            self.T[iT].EeV          = np.zeros((self.NGroups))
            self.T[iT].QRatio       = np.zeros((self.NGroups))

            for iLevel in range(NLevels):
                if (LevelWrite_Flg[iLevel]):
                    self.T[iT].Q[self.Mapping[iLevel]]   = self.T[iT].Q[self.Mapping[iLevel]]   +                    self.T[iT].Expvec[iLevel]
                    self.T[iT].Q0[self.Mapping[iLevel]]  = self.T[iT].Q0[self.Mapping[iLevel]]  +                    self.T[iT].ExpvecThermo[iLevel]
                    self.T[iT].EeV[self.Mapping[iLevel]] = self.T[iT].EeV[self.Mapping[iLevel]] + LevelEeV[iLevel] * self.T[iT].Expvec[iLevel]

            self.T[iT].EeV    = self.T[iT].EeV / self.T[iT].Q
            self.T[iT].QRatio = self.T[iT].Q   / np.sum(self.T[iT].Q)
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Write_InitialConditions( self, InputData, Syst, Temp, In_Flg ):
        print('    [Molecule.py - Write_InitialConditions]: Writing Initial Mole Fractions for Grouped Molecule: ' + self.CFDCompName )

        if (In_Flg == 1):
            GroupsSuffix = InputData.Kin.GroupsInSuffix
        else:
            GroupsSuffix = InputData.Kin.GroupsOutSuffix

        mkdirs(    InputData.Kin.WriteFldr + '/thermo/' ) 
        mkdirs(    InputData.Kin.WriteFldr + '/thermo/' + Syst.NameLong + GroupsSuffix )        
        TempFldr = InputData.Kin.WriteFldr + '/thermo/' + Syst.NameLong + GroupsSuffix

        MoleFracsFile = TempFldr + '/' + self.CFDCompName + '_InitialMoleFracs_T' + str(int(self.T[0].Value)) + 'K.dat' 
        csvmole       = open(MoleFracsFile, 'w')
        Line          = '# Percentage of ' + self.CFDCompName + ' Contained in Each Group\n'
        csvmole.write(Line)

        for iGroup in range(self.NGroups):
            Line     = '%.10e\n' % float(self.T[0].QRatio[iGroup])
            csvmole.write(Line)

        csvmole.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Write_ThermoProperties( self, InputData, Syst, Temp, In_Flg ):

        if (In_Flg == 1):
            GroupsSuffix = InputData.Kin.GroupsInSuffix
        else:
            GroupsSuffix = InputData.Kin.GroupsOutSuffix

        mkdirs(    InputData.Kin.WriteFldr + '/thermo/' ) 
        mkdirs(    InputData.Kin.WriteFldr + '/thermo/' + Syst.NameLong + GroupsSuffix )      
        TempFldr = InputData.Kin.WriteFldr + '/thermo/' + Syst.NameLong + GroupsSuffix

        for iT in range(1, self.NTs):
            print('      [Molecule.py - Write_ThermoProperties]: Writing Thermo File for Grouped Molecule: ' + self.CFDCompName + ' at T = ' + str(int(self.TVec[iT])) + ' K' )

            PathToFile = TempFldr + '/' + self.CFDCompName + '_T' + str(int(self.TVec[iT])) + 'K'


            if (InputData.Kin.WriteFormat == 'PLATO'):
        
                PathToFileOrig = TempFldr + '/../' + self.CFDCompName + '_Format'
                if not os.path.isfile(PathToFileOrig):
                    Syst.Molecule[self.ToMol].Create_ThermoFile_Format( PathToFileOrig )

                DestTemp       = shutil.copyfile(PathToFileOrig, PathToFile)

                with open(PathToFile, 'a') as f:
                    Line = 'NB_ENERGY_LEVELS = ' + str(self.NGroups) + '\n'
                    f.write(Line)
                    np.savetxt(f, np.transpose(np.array([self.T[iT].Q0, self.T[iT].Q0*0.0])), fmt='%.8e    %.8e')
                f.close()


            elif (InputData.Kin.WriteFormat == 'csv'):
                
                with open(PathToFile, 'w') as f:
                    f.write('# $Q_I$,$E_I~[eV]$\n')
                    np.savetxt(f, np.transpose(np.array([self.T[iT].Q0, self.T[iT].EeV])), fmt='%e,%e')
                f.close()


            elif (InputData.Kin.WriteFormat == 'custom'):
                print('    [System.py - Write_ThermoProperties]: Customize the File by Specifying the Desired Format' )


            else:
                raise NameError("    [System.py - Write_ThermoProperties]: InputData.Kin.WriteFormat must be specified. Select 'PLATO', 'csv', or 'custom'")

    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Save_PartFuncsAndEnergiesAtT_HDF5( self, InputData, Syst, TTra, TInt, iT, In_Flg ):
        print('        [Molecule.py - Save_PartFuncsAndEnergiesAtT_HDF5]: Saving Data')

        if (In_Flg == 1):
            GroupsSuffix = InputData.Kin.GroupsInSuffix
        else:
            GroupsSuffix = InputData.Kin.GroupsOutSuffix


        PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        f          = h5py.File(PathToFile, 'a')

        TStr       = 'T_' + str(int(TTra)) + '_' + str(int(TInt))
        TempStr1   = TStr + '/' + self.CFDCompName + GroupsSuffix
        print('        [Molecule.py - Save_PartFuncsAndEnergiesAtT_HDF5]: Looking in ' + TempStr1)
        TempStr2   = TempStr1 + '/GroupEeV'     
        if TempStr2 in f.keys():
            print('        [Molecule.py - Save_PartFuncsAndEnergiesAtT_HDF5]: Overwriting the Level Properties in the HDF5 File')

            grp       = f[TempStr1]
            Data      = grp["GroupEeV"]
            Data[...] = self.T[iT-1].EeV
            Data      = grp["GroupQ"]
            Data[...] = self.T[iT-1].Q
            Data      = grp["GroupQ0"]
            Data[...] = self.T[iT-1].Q0
            Data      = grp["GroupQRatio"]
            Data[...] = self.T[iT-1].QRatio

        else:
            print('        [Molecule.py - Save_PartFuncsAndEnergiesAtT_HDF5]: Level Properties NOT Found in the HDF5 File. Saving them for the FIRST Time')

            if (not (TempStr1 in f.keys())):
                grp       = f.create_group(TempStr1)
            else:
                grp       = f[TempStr1]
            GroupEeV     = grp.create_dataset("GroupEeV",    data=self.T[iT-1].EeV,    compression="gzip", compression_opts=9)
            GroupQ       = grp.create_dataset("GroupQ",      data=self.T[iT-1].Q,      compression="gzip", compression_opts=9)
            GroupQ0      = grp.create_dataset("GroupQ0",     data=self.T[iT-1].Q0,     compression="gzip", compression_opts=9)
            GroupQRatio  = grp.create_dataset("GroupQRatio", data=self.T[iT-1].QRatio, compression="gzip", compression_opts=9)

        f.close()
    # ...........................................................................................................................



# --------------------------------------------------------------------------------------------------------------------- CLASS ---





# ===================================================================================================================== CLASS ===
class t_properties(object):

    def __init__(self):

        self.LevelEeV = 0.0
        self.QRatio   = 0.0
        self.Q        = 0.0
        self.QExp     = 0.0
# --------------------------------------------------------------------------------------------------------------------- CLASS ---





# ===================================================================================================================== CLASS ===
class molecule(object):

    def __init__(self, NTTran):

        self.Name              = ""
        self.DissEn            = np.array([1, 1, 1])
        self.DegeneracyFactor  = 0
        self.Mu                = 0

        self.Nins             = 0
        self.KinMthd           = ""

        self.CFDCompName       = ""

        self.StartLevel        = 0
        self.FinalLevel        = 0

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
        self.LevelToGroupIn      = 0

        self.NTTran            = NTTran
        self.T                 = [t_properties() for iTTran in range(NTTran)]



    # ***************************************************************************************************************************
    def Initialize( self, InputData, Syst, Temp, iMol ):

        self.CFDCompName = Syst.CFDComp[Syst.MolToCFDComp[iMol]].Name


        GetProperties       = getattr( MolecularProperties, 'GetProperties_' + self.CFDCompName )
        Syst.Molecule[iMol] = GetProperties( self )


        print('    [Molecule.py - Initialize]: Reading the Levels for the Molecule ' + self.Name )
        self.Read_Levels( InputData, Syst )


        self.KinMthdIn  = InputData.Kin.MolResolutionIn[iMol]
        if ( self.KinMthdIn  == 'StS' ):  
            self.EqNStatesIn = self.NLevels 
            print('    [Molecule.py - Initialize]:   We are Starting from a StS Molecule ' + self.Name + '. Nb of Levels = ' + str(self.EqNStatesIn) )
        elif ( self.KinMthdIn  == 'VSM' ):  
            self.EqNStatesIn = self.Nvqn
            print('    [Molecule.py - Initialize]:   We are Starting from a Vibrational Specific Molecule ' + self.Name + '. Nb of Groups = ' + str(self.EqNStatesIn) )
        elif ( self.KinMthdIn  == 'CGM' ): 
            self.EqNStatesIn = InputData.NGroupsIn[iMol]
            print('    [Molecule.py - Initialize]:   We are Starting from a Coarse-Grained Molecule ' + self.Name + '. Nb of Groups = ' + str(self.EqNStatesIn) )

        self.KinMthdOut  = InputData.Kin.MolResolutionOut[iMol]
        if ( self.KinMthdOut  == 'StS' ):  
            self.EqNStatesOut = self.NLevels 
            print('    [Molecule.py - Initialize]:   We are Ending with a StS Molecule ' + self.Name + '. Nb of Levels = ' + str(self.EqNStatesOut) )
        elif ( self.KinMthdOut  == 'VSM' ):  
            self.EqNStatesOut = self.Nvqn
            print('    [Molecule.py - Initialize]:   We are Ending with a Vibrational Specific Molecule ' + self.Name + '. Nb of Groups = ' + str(self.EqNStatesOut) )
        elif ( self.KinMthdOut  == 'CGM' ): 
            self.EqNStatesOut = InputData.NGroupsOut[iMol]
            print('    [Molecule.py - Initialize]:   We are Ending with a Coarse-Grained Molecule ' + self.Name + '. Nb of Groups = ' + str(self.EqNStatesOut) )
             

        self.Read_qnsEnBin( InputData, Syst, Temp )

        
        if ( (not self.KinMthdIn  == 'StS') ):            
            print('    [Molecule.py - Initialize]:     Initializing the Intial Grouped Molecule')
            self.GroupsIn = groupedmolecule(Temp, InputData.Kin.GroupsInPathsToMapping[iMol], InputData.T0, Syst.NProcTypes, self.Name, self.CFDCompName, self.KinMthdIn, iMol  )
            self.GroupsIn.Initialize( InputData, Syst, Temp, self.NLevels, self.Levelvqn, self.LevelEeV, self.Levelg, self.LevelWrite_Flg, 1 )
            
            for iT in Temp.iTVec:
                self.T[iT-1].EqEeV0In = self.GroupsIn.T[iT-1].EeV

        else:
            for iT in Temp.iTVec:
                self.T[iT-1].EqEeV0In = self.LevelEeV0


        if (not self.KinMthdOut == 'StS'):   
            print('    [Molecule.py - Initialize]:      Initializing the Final Grouped Molecule')
            self.GroupsOut = groupedmolecule(Temp, InputData.Kin.GroupsOutPathsToMapping[iMol], InputData.T0, Syst.NProcTypes, self.Name, self.CFDCompName, self.KinMthdOut, iMol )        
            self.GroupsOut.Initialize( InputData, Syst, Temp, self.NLevels, self.Levelvqn, self.LevelEeV, self.Levelg, self.LevelWrite_Flg, 0 )

            self.GroupsOut.Write_InitialConditions( InputData, Syst, Temp, 0 )
            
            self.GroupsOut.Write_ThermoProperties( InputData, Syst, Temp, 0 )

            for iT in Temp.iTVec:
                self.T[iT-1].EqEeV0Out = self.GroupsOut.T[iT-1].EeV

        else:
            print('    [Molecule.py - Initialize]:    We Desire a StS Molecule ' + self.Name )  

            self.Write_InitialConditions( InputData, Syst, Temp )
            
            self.Write_ThermoProperties( InputData, Syst, Temp )

            for iT in Temp.iTVec:
                self.T[iT-1].EqEeV0Out = self.LevelEeV0


        if (InputData.Kin.PackUnpackDiss_Flg):
            print('    [Molecule.py - Initialize]:    We Desire To Pack the StS Dissociation Rates And then Unpacking Them. Initializing.')
            self.PackUnpack = groupedmolecule( Temp, InputData.Kin.PackUnpackPathsToMapping[iMol], InputData.T0, Syst.NProcTypes, self.Name, self.CFDCompName, InputData.Kin.PackUnpackType[iMol], iMol  )
            self.PackUnpack.Initialize( InputData, Syst, Temp, self.NLevels, self.Levelvqn, self.LevelEeV, self.Levelg, self.LevelWrite_Flg, -1 )


        print('    [Molecule.py - Initialize]: Done Initializing the Molecule ' + self.Name )
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Read_Levels( self, InputData, Syst ):
        print('      [Molecule.py - Read_Levels]: Reading Level Properties')

        PathToFile  = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        HDF5Exist_Flg = path.exists(PathToFile)
        if (HDF5Exist_Flg):
            f = h5py.File(PathToFile, 'a')
            TempStr     = self.Name + '/LevelrMin'
            PresentFlg  = TempStr in f.keys()
            f.close()
        else:
            f = {'key': 'value'}
            PresentFlg = False
        if (PresentFlg): print('      [Molecule.py - Read_Levels]:   Found Data for Molecule ' + self.Name + ' in the HDF5 File')

        if ( (PresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):
            print('      [Molecule.py - Read_Levels]:     Loading Data for Molecule ' + self.Name + ' from the HDF5 File')
            self.Load_Levels_HDF5(Syst)

        else:
            PathToFile = Syst.PathToReadFldr + '/' + self.Name + '/levels_cut.inp'
            print('      [Molecule.py - Read_Levels]:     Reading File: ', PathToFile)
            
            Data = pandas.read_csv(PathToFile, header=None, skiprows=15, delimiter=r"\s+")
            Data = Data.apply(pandas.to_numeric, errors='coerce')

            self.Levelvqn   = np.array(Data.values[:,0],  dtype=np.int16  )
            self.Leveljqn   = np.array(Data.values[:,1],  dtype=np.int16  )
            self.LevelEEh   = np.array(Data.values[:,2],  dtype=np.float64)
            self.LevelEgam  = np.array(Data.values[:,3],  dtype=np.float64)
            self.LevelrMin  = np.array(Data.values[:,4],  dtype=np.float64)
            self.LevelrMax  = np.array(Data.values[:,5],  dtype=np.float64)
            self.LevelVMin  = np.array(Data.values[:,6],  dtype=np.float64)
            self.LevelVMax  = np.array(Data.values[:,7],  dtype=np.float64)
            self.LevelTau   = np.array(Data.values[:,8],  dtype=np.float64)
            self.LevelrIn   = np.array(Data.values[:,9],  dtype=np.float64)
            self.LevelrOut  = np.array(Data.values[:,10], dtype=np.float64)
            
            self.LevelEeV   = self.LevelEEh * Hartree_To_eV 
            self.NLevels    = np.size(self.LevelEeV,0)

            if (InputData.HDF5.Save_Flg):
                print('      [Molecule.py - Read_Levels]:       Saving Data for Molecule ' + self.Name + ' in the HDF5 File')
                self.Save_Levels_HDF5( Syst )
        
        self.LevelEeV0  = (self.LevelEEh - np.amin(self.LevelEEh )) * Hartree_To_eV 
        self.DissEeV    = (                np.amin(self.LevelVMax)) * Hartree_To_eV 
        self.Nvqn       = np.max(self.Levelvqn)+1
        self.Njqn       = np.max(self.Leveljqn)+1


        self.LevelWrite_Flg  = np.zeros((self.NLevels), dtype=bool)
        self.LevelNewMapping = np.zeros((self.NLevels), dtype=np.int64) - 1
        jLevel = -1
        for iLevel in range (self.NLevels):
            if ( (self.LevelEeV[iLevel] < self.DissEeV) and (InputData.Kin.WriteQB_IntFlg != 1) ) or ( (self.LevelEeV[iLevel] > self.DissEeV) and (InputData.Kin.WriteQB_IntFlg != 0) ):
                jLevel                       = jLevel + 1
                self.LevelWrite_Flg[iLevel]  = True
                self.LevelNewMapping[iLevel] = jLevel 
        self.NLevelsNewMapping = jLevel + 1 


        self.Compute_ERot()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Read_qnsEnBin( self, InputData, Syst, Temp ):
        print('      [Molecule.py - Read_qnsEnBin]: Reading Level To Group Mapping')

       
        PathToFile  = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        f           = h5py.File(PathToFile, 'a')
        TempStr     = self.Name + '/Levelg'
        PresentFlg  = TempStr in f.keys()
        f.close()
        if (PresentFlg): print('      [Molecule.py - Read_qnsEnBin]:   Found Data for Molecule ' + self.Name + ' in the HDF5 File')

        if ( (PresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):
            print('      [Molecule.py - Read_qnsEnBin]:   Loading Data for Molecule ' + self.Name + ' from HDF5 File')
            self.Load_qnsEnBin_HDF5( InputData, Syst )

        else:
            if (InputData.OldVersion_IntFlg > 0):
                PathToFile = Syst.PathToReadFldr + '/' + self.Name + '/' + self.Name + '_' + str(self.EqNStatesIn) + '/qnsEnBin.dat'
                Data = pandas.read_csv(PathToFile, header=None, skiprows=1, delimiter=r"\s+")
            else:
                PathToFile = Syst.PathToReadFldr + '/' + self.Name + '/Bins_' + str(self.EqNStatesIn) + '/QNsEnBin.csv'
                Data = pandas.read_csv(PathToFile, header=None, skiprows=1)
            print('      [Molecule.py - Read_qnsEnBin]:   Reading File: ', PathToFile)
            
            Data = Data.apply(pandas.to_numeric, errors='coerce')

            self.Levelg         = np.array(Data.values[:,4], dtype=np.float64)
            self.LevelToGroupIn = np.array(Data.values[:,5], dtype=np.int32  )

            if (InputData.HDF5.Save_Flg):
                print('      [Molecule.py - Read_qnsEnBin]:   Saving Data for Molecule ' + self.Name + ' in the HDF5 File')
                self.Save_qnsEnBin_HDF5( InputData, Syst )


        for iT in Temp.iTVec:
            TTra = Temp.TranVec[iT-1]
            TInt = TTra

            TStr       = 'T_' + str(int(TTra)) + '_' + str(int(TInt))
            TempStr    = TStr + '/' + self.CFDCompName
            PresentFlg = TempStr in f.keys()

            if ( (PresentFlg) and (not InputData.HDF5.ForceReadDat_Flg) ):

                print('      [Molecule.py - Read_qnsEnBin]:   Loading Data for Temperature Nb ' + str(iT) + ' (T = ' + str(int(TTra)) + 'K) in the HDF5 File')
                self.Load_PartFuncsAndEnergiesAtT_HDF5( Syst, iT, TTra, TInt )
            
            else:
                    
                #self.T[iT-1].LevelEeV    = self.LevelEEh * Hartree_To_eV
                self.T[iT-1].LevelQE      = np.exp( - self.LevelEeV * Ue / (TTra * UKb) )
                self.T[iT-1].LevelQERatio = self.T[iT-1].LevelQE / np.sum(self.T[iT-1].LevelQE)
                self.T[iT-1].LevelQ       = self.Levelg * np.exp( - self.LevelEeV * Ue / (TTra * UKb) )
                self.T[iT-1].LevelQRatio  = self.T[iT-1].LevelQ / np.sum(self.T[iT-1].LevelQ)

                if (InputData.HDF5.Save_Flg):
                    print('      [Molecule.py - Read_qnsEnBin]:   Saving Data in the HDF5 File')
                    self.Save_PartFuncsAndEnergiesAtT_HDF5( Syst, TTra, TInt, iT )
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Load_Levels_HDF5( self, Syst ):
        print('        [Molecule.py - Load_Levels_HDF5]: Loading Data')

        PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        HDF5Exist_Flg = path.exists(PathToFile)
        if (HDF5Exist_Flg):
            f = h5py.File(PathToFile, 'a')
        else:
            f = {'key': 'value'}

        TempStr = self.CFDCompName   
        grp     = f[TempStr]

        Data                          = grp["NLevels"]
        self.NLevels   = Data[...][0]
        Data                          = grp["Levelvqn"]
        self.Levelvqn  = Data[...]
        Data                          = grp["Leveljqn"]
        self.Leveljqn  = Data[...]
        Data                          = grp["LevelEEh"]
        self.LevelEEh  = Data[...]
        Data                          = grp["LevelEeV"]
        self.LevelEeV  = Data[...]
        Data                          = grp["LevelEgam"]
        self.LevelEgam = Data[...]
        Data                          = grp["LevelrMin"]
        self.LevelrMin = Data[...]
        Data                          = grp["LevelrMax"]
        self.LevelrMax = Data[...]
        Data                          = grp["LevelVMin"]
        self.LevelVMin = Data[...]
        Data                          = grp["LevelVMax"]
        self.LevelVMax = Data[...]
        Data                          = grp["LevelTau"]
        self.LevelTau  = Data[...]
        Data                          = grp["LevelrIn"]
        self.LevelrIn  = Data[...]
        Data                          = grp["LevelrOut"]
        self.LevelrOut = Data[...]

        f.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Load_qnsEnBin_HDF5( self, InputData, Syst ):
        print('        [Molecule.py - Load_qnsEnBin_HDF5]: Loading Data')

        PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        HDF5Exist_Flg = path.exists(PathToFile)
        if (HDF5Exist_Flg):
            f = h5py.File(PathToFile, 'a')
        else:
            f = {'key': 'value'}

        TempStr = self.CFDCompName
        grp     = f[TempStr]

        Data                = grp["Levelg"]
        self.Levelg         = Data[...]
        StrTemp             = "LevelToGroupIn" + InputData.Kin.GroupsInSuffix
        Data                = grp[StrTemp]
        self.LevelToGroupIn = Data[...]

        f.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Load_PartFuncsAndEnergiesAtT_HDF5( self, Syst, iT, TTra, TInt ):
        print('        [Molecule.py - Load_PartFuncsAndEnergiesAtT_HDF5]: Loading Data')

        PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        f          = h5py.File(PathToFile, 'a')

        TStr    = 'T_' + str(int(TTra)) + '_' + str(int(TInt))
        TempStr = TStr + '/' + self.CFDCompName
        grp     = f[TempStr]

        Data                                 = grp["LevelQRatio"]
        self.T[iT-1].LevelQRatio  = Data[...]
        Data                                 = grp["LevelQ"]
        self.T[iT-1].LevelQ       = Data[...] 
        Data                                 = grp["LevelQERatio"]
        self.T[iT-1].LevelQERatio = Data[...]
        Data                                 = grp["LevelQE"]
        self.T[iT-1].LevelQE      = Data[...] 

        f.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Save_Levels_HDF5( self, Syst ):
        print('        [Molecule.py - Save_Levels_HDF5]: Saving Data')

        PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        HDF5Exist_Flg = path.exists(PathToFile)
        if (HDF5Exist_Flg):
            f = h5py.File(PathToFile, 'a')
        else:
            f = h5py.File(PathToFile, 'w')

        TempStr1   = self.CFDCompName
        TempStr2   = TempStr1 + '/LevelrMin'
        if TempStr2 in f.keys():
            print('        [Molecule.py - Save_Levels_HDF5]:   Overwriting the Level Properties in the HDF5 File')

            grp       = f[TempStr1]
            Data      = grp["NLevels"]
            Data[...] = np.array([self.NLevels])
            Data      = grp["Levelvqn"]
            Data[...] = self.Levelvqn
            Data      = grp["Leveljqn"]
            Data[...] = self.Leveljqn
            Data      = grp["LevelEEh"]
            Data[...] = self.LevelEEh
            Data      = grp["LevelEeV"]
            Data[...] = self.LevelEeV
            Data      = grp["LevelEgam"]
            Data[...] = self.LevelEgam
            Data      = grp["LevelrMin"]
            Data[...] = self.LevelrMin
            Data      = grp["LevelrMax"]
            Data[...] = self.LevelrMax
            Data      = grp["LevelVMin"]
            Data[...] = self.LevelVMin
            Data      = grp["LevelVMax"]
            Data[...] = self.LevelVMax
            Data      = grp["LevelTau"]
            Data[...] = self.LevelTau
            Data      = grp["LevelrIn"]
            Data[...] = self.LevelrIn
            Data      = grp["LevelrOut"]
            Data[...] = self.LevelrOut

        else:
            print('        [Molecule.py - Save_Levels_HDF5]:   Level Properties NOT Found in the HDF5 File. Saving them for the FIRST Time')

            if (not (TempStr1 in f.keys())):
                grp       = f.create_group(TempStr1)
            else:
                grp       = f[TempStr1]
            NLevels   = grp.create_dataset("NLevels",   data=np.array([self.NLevels]), compression="gzip", compression_opts=9)
            Levelvqn  = grp.create_dataset("Levelvqn",  data=self.Levelvqn,            compression="gzip", compression_opts=9)
            Leveljqn  = grp.create_dataset("Leveljqn",  data=self.Leveljqn,            compression="gzip", compression_opts=9)
            LevelEEh  = grp.create_dataset("LevelEEh",  data=self.LevelEEh,            compression="gzip", compression_opts=9)
            LevelEeV  = grp.create_dataset("LevelEeV",  data=self.LevelEeV,            compression="gzip", compression_opts=9)
            LevelEgam = grp.create_dataset("LevelEgam", data=self.LevelEgam,           compression="gzip", compression_opts=9)
            LevelrMin = grp.create_dataset("LevelrMin", data=self.LevelrMin,           compression="gzip", compression_opts=9)
            LevelrMax = grp.create_dataset("LevelrMax", data=self.LevelrMax,           compression="gzip", compression_opts=9)
            LevelVMin = grp.create_dataset("LevelVMin", data=self.LevelVMin,           compression="gzip", compression_opts=9)
            LevelVMax = grp.create_dataset("LevelVMax", data=self.LevelVMax,           compression="gzip", compression_opts=9)
            LevelTau  = grp.create_dataset("LevelTau",  data=self.LevelTau,            compression="gzip", compression_opts=9)
            LevelrIn  = grp.create_dataset("LevelrIn",  data=self.LevelrIn,            compression="gzip", compression_opts=9)
            LevelrOut = grp.create_dataset("LevelrOut", data=self.LevelrOut,           compression="gzip", compression_opts=9)

        f.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Save_qnsEnBin_HDF5( self, InputData, Syst ):
        print('        [Molecule.py - Save_qnsEnBin_HDF5]: Saving Data')

        PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        #HDF5Exist_Flg = path.exists(PathToFile)
        #if (HDF5Exist_Flg):
        f = h5py.File(PathToFile, 'a')
        #else:
        #    f = {'key': 'value'}
       
        StrTemp   = "LevelToGroupIn" + InputData.Kin.GroupsInSuffix

        TempStr1   = self.CFDCompName
        TempStr2   = TempStr1 + '/Levelg'
        if TempStr2 in f.keys():
            print('        [Molecule.py - Save_qnsEnBin_HDF5]: Overwriting the Level Properties in the HDF5 File')

            grp       = f[TempStr1]
            Data      = grp["Levelg"]
            Data[...] = self.Levelg
            Data      = grp[StrTemp]
            Data[...] = self.LevelToGroupIn

        else:
            print('        [Molecule.py - Save_qnsEnBin_HDF5]: Level Properties NOT Found in the HDF5 File. Saving them for the FIRST Time')
            if (not (TempStr1 in f.keys())):
                grp        = f.create_group(TempStr1)
            else:
                grp        = f[TempStr1]
            Levelg         = grp.create_dataset("Levelg", data=self.Levelg,         compression="gzip", compression_opts=9)
            LevelToGroupIn = grp.create_dataset(StrTemp,  data=self.LevelToGroupIn, compression="gzip", compression_opts=9)

        f.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Save_PartFuncsAndEnergiesAtT_HDF5( self, Syst, TTra, TInt, iT ):
        print('        [Molecule.py - Save_PartFuncsAndEnergiesAtT_HDF5]: Saving Data')

        PathToFile = Syst.PathToHDF5 + '/' + Syst.NameLong + '.hdf5'
        #HDF5Exist_Flg = path.exists(PathToFile)
        #if (HDF5Exist_Flg):
        f = h5py.File(PathToFile, 'a')
        #else:
        #    f = {'key': 'value'}

        TStr     = 'T_' + str(int(TTra)) + '_' + str(int(TInt))
        TempStr1 = TStr + '/' + self.CFDCompName + '/'
        print('        [Molecule.py - Save_PartFuncsAndEnergiesAtT_HDF5]: Looking into ' + TempStr1 )
        TempStr2 = TempStr1 + '/LevelQ'     
        if TempStr2 in f.keys():
            print('        [Molecule.py - Save_PartFuncsAndEnergiesAtT_HDF5]: Overwriting the Level Properties in the HDF5 File')

            grp       = f[TempStr1]
            Data      = grp["LevelQ"]
            Data[...] = self.T[iT-1].LevelQ
            Data      = grp["LevelQRatio"]
            Data[...] = self.T[iT-1].LevelQRatio
            Data      = grp["LevelQE"]
            Data[...] = self.T[iT-1].LevelQE
            Data      = grp["LevelQERatio"]
            Data[...] = self.T[iT-1].LevelQERatio

        else:
            print('        [Molecule.py - Save_PartFuncsAndEnergiesAtT_HDF5]: Level Properties NOT Found in the HDF5 File. Saving them for the FIRST Time')

            if (not (TempStr1 in f.keys())):
                grp      = f.create_group(TempStr1)
            else:
                grp      = f[TempStr1]
            LevelQ       = grp.create_dataset("LevelQ",       data=self.T[iT-1].LevelQ,       compression="gzip", compression_opts=9)
            LevelQRatio  = grp.create_dataset("LevelQRatio",  data=self.T[iT-1].LevelQRatio,  compression="gzip", compression_opts=9)
            LevelQE      = grp.create_dataset("LevelQE",      data=self.T[iT-1].LevelQE,      compression="gzip", compression_opts=9)
            LevelQERatio = grp.create_dataset("LevelQERatio", data=self.T[iT-1].LevelQERatio, compression="gzip", compression_opts=9)

        f.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Compute_ERot( self ):
        print('        [Molecule.py - Compute_ERot]: Computing Vibrational and Rotational Energies from the Ro-Vibrational Ones (By Prioritizing Vibration)')

        self.vEeV              = np.zeros(np.max(self.Levelvqn)+1)
        self.LevelEVib         = np.zeros(self.NLevels)
        self.LevelERot         = np.zeros(self.NLevels)
        for iLevel in range(self.NLevels):
            if (self.Leveljqn[iLevel] == 0):
                self.vEeV[int(self.Levelvqn[iLevel])] = self.LevelEeV[iLevel]
            self.LevelEVib[iLevel] = self.vEeV[int(self.Levelvqn[iLevel])]
            self.LevelERot[iLevel] = self.LevelEeV[iLevel] - self.LevelEVib[iLevel]
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Write_InitialConditions( self, InputData, Syst, Temp ):

        print('    [Molecule.py - Write_InitialConditions]: Computing Ratios of Initial Mole Fractions for Molecule: ' + self.CFDCompName )
        self.Compute_QRatio0(Temp.T0)

        mkdirs(    InputData.Kin.WriteFldr + '/thermo/' ) 
        mkdirs(    InputData.Kin.WriteFldr + '/thermo/' + Syst.NameLong + InputData.Kin.GroupsOutSuffix )   
        TempFldr = InputData.Kin.WriteFldr + '/thermo/' + Syst.NameLong + InputData.Kin.GroupsOutSuffix
        
        MoleFracsFile = TempFldr + '/' + self.CFDCompName + '_InitialMoleFracs_T' + str(int(Temp.T0)) + 'K.dat' 
        print('    [Molecule.py - Write_InitialConditions]: Writing Initial Mole Fractions for Molecule: ' + self.CFDCompName )
        csvmole       = open(MoleFracsFile, 'w')
        Line          = '# Percentage of ' + self.CFDCompName + ' Contained in Each Level\n'
        csvmole.write(Line)

        for iLevel in range(self.NLevels):
            if (self.LevelWrite_Flg[iLevel]):
                jLevel = self.LevelNewMapping[iLevel]
                Line     = '%.10e\n' % float(np.maximum( self.QRatio0[iLevel], 1.e-99+self.QRatio0[iLevel]*0.0 ))
                csvmole.write(Line)

        csvmole.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Write_ThermoProperties( self, InputData, Syst, Temp ):

        mkdirs(    InputData.Kin.WriteFldr + '/thermo/' ) 
        mkdirs(    InputData.Kin.WriteFldr + '/thermo/' + Syst.NameLong + InputData.Kin.GroupsOutSuffix )   
        TempFldr = InputData.Kin.WriteFldr + '/thermo/' + Syst.NameLong + InputData.Kin.GroupsOutSuffix

        for iT in Temp.iTVec:
            print('      [Molecule.py - Write_ThermoProperties]: Writing Thermo File for Molecule: ' + self.CFDCompName + ' at T = ' + str(int(Temp.TranVec[iT-1])) + ' K' )

            PathToFile     = TempFldr + '/'    + self.CFDCompName + '_T' + str(int(Temp.TranVec[iT-1])) + 'K'


            if (InputData.Kin.WriteFormat == 'PLATO'):

                PathToFileOrig = TempFldr + '/../' + self.CFDCompName + '_Format'
                if not os.path.isfile(PathToFileOrig):
                    self.Create_ThermoFile_Format( PathToFileOrig )

                DestTemp  = shutil.copyfile(PathToFileOrig, PathToFile)
                print('      [Molecule.py - Write_ThermoProperties]: Copied File to: ', DestTemp)

                with open(PathToFile, 'a') as f:
                    Line = 'NB_ENERGY_LEVELS = ' + str(self.NLevelsNewMapping) + '\n'
                    f.write(Line)
                    #np.savetxt(f, np.transpose(np.array([self.Levelg, self.LevelEeV0])), fmt='%.8e    %.8e')
                f.close()

                csvmole = open(PathToFile, 'a')
                for iLevel in range(self.NLevels):
                    if (self.LevelWrite_Flg[iLevel]):
                        jLevel = self.LevelNewMapping[iLevel]
                        Line   = '%.8e    %.8e\n' % ( float(self.Levelg[iLevel]), float(self.LevelEeV0[iLevel]) )
                        csvmole.write(Line)
                csvmole.close()


            elif (InputData.Kin.WriteFormat == 'csv'):
                
                csvmole = open(PathToFile, 'w')
                csvmole.write('# $g_i$,$E_i~[eV]$\n')
                for iLevel in range(self.NLevels):
                    if (self.LevelWrite_Flg[iLevel]):
                        jLevel = self.LevelNewMapping[iLevel]
                        Line   = '%e,%e\n' % ( float(self.Levelg[iLevel]), float(self.LevelEeV0[iLevel]) )
                        csvmole.write(Line)
                csvmole.close()


            elif (InputData.Kin.WriteFormat == 'custom'):
                print('    [System.py - Write_ThermoProperties]: Customize the File by Specifying the Desired Format' )


            else: 
                raise NameError("    [System.py - Write_ThermoProperties]: InputData.Kin.WriteFormat must be specified. Select 'PLATO', 'csv', or 'custom'")

    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Create_ThermoFile_Format( self, PathToFileOrig ):
        print('    [Molecule.py - Create_ThermoFile_Format]: Creating Thermo File Format for Molecule: ' + self.CFDCompName )

        ThermoFile = PathToFileOrig
        csvthemo   = open(ThermoFile, 'w')

        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('#Species file for ' + self.CFDCompName +  '\n')
        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('NAME = '            + self.CFDCompName +  '\n')
        csvthemo.write('                                           \n')
        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('#Molecular mass [kg/mol]                   \n')
        csvthemo.write('MOLAR_MASS = %.8e\n' % self.MolarMass         )
        csvthemo.write('                                           \n')
        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('#Number of atoms                           \n')
        csvthemo.write('NB_ATOMS = 2                               \n')
        csvthemo.write('                                           \n')
        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('#Chemical elements (symbols and quantities)\n')
        csvthemo.write('ELEM_NAME = '  + self.ElementNames +      '\n')
        csvthemo.write('ELEM_QUANT = ' + self.ElementQnts  +      '\n')
        csvthemo.write('                                           \n')
        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('#Formation energy [J/mol]                  \n')
        csvthemo.write('EF = %.8e\n' % self.FormationE                )
        csvthemo.write('                                           \n')
        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('#Electric charge                           \n')
        csvthemo.write('CHARGE = 0                                 \n')
        csvthemo.write('                                           \n')
        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('#Linearity                                 \n')
        csvthemo.write('LIN = %d\n' % self.LinFactor                  )
        csvthemo.write('                                           \n')
        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('#Symmetry factor                           \n')
        csvthemo.write('SYM = %.2e\n' % self.SymmFactor               )
        csvthemo.write('                                           \n')
        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('#Thermal model for internal energy         \n')
        csvthemo.write('MODEL = NONE                               \n')
        csvthemo.write('                                           \n')
        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('#Characteristic rotational temperature [K] \n')
        csvthemo.write('THETA_ROT = %.8e\n' % self.ThetaRot           )
        csvthemo.write('                                           \n')
        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('#Characteristic vibrational temperature [K]\n')
        csvthemo.write('THETA_VIB = %.8e\n' % self.ThetaVib           )
        csvthemo.write('                                           \n')
        csvthemo.write('#------------------------------------------\n')
        csvthemo.write('#Levels degeneracies and energies [eV]     \n')

        csvthemo.close()
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Compute_QRatio0( self, T0 ):
        print('      [Molecule.py - Compute_QRatio0]: Computing Molecule Distribution Function at Initial Temperature')

        self.QRatio0 = self.Levelg  * np.exp( - self.LevelEeV0 * Ue / (T0 * UKb) )
        self.QRatio0 = self.QRatio0 / np.sum(self.QRatio0)
    # ...........................................................................................................................



    # ***************************************************************************************************************************
    def Initialize_WindAvrg_Objects( self, InputData, Syst):
        print('      [Molecule.py - Initialize_WindAvrg_Objects]: Initializing Objects for Kinetics Window-Averaging')
        
        self.WindAvrgDet   = int( (2*InputData.Kin.WindAvrgJs+1) * (2*InputData.Kin.WindAvrgVs+1) )
        self.WindAvrgMat   = np.zeros( (self.NLevels,  self.WindAvrgDet), dtype=np.int64)
        self.WindAvrgFound = np.zeros( (self.NLevels,  1), dtype=np.int64 )

        for iLevel in range(self.NLevels):
            iFound = -1
            for jLevel in range(self.NLevels):
                if (np.absolute(self.Levelvqn[iLevel] - self.Levelvqn[jLevel]) <= InputData.Kin.WindAvrgVs) and  (np.absolute(self.Leveljqn[iLevel] - self.Leveljqn[jLevel]) <= InputData.Kin.WindAvrgJs):
                    iFound     = iFound  + 1
                    self.WindAvrgMat[iLevel,iFound] = jLevel
            self.WindAvrgFound[iLevel]              = iFound 
        
        print('      [Molecule.py - Initialize_WindAvrg_Objects]:   Done Initializing Objects for Kinetics Window-Averaging\n')
    # ...........................................................................................................................



# --------------------------------------------------------------------------------------------------------------------- CLASS ---