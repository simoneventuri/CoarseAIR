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

from Atom      import atom
from Molecule  import molecule
from Pair      import pair
from CFDComp   import cfdcomp
from Processes import processes
from QSS       import qss


class t_properties(object):

    def __init__(self, NProcTypes):

        self.Proc     = [processes() for iProc in range(4)]
        self.ProcExch = [processes() for iProc in range(NProcTypes-2)]
        self.ProcTot  = [processes() for iProc in range(NProcTypes)]
        self.QSS      = qss(NProcTypes)


class arrhenius(object):

    def __init__(self, NProcTypes):

        self.Proc    = [processes() for iProc in range(NProcTypes)]
        self.MinRate = 0.0
        self.MinNbTs = 3


class system(object):

    def __init__(self, SystNameLong, SystName, NAtoms, NMolecules, NPairs, NCFDComp, NTTran, NProcTypes):

        self.NameLong     = SystNameLong
        self.Name         = SystName

        self.NAtoms       = NAtoms
        self.Atom         = [atom() for iA in range(NAtoms)]

        self.NMolecules   = NMolecules
        self.Molecule     = [molecule(NTTran) for iMol in range(NMolecules)]

        self.NPairs       = NPairs
        self.Pair         = [pair() for iP in range(NPairs)]

        self.NCFDComp     = NCFDComp
        self.CFDComp      = [cfdcomp() for iComp in range(NCFDComp)]
        self.MolToCFDComp = 0

        self.PathToFolder = ''
        self.PathToHDF5   = ''

        self.NProcTypes   = NProcTypes

        self.T            = [t_properties(NProcTypes) for iT in range(NTTran)]

        self.Arr          = arrhenius(NProcTypes)

        self.Proc         = [processes() for iProc in range(NProcTypes)]

        self.RatesTh      = np.zeros((NTTran, NProcTypes))
        self.RatesQSS     = np.zeros((NTTran, NProcTypes))

        self.MolToCFDComp = np.zeros((NMolecules), dtype=np.int8)