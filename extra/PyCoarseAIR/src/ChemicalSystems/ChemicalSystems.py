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

from System import system


################################################################################################################
### CO + O System (From NASA Ames, Dr. D. Schwenke)
def CO2_NASA_Upload( Temp ):   

    SystNameLong = 'CO2_NASA'
    SystName      = 'CO2'          
    

    NAtoms        = 3
    NMolecules    = 2
    NDistMolecules= 2
    NPairs        = 3
    NCFDComp      = 4
    NProcTypes    = 4
    Syst          = system(SystNameLong, SystName, NAtoms, NMolecules, NDistMolecules, NPairs, NCFDComp, Temp.NTran, NProcTypes)


    Syst.Atom[0].Name  = 'C'
    Syst.Atom[1].Name  = 'O'
    Syst.Atom[2].Name  = 'O'

    Syst.Atom[0].Color = np.array([0, 0, 0])
    Syst.Atom[1].Color = np.array([0, 0, 1])
    Syst.Atom[2].Color = np.array([0, 0, 1])

    Syst.Atom[0].Size  = 100
    Syst.Atom[1].Size  = 200
    Syst.Atom[2].Size  = 200

    Syst.Atom[0].Mass  = 21868.661757
    Syst.Atom[1].Mass  = 29148.94559
    Syst.Atom[2].Mass  = 29148.94559



    Syst.Molecule[0].Name             = 'CO'
    Syst.Molecule[1].Name             = 'O2'

    Syst.Molecule[0].DissEn           = 0.0
    Syst.Molecule[1].DissEn           = 0.0

    Syst.Molecule[0].DegeneracyFactor = 1
    Syst.Molecule[1].DegeneracyFactor = 6

    Syst.Molecule[0].Mu               = 28.0104e-3 
    Syst.Molecule[1].Mu               = 31.9988e-3

    #Syst.Molecule[0].NLevels          = 13521
    #Syst.Molecule[1].NLevels          = 6078


    Syst.Pair[0].Name  = 'CO'
    Syst.Pair[1].Name  = 'CO'
    Syst.Pair[2].Name  = 'O2'

    Syst.Pair[0].ToMol = 0
    Syst.Pair[1].ToMol = 0
    Syst.Pair[2].ToMol = 1

    Syst.Pair[0].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[1].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[2].Color = np.array([0, 0, 256])  / 256


    Syst.CFDComp[0].Name   = 'C'
    Syst.CFDComp[1].Name   = 'O'
    Syst.CFDComp[2].Name   = 'CO'
    Syst.CFDComp[3].Name   = 'O2'

    Syst.CFDComp[0].ToMol   = -1
    Syst.CFDComp[1].ToMol   = -1
    Syst.CFDComp[2].ToMol   =  0
    Syst.CFDComp[3].ToMol   =  1

    Syst.CFDComp[0].Mass    = Syst.Atom[0].Mass
    Syst.CFDComp[1].Mass    = Syst.Atom[1].Mass
    Syst.CFDComp[2].Mass    = Syst.Atom[0].Mass+Syst.Atom[1].Mass
    Syst.CFDComp[3].Mass    = 2.0*Syst.Atom[2].Mass

    Syst.CFDComp[0].Deg     = 1
    Syst.CFDComp[1].Deg     = 9
    Syst.CFDComp[2].Deg     = 1
    Syst.CFDComp[3].Deg     = 1

    Syst.CFDComp[0].Color   = np.array([ 102, 102, 102]) / 256
    Syst.CFDComp[1].Color   = np.array([   0,   0,   0]) / 256
    Syst.CFDComp[2].Color   = np.array([ 204,   0,   0]) / 256
    Syst.CFDComp[3].Color   = np.array([   0,   0, 234]) / 256

    Syst.CFDComp[0].LineStyle = ':'
    Syst.CFDComp[1].LineStyle = '-.'
    Syst.CFDComp[2].LineStyle = '-'
    Syst.CFDComp[3].LineStyle = '--'

    Syst.CFDComp[0].RxLxIdx = -1
    Syst.CFDComp[1].RxLxIdx = -1
    Syst.CFDComp[2].RxLxIdx =  1
    Syst.CFDComp[3].RxLxIdx =  0

    Syst.MolToCFDComp       = np.array([2, 3])

    Syst.ExchtoMol          = np.array([0,1])
    Syst.ExchtoAtom         = np.array([2,0])

    Syst.ColPartToComp      = 0

    return Syst
################################################################################################################


################################################################################################################
### CO + O System (From NASA Ames, Dr. D. Schwenke)
def O2C_NASA_Upload( Temp ):   

    SystNameLong  = 'O2C_NASA'        
    SystName      = 'CO2'          
    

    NAtoms        = 3
    NMolecules    = 2
    NDistMolecules= 2
    NPairs        = 3
    NCFDComp      = 4
    NProcTypes    = 3
    Syst          = system(SystNameLong, SystName, NAtoms, NMolecules, NDistMolecules, NPairs, NCFDComp, Temp.NTran, NProcTypes)


    Syst.Atom[0].Name  = 'O'
    Syst.Atom[1].Name  = 'O'
    Syst.Atom[2].Name  = 'C'

    Syst.Atom[0].Color = np.array([0, 0, 1])
    Syst.Atom[1].Color = np.array([0, 0, 1])
    Syst.Atom[2].Color = np.array([0, 0, 0])

    Syst.Atom[0].Size  = 200
    Syst.Atom[1].Size  = 200
    Syst.Atom[2].Size  = 100

    Syst.Atom[0].Mass  = 29148.94559
    Syst.Atom[1].Mass  = 29148.94559
    Syst.Atom[2].Mass  = 21868.661757


    Syst.Molecule[0].Name             = 'O2'
    Syst.Molecule[1].Name             = 'CO'

    Syst.Molecule[0].DissEn           = 0.0
    Syst.Molecule[1].DissEn           = 0.0

    Syst.Molecule[0].DegeneracyFactor = 6
    Syst.Molecule[1].DegeneracyFactor = 1

    Syst.Molecule[0].Mu               = 31.9988e-3
    Syst.Molecule[1].Mu               = 28.0104e-3 

    #Syst.Molecule[0].NLevels          =  6078
    #Syst.Molecule[1].NLevels          = 13521


    Syst.Pair[0].Name  = 'O2'
    Syst.Pair[1].Name  = 'CO'
    Syst.Pair[2].Name  = 'CO'

    Syst.Pair[0].ToMol = 0
    Syst.Pair[1].ToMol = 1
    Syst.Pair[2].ToMol = 1

    Syst.Pair[0].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[1].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[2].Color = np.array([0, 0, 256])  / 256


    Syst.CFDComp[0].Name   = 'C'
    Syst.CFDComp[1].Name   = 'O'
    Syst.CFDComp[2].Name   = 'CO'
    Syst.CFDComp[3].Name   = 'O2'

    Syst.CFDComp[0].ToMol   = -1
    Syst.CFDComp[1].ToMol   = -1
    Syst.CFDComp[2].ToMol   =  0
    Syst.CFDComp[3].ToMol   =  1

    Syst.CFDComp[0].Mass    = Syst.Atom[2].Mass
    Syst.CFDComp[1].Mass    = Syst.Atom[0].Mass
    Syst.CFDComp[2].Mass    = Syst.Atom[1].Mass+Syst.Atom[2].Mass
    Syst.CFDComp[3].Mass    = 2.0*Syst.Atom[0].Mass

    Syst.CFDComp[0].Deg     = 1
    Syst.CFDComp[1].Deg     = 9
    Syst.CFDComp[2].Deg     = 1
    Syst.CFDComp[3].Deg     = 1

    Syst.CFDComp[0].Color   = np.array([ 102, 102, 102]) / 256
    Syst.CFDComp[1].Color   = np.array([   0,   0,   0]) / 256
    Syst.CFDComp[2].Color   = np.array([ 204,   0,   0]) / 256
    Syst.CFDComp[3].Color   = np.array([   0,   0, 234]) / 256

    Syst.CFDComp[0].LineStyle = ':'
    Syst.CFDComp[1].LineStyle = '-.'
    Syst.CFDComp[2].LineStyle = '-'
    Syst.CFDComp[3].LineStyle = '--'

    Syst.CFDComp[0].RxLxIdx = -1
    Syst.CFDComp[1].RxLxIdx = -1
    Syst.CFDComp[2].RxLxIdx =  1
    Syst.CFDComp[3].RxLxIdx =  0

    Syst.MolToCFDComp       = np.array([3, 2])

    Syst.ExchtoMol          = np.array([1])
    Syst.ExchtoAtom         = np.array([0])

    Syst.ColPartToComp      = 1

    return Syst
################################################################################################################


################################################################################################################
### CH + N System (From UIUC)
def CHN_UIUC_Upload( Temp ):   

    SystNameLong  = 'CHN_UIUC' 
    SystName      = 'CHN'          
    

    NAtoms        = 3
    NMolecules    = 3
    NDistMolecules= 3
    NPairs        = 3
    NCFDComp      = 6
    NProcTypes    = 4
    Syst          = system(SystNameLong, SystName, NAtoms, NMolecules, NDistMolecules, NPairs, NCFDComp, Temp.NTran, NProcTypes)


    Syst.Atom[0].Name  = 'C'
    Syst.Atom[1].Name  = 'H'
    Syst.Atom[2].Name  = 'N'

    Syst.Atom[0].Color = np.array([0, 0, 0])
    Syst.Atom[1].Color = np.array([0, 0, 1])
    Syst.Atom[2].Color = np.array([0, 1, 0])

    Syst.Atom[0].Size  = 200
    Syst.Atom[1].Size  = 100
    Syst.Atom[2].Size  = 150

    Syst.Atom[0].Mass  = 21868.661757
    Syst.Atom[1].Mass  = 1835.0397616
    Syst.Atom[2].Mass  = 25519.042285


    Syst.Molecule[0].Name             = 'CH'
    Syst.Molecule[1].Name             = 'CN'
    Syst.Molecule[2].Name             = 'HN'

    Syst.Molecule[0].DissEn           = 0.0
    Syst.Molecule[1].DissEn           = 0.0
    Syst.Molecule[2].DissEn           = 0.0

    Syst.Molecule[0].DegeneracyFactor = 2
    Syst.Molecule[1].DegeneracyFactor = 2
    Syst.Molecule[1].DegeneracyFactor = 3

    Syst.Molecule[0].Mu               = 15.01454e-3 
    Syst.Molecule[1].Mu               = 13.01854e-3
    Syst.Molecule[2].Mu               = 26.0174e-3

    #Syst.Molecule[0].NLevels          = 0
    #Syst.Molecule[1].NLevels          = 0
    #Syst.Molecule[2].NLevels          = 0

    Syst.Pair[0].Name  = 'CH'
    Syst.Pair[1].Name  = 'CN'
    Syst.Pair[2].Name  = 'HN'

    Syst.Pair[0].ToMol = 0
    Syst.Pair[1].ToMol = 1
    Syst.Pair[2].ToMol = 2

    Syst.Pair[0].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[1].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[2].Color = np.array([0, 0, 256])  / 256


    Syst.CFDComp[0].Name   = 'C'
    Syst.CFDComp[1].Name   = 'H'
    Syst.CFDComp[2].Name   = 'N'
    Syst.CFDComp[1].Name   = 'CH'
    Syst.CFDComp[1].Name   = 'CN'
    Syst.CFDComp[3].Name   = 'HN'

    Syst.CFDComp[0].ToMol   = -1
    Syst.CFDComp[1].ToMol   = -1
    Syst.CFDComp[2].ToMol   = -1
    Syst.CFDComp[3].ToMol   =  0
    Syst.CFDComp[4].ToMol   =  1
    Syst.CFDComp[5].ToMol   =  2

    Syst.CFDComp[0].Mass    = Syst.Atom[0].Mass
    Syst.CFDComp[1].Mass    = Syst.Atom[1].Mass
    Syst.CFDComp[2].Mass    = Syst.Atom[2].Mass
    Syst.CFDComp[3].Mass    = Syst.Atom[0].Mass + Syst.Atom[1].Mass
    Syst.CFDComp[4].Mass    = Syst.Atom[0].Mass + Syst.Atom[2].Mass
    Syst.CFDComp[5].Mass    = Syst.Atom[1].Mass + Syst.Atom[2].Mass

    Syst.CFDComp[0].Deg     = 1
    Syst.CFDComp[1].Deg     = 1
    Syst.CFDComp[2].Deg     = 1
    Syst.CFDComp[3].Deg     = 1
    Syst.CFDComp[4].Deg     = 1
    Syst.CFDComp[5].Deg     = 1


    Syst.CFDComp[0].Color   = np.array([ 102, 102, 102]) / 256
    Syst.CFDComp[1].Color   = np.array([   0, 153, 102]) / 256
    Syst.CFDComp[2].Color   = np.array([ 204,   0,   0]) / 256
    Syst.CFDComp[3].Color   = np.array([   0,   0, 234]) / 256
    Syst.CFDComp[4].Color   = np.array([ 204,   0,   0]) / 256
    Syst.CFDComp[5].Color   = np.array([   0,   0, 234]) / 256

    Syst.CFDComp[0].LineStyle = ':'
    Syst.CFDComp[1].LineStyle = '-.'
    Syst.CFDComp[2].LineStyle = '--'
    Syst.CFDComp[3].LineStyle = '-.'
    Syst.CFDComp[4].LineStyle = '-'
    Syst.CFDComp[5].LineStyle = '--'

    Syst.CFDComp[0].RxLxIdx = -1
    Syst.CFDComp[1].RxLxIdx = -1
    Syst.CFDComp[2].RxLxIdx = -1
    Syst.CFDComp[3].RxLxIdx =  1
    Syst.CFDComp[4].RxLxIdx =  1
    Syst.CFDComp[5].RxLxIdx =  1

    Syst.MolToCFDComp       = np.array([3, 4, 5])

    Syst.ColPartToComp      = 2

    return Syst
################################################################################################################



################################################################################################################
### O2 + O System (From UMN, Prof. D. Truhlar)
def O3_UMN_Upload( Temp ):   

    SystNameLong  = 'O3_UMN' 
    SystName      = 'O3'          
    

    NAtoms        = 3
    NMolecules    = 1
    NDistMolecules= 1
    NPairs        = 3
    NCFDComp      = 2
    NProcTypes    = 3
    Syst          = system(SystNameLong, SystName, NAtoms, NMolecules, NDistMolecules, NPairs, NCFDComp, Temp.NTran, NProcTypes)


    Syst.Atom[0].Name  = 'O'
    Syst.Atom[1].Name  = 'O'
    Syst.Atom[2].Name  = 'O'

    Syst.Atom[0].Color = np.array([0, 0, 1])
    Syst.Atom[1].Color = np.array([0, 0, 1])
    Syst.Atom[2].Color = np.array([0, 0, 1])

    Syst.Atom[0].Size  = 200
    Syst.Atom[1].Size  = 200
    Syst.Atom[2].Size  = 200

    Syst.Atom[0].Mass  = 29148.94559
    Syst.Atom[1].Mass  = 29148.94559
    Syst.Atom[2].Mass  = 29148.94559


    Syst.Molecule[0].Name             = 'O2'
    Syst.Molecule[0].DissEn           = 5.2
    Syst.Molecule[0].DegeneracyFactor = 6
    Syst.Molecule[0].Mu               = 31.9988e-3
    #Syst.Molecule[0].NLevels          = 6115


    Syst.Pair[0].Name  = 'O2'
    Syst.Pair[1].Name  = 'O2'
    Syst.Pair[2].Name  = 'O2'

    Syst.Pair[0].ToMol = 0
    Syst.Pair[1].ToMol = 0
    Syst.Pair[2].ToMol = 0

    Syst.Pair[0].Color = np.array([0, 0, 256]) / 256
    Syst.Pair[1].Color = np.array([0, 0, 256]) / 256
    Syst.Pair[2].Color = np.array([0, 0, 256]) / 256


    Syst.CFDComp[0].Name   = 'O'
    Syst.CFDComp[1].Name   = 'O2'

    Syst.CFDComp[0].ToMol   = -1
    Syst.CFDComp[1].ToMol   =  0

    Syst.CFDComp[0].Mass    = Syst.Atom[0].Mass
    Syst.CFDComp[1].Mass    = Syst.Atom[1].Mass+Syst.Atom[2].Mass

    Syst.CFDComp[0].Deg     = 9
    Syst.CFDComp[1].Deg     = 1

    Syst.CFDComp[0].Color   = np.array([   0,   0,   0]) / 256
    Syst.CFDComp[1].Color   = np.array([   0,   0, 234]) / 256

    Syst.CFDComp[0].LineStyle = ':'
    Syst.CFDComp[1].LineStyle = '-'

    Syst.CFDComp[0].RxLxIdx = -1
    Syst.CFDComp[1].RxLxIdx =  0

    Syst.MolToCFDComp       = np.array([1])


    Syst.ExchtoMol          = np.array([0])
    Syst.ExchtoAtom         = np.array([2])

    Syst.ColPartToComp      = 0

    return Syst
################################################################################################################


################################################################################################################
### N2 + N System (From NASA Ames, Dr. R. L. Jaffe & Dr. D. Schwenke)
def N3_NASA_Upload( Temp ):   

    SystNameLong  = 'N3_NASA'
    SystName      = 'N3'          
    

    NAtoms        = 3
    NMolecules    = 1
    NDistMolecules= 1
    NPairs        = 3
    NCFDComp      = 2
    NProcTypes    = 2
    Syst          = system(SystNameLong, SystName, NAtoms, NMolecules, NDistMolecules, NPairs, NCFDComp, Temp.NTran, NProcTypes)


    Syst.Atom[0].Name  = 'N'
    Syst.Atom[1].Name  = 'N'
    Syst.Atom[2].Name  = 'N'

    Syst.Atom[0].Color = np.array([0, 0, 1])
    Syst.Atom[1].Color = np.array([0, 0, 1])
    Syst.Atom[2].Color = np.array([0, 0, 1])

    Syst.Atom[0].Size  = 200
    Syst.Atom[1].Size  = 200
    Syst.Atom[2].Size  = 200

    Syst.Atom[0].Mass  = 25526.04298
    Syst.Atom[1].Mass  = 25526.04298
    Syst.Atom[2].Mass  = 25526.04298


    Syst.Molecule[0].Name             = 'N2'
    Syst.Molecule[0].DissEn           = 0.0
    Syst.Molecule[0].DegeneracyFactor = 1
    Syst.Molecule[0].Mu               = 31.9988e-3
    #Syst.Molecule[0].NLevels          = 9092


    Syst.Pair[0].Name  = 'N2'
    Syst.Pair[1].Name  = 'N2'
    Syst.Pair[2].Name  = 'N2'

    Syst.Pair[0].ToMol = 0
    Syst.Pair[1].ToMol = 0
    Syst.Pair[2].ToMol = 0

    Syst.Pair[0].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[1].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[2].Color = np.array([17, 17, 17]) / 256
    
    Syst.CFDComp[0].Name   = 'N'
    Syst.CFDComp[1].Name   = 'N2'

    Syst.CFDComp[0].ToMol   = -1
    Syst.CFDComp[1].ToMol   =  0

    Syst.CFDComp[0].Mass    = Syst.Atom[0].Mass
    Syst.CFDComp[1].Mass    = 2.0*Syst.Atom[2].Mass

    Syst.CFDComp[0].Deg     = 1
    Syst.CFDComp[1].Deg     = 9

    Syst.CFDComp[0].Color   = np.array([ 102, 102, 102]) / 256
    Syst.CFDComp[1].Color   = np.array([   0,   0,   0]) / 256

    Syst.CFDComp[0].LineStyle = ':'
    Syst.CFDComp[1].LineStyle = '-.'

    Syst.CFDComp[0].RxLxIdx = -1
    Syst.CFDComp[1].RxLxIdx =  1

    Syst.MolToCFDComp       = np.array([1])

    Syst.ExchtoMol          = np.array([0])
    Syst.ExchtoAtom         = np.array([2])

    Syst.ColPartToComp      = 0

    return Syst
################################################################################################################


################################################################################################################
### N2 + N2 System (From NASA Ames, Dr. D. Schwenke)
def N4_NASA_Upload( Temp ):   

    SystNameLong = 'N4_NASA'
    SystName      = 'N4'          
    

    NAtoms        = 4
    NMolecules    = 1
    NDistMolecules= 1
    NPairs        = 6
    NCFDComp      = 2
    NProcTypes    = 3
    Syst          = system(SystNameLong, SystName, NAtoms, NMolecules, NDistMolecules, NPairs, NCFDComp, Temp.NTran, NProcTypes)


    Syst.Atom[0].Name  = 'N'
    Syst.Atom[1].Name  = 'N'
    Syst.Atom[2].Name  = 'N'
    Syst.Atom[3].Name  = 'N'

    Syst.Atom[0].Color = np.array([0, 0, 1])
    Syst.Atom[1].Color = np.array([0, 0, 1])
    Syst.Atom[2].Color = np.array([0, 0, 1])
    Syst.Atom[3].Color = np.array([0, 0, 1])

    Syst.Atom[0].Size  = 200
    Syst.Atom[1].Size  = 200
    Syst.Atom[2].Size  = 200
    Syst.Atom[3].Size  = 200

    Syst.Atom[0].Mass  = 25526.04298
    Syst.Atom[1].Mass  = 25526.04298
    Syst.Atom[2].Mass  = 25526.04298
    Syst.Atom[3].Mass  = 25526.04298


    Syst.Molecule[0].Name             = 'N2'
    Syst.Molecule[0].DissEn           = 0.0
    Syst.Molecule[0].DegeneracyFactor = 1
    Syst.Molecule[0].Mu               = 31.9988e-3
    #Syst.Molecule[0].NLevels          = 9092


    Syst.Pair[0].Name  = 'N2'
    Syst.Pair[1].Name  = 'N2'
    Syst.Pair[2].Name  = 'N2'
    Syst.Pair[3].Name  = 'N2'
    Syst.Pair[4].Name  = 'N2'
    Syst.Pair[5].Name  = 'N2'

    Syst.Pair[0].ToMol = 0
    Syst.Pair[1].ToMol = 0
    Syst.Pair[2].ToMol = 0
    Syst.Pair[3].ToMol = 0
    Syst.Pair[4].ToMol = 0
    Syst.Pair[5].ToMol = 0

    Syst.Pair[0].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[1].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[2].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[3].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[4].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[5].Color = np.array([17, 17, 17]) / 256

    Syst.CFDComp[0].Name   = 'N'
    Syst.CFDComp[1].Name   = 'N2'

    Syst.CFDComp[0].ToMol   = -1
    Syst.CFDComp[1].ToMol   =  0

    Syst.CFDComp[0].Mass    = Syst.Atom[0].Mass
    Syst.CFDComp[1].Mass    = 2.0*Syst.Atom[2].Mass

    Syst.CFDComp[0].Deg     = 1
    Syst.CFDComp[1].Deg     = 9

    Syst.CFDComp[0].Color   = np.array([ 102, 102, 102]) / 256
    Syst.CFDComp[1].Color   = np.array([   0,   0,   0]) / 256

    Syst.CFDComp[0].LineStyle = ':'
    Syst.CFDComp[1].LineStyle = '-.'

    Syst.CFDComp[0].RxLxIdx = -1
    Syst.CFDComp[1].RxLxIdx =  1

    Syst.MolToCFDComp       = np.array([1])

    Syst.ExchtoMol          = np.array([0,1])
    Syst.ExchtoAtom         = np.array([2,0])

    Syst.ColPartToComp      = 0

    return Syst
################################################################################################################


################################################################################################################
### N2 + N2 System (From NASA Ames, Dr. D. Schwenke)
def NaNbNcNd_NASA_Upload( Temp ):   

    SystNameLong  = 'NaNbNcNd_NASA'
    SystName      = 'NaNbNcNd'          

    NAtoms        = 4
    NMolecules    = 6
    NDistMolecules= 1
    NPairs        = 6
    NCFDComp      = 2
    NProcTypes    = 3
    Syst          = system(SystNameLong, SystName, NAtoms, NMolecules, NDistMolecules, NPairs, NCFDComp, Temp.NTran, NProcTypes)
    Syst.SymmFlg  = True

    Syst.Atom[0].Name  = 'N'
    Syst.Atom[1].Name  = 'N'
    Syst.Atom[2].Name  = 'N'
    Syst.Atom[3].Name  = 'N'

    Syst.Atom[0].Color = np.array([0, 0, 1])
    Syst.Atom[1].Color = np.array([0, 0, 1])
    Syst.Atom[2].Color = np.array([0, 0, 1])
    Syst.Atom[3].Color = np.array([0, 0, 1])

    Syst.Atom[0].Size  = 200
    Syst.Atom[1].Size  = 200
    Syst.Atom[2].Size  = 200
    Syst.Atom[3].Size  = 200

    Syst.Atom[0].Mass  = 25526.04298
    Syst.Atom[1].Mass  = 25526.04298
    Syst.Atom[2].Mass  = 25526.04298
    Syst.Atom[3].Mass  = 25526.04298


    Syst.Molecule[0].Name             = 'NaNb'
    Syst.Molecule[1].Name             = 'NaNc'
    Syst.Molecule[2].Name             = 'NaNd'
    Syst.Molecule[3].Name             = 'NbNc'
    Syst.Molecule[4].Name             = 'NbNd'
    Syst.Molecule[5].Name             = 'NcNd'

    Syst.Molecule[0].DissEn           = 0.0
    Syst.Molecule[0].DegeneracyFactor = 1
    Syst.Molecule[0].Mu               = 31.9988e-3
    #Syst.Molecule[0].NLevels          = 9092

    Syst.Molecule[1].DissEn           = 0.0
    Syst.Molecule[1].DegeneracyFactor = 1
    Syst.Molecule[1].Mu               = 31.9988e-3
    #Syst.Molecule[1].NLevels          = 9092

    Syst.Molecule[2].DissEn           = 0.0
    Syst.Molecule[2].DegeneracyFactor = 1
    Syst.Molecule[2].Mu               = 31.9988e-3
    #Syst.Molecule[2].NLevels          = 9092

    Syst.Molecule[3].DissEn           = 0.0
    Syst.Molecule[3].DegeneracyFactor = 1
    Syst.Molecule[3].Mu               = 31.9988e-3
    #Syst.Molecule[3].NLevels          = 9092

    Syst.Molecule[4].DissEn           = 0.0
    Syst.Molecule[4].DegeneracyFactor = 1
    Syst.Molecule[4].Mu               = 31.9988e-3
    #Syst.Molecule[4].NLevels          = 9092

    Syst.Molecule[5].DissEn           = 0.0
    Syst.Molecule[5].DegeneracyFactor = 1
    Syst.Molecule[5].Mu               = 31.9988e-3
    #Syst.Molecule[5].NLevels          = 9092

    Syst.Pair[0].Name  = 'N2'
    Syst.Pair[1].Name  = 'N2'
    Syst.Pair[2].Name  = 'N2'
    Syst.Pair[3].Name  = 'N2'
    Syst.Pair[4].Name  = 'N2'
    Syst.Pair[5].Name  = 'N2'

    Syst.Pair[0].ToMol = 0
    Syst.Pair[1].ToMol = 1
    Syst.Pair[2].ToMol = 2
    Syst.Pair[3].ToMol = 3
    Syst.Pair[4].ToMol = 4
    Syst.Pair[5].ToMol = 5

    Syst.Pair[0].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[1].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[2].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[3].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[4].Color = np.array([17, 17, 17]) / 256
    Syst.Pair[5].Color = np.array([17, 17, 17]) / 256

    Syst.CFDComp[0].Name   = 'N'
    Syst.CFDComp[1].Name   = 'N2'

    Syst.CFDComp[0].ToMol      = -1
    Syst.CFDComp[1].ToMol      =  0
    Syst.CFDComp[1].ToOppAtoms = np.array([0,0])

    Syst.CFDComp[0].Mass    = Syst.Atom[0].Mass
    Syst.CFDComp[1].Mass    = 2.0*Syst.Atom[2].Mass

    Syst.CFDComp[0].Deg     = 1
    Syst.CFDComp[1].Deg     = 9

    Syst.CFDComp[0].Color   = np.array([ 102, 102, 102]) / 256
    Syst.CFDComp[1].Color   = np.array([   0,   0,   0]) / 256

    Syst.CFDComp[0].LineStyle = ':'
    Syst.CFDComp[1].LineStyle = '-.'

    Syst.CFDComp[0].RxLxIdx   = -1
    Syst.CFDComp[1].RxLxIdx   =  1

    Syst.MolToCFDComp[0]      = 1
    Syst.MolToCFDComp[1]      = 1
    Syst.MolToCFDComp[2]      = 1
    Syst.MolToCFDComp[3]      = 1
    Syst.MolToCFDComp[4]      = 1
    Syst.MolToCFDComp[5]      = 1

    Syst.ExchtoMol          = np.array([[0,0]])

    Syst.ColPartToComp      = 0

    return Syst
################################################################################################################