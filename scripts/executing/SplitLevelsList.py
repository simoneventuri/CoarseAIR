#!/usr/bin/python

from __future__ import print_function

import pandas
#import numpy
#import math
import sys


    
def read_levels(PathToLabeledInput):

    #print(('    [SplitTrajsPESs.py]: Loading Trajectories from File: ' + PathToLabeledInput ))

    LevelsData = pandas.read_csv(PathToLabeledInput, header=None, skiprows=15, delimiter=r"\s+")
    LevelsData = LevelsData.apply(pandas.to_numeric, errors='coerce')
    LevelsData = LevelsData.values

    return LevelsData


def read_mapping(PathToLabeledInput):

    #print(('    [SplitTrajsPESs.py]: Loading Trajectories from File: ' + PathToLabeledInput ))

    MappingData = pandas.read_csv(PathToLabeledInput, header=None, skiprows=1)
    MappingData = MappingData.apply(pandas.to_numeric, errors='coerce')
    MappingData = MappingData.values

    return MappingData
    

######################################################################################################################################
### RUNNING
######################################################################################################################################
if __name__ == '__main__':

    # print('Number of arguments:', len(sys.argv), 'arguments.')
    # print('Argument List:', str(sys.argv))

    LocalDir       = sys.argv[1]
    print('      [SplitLevelsList.py]: LocalDir       = ', LocalDir)
    OutputDir      = sys.argv[2]
    print('      [SplitLevelsList.py]: OutputDir      = ', OutputDir)
    System         = sys.argv[3]
    print('      [SplitLevelsList.py]: System         = ', System)
    NInitMolecules = int(sys.argv[4])
    print('      [SplitLevelsList.py]: NInitMolecules = ', NInitMolecules)


    iTemp = 4
    if (NInitMolecules == 1):

        iTemp += 1
        Molecules_Name = sys.argv[iTemp]  
        print('      [SplitLevelsList.py]: Molecules_Name = ', Molecules_Name)

        iTemp += 1
        NBins          = int(sys.argv[iTemp])
        print('      [SplitLevelsList.py]: NBins          = ', NBins)

        iTemp += 1
        BinOI          = int(sys.argv[iTemp])
        print('      [SplitLevelsList.py]: BinOI          = ', BinOI)

        FileName    = OutputDir +  '/' + System + '/' + Molecules_Name + '/Bins_' + str(NBins)  +  '/QNsEnBin.csv'
        MappingData = read_mapping(FileName)
        Bins        = MappingData[:,5]

        FileName   = OutputDir +  '/' + System + '/' + Molecules_Name + '/levels_cut.inp'
        LevelsData = read_levels(FileName)
        NLevels    = LevelsData.shape[0]

        FileName = LocalDir + '/levels_' + Molecules_Name + '.inp'
        with open(FileName, 'w') as csvLevels:
            csvLevels.write('#================================================================================================================================================\n')
            csvLevels.write('#\n')                 
            csvLevels.write('#================================================================================================================================================\n')
            csvLevels.write('# vqn, jqn,         E[Eh],      EGam[au],      rMin[a0],      rMax[a0],      VMin[Eh],      VMax[Eh],       Tau[au],       rIn[a0],      rOut[a0]\n')
            for iLevel in range(NLevels):
                if (int(Bins[iLevel]) == BinOI):
                    csvLevels.write(' %4d,%4d,%14.7E,%14.7E,%14.7E,%14.7E,%14.7E,%14.7E,%14.7E,%14.7E,%14.7E\n' % (LevelsData[iLevel,0], LevelsData[iLevel,1],  LevelsData[iLevel,2], LevelsData[iLevel,3], LevelsData[iLevel,4], LevelsData[iLevel,5],  LevelsData[iLevel,6], LevelsData[iLevel,7], LevelsData[iLevel,8], LevelsData[iLevel,9],  LevelsData[iLevel,10]) )
        csvLevels.close()

    else:

        for iMol in range (NInitMolecules):

            iTemp += 1
            Molecules_Name = sys.argv[iTemp]  
            print('      [SplitLevelsList.py]: Molecules_Name = ', Molecules_Name)

            iTemp += 1
            NBins          = int(sys.argv[iTemp])
            print('      [SplitLevelsList.py]: NBins          = ', NBins)

            iTemp += 1
            BinOI          = int(sys.argv[iTemp])
            print('      [SplitLevelsList.py]: BinOI          = ', BinOI)

            FileName    = OutputDir +  '/' + System + '/' + Molecules_Name + '/Bins_' + str(NBins)  +  '/QNsEnBin.csv'
            MappingData = read_mapping(FileName)
            Bins        = MappingData[:,5]

            FileName   = OutputDir +  '/' + System + '/' + Molecules_Name + '/levels_cut.inp'
            LevelsData = read_levels(FileName)
            NLevels    = LevelsData.shape[0]

            FileName = LocalDir + '/levels_' + Molecules_Name + str(iMol+1) + '.inp'
            with open(FileName, 'w') as csvLevels:
                csvLevels.write('#================================================================================================================================================\n')
                csvLevels.write('#\n')                 
                csvLevels.write('#================================================================================================================================================\n')
                csvLevels.write('# vqn, jqn,         E[Eh],      EGam[au],      rMin[a0],      rMax[a0],      VMin[Eh],      VMax[Eh],       Tau[au],       rIn[a0],      rOut[a0]\n')
                for iLevel in range(NLevels):
                    if (int(Bins[iLevel]) == BinOI):
                        csvLevels.write(' %4d,%4d,%14.7E,%14.7E,%14.7E,%14.7E,%14.7E,%14.7E,%14.7E,%14.7E,%14.7E\n' % (LevelsData[iLevel,0], LevelsData[iLevel,1],  LevelsData[iLevel,2], LevelsData[iLevel,3], LevelsData[iLevel,4], LevelsData[iLevel,5],  LevelsData[iLevel,6], LevelsData[iLevel,7], LevelsData[iLevel,8], LevelsData[iLevel,9],  LevelsData[iLevel,10]) )
            csvLevels.close()

    