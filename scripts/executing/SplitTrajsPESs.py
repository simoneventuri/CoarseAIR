from __future__ import print_function

import pandas
#import numpy
#import math
import sys

    
def read_file(PathToLabeledInput):

    #print(('    [SplitTrajsPESs.py]: Loading Trajectories from File: ' + PathToLabeledInput ))

    TrajData = pandas.read_csv(PathToLabeledInput, header=None, skiprows=1, delim_whitespace=True)
    #xDatay = pandas.read_csv(PathToData, header=0)
    TrajData = TrajData.apply(pandas.to_numeric, errors='coerce')
    #xData  = numpy.transpose(xxData.values)
    TrajData = TrajData.values

    return TrajData
    

######################################################################################################################################
### RUNNING
######################################################################################################################################
if __name__ == '__main__':

    FolderPath = str(sys.argv[1])
    #print('    [SplitTrajsPESs.py]: Folder Path = ', FolderPath )

    iPESStart  = int(sys.argv[2])
    #print('    [SplitTrajsPESs.py]: iPESStart = ', iPESStart )

    TrajData = read_file(FolderPath + '/trajectories.csv')
    TrajVec  = TrajData[:,0].astype(int)
    iPESVec  = TrajData[:,1].astype(int)
    bMaxVec  = TrajData[:,2]
    bVec     = TrajData[:,3]
    jiVec    = TrajData[:,4]
    viVec    = TrajData[:,5]
    AiVec    = TrajData[:,6]
    jfVec    = TrajData[:,7]
    vfVec    = TrajData[:,8]
    AfVec    = TrajData[:,9]

    #print('    [SplitTrajsPESs.py]: Found ', max(iPESVec), ' different PESs in the Trajectories File')
    for iPES in range(max(iPESVec)):
        #print('    [SplitTrajsPESs.py]: iPES = ', iPES+1, ' Done')

        PathToFile = FolderPath + '/trajectories.csv.' + str(iPES+iPESStart) 
        print(PathToFile)
        with open(PathToFile, 'w') as the_file:
            the_file.write('#  iTraj   iPES          bmax           b_i              j_i              v_i            arr_i              j_f              v_f            arr_f\n')

            for iTraj in range(len(iPESVec)):
                if (iPES+1 == iPESVec[iTraj]):
                    the_file.write('{:8d}{:7d}{:14.5E}{:14.5E}{:17.8E}{:17.8E}{:17.8E}{:17.8E}{:17.8E}{:17.8E}\n'.format(TrajVec[iTraj], iPESVec[iTraj]+iPESStart-1, bMaxVec[iTraj], bVec[iTraj], jiVec[iTraj], viVec[iTraj], AiVec[iTraj], jfVec[iTraj], vfVec[iTraj], AfVec[iTraj]))
