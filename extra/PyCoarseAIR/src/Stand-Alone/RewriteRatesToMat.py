import sys
import os


######################################################################################################################################
### RUNNING
######################################################################################################################################
if __name__ == '__main__':


    #WORKSPACEFldr = str(sys.argv[1])
    #WORKSPACEFldr = '/Users/sventuri/WORKSPACE'
    WORKSPACEFldr = '/nobackup/sventuri/'
    print('\n    WORKSPACEFldr Folder Path = ', WORKSPACEFldr, '\n')

    #LevelsFldr = str(sys.argv[2])
    LevelsFldr = WORKSPACEFldr + '/CoarseAIR/run_LEPS/Test/N3/N2/'
    print('    Levels Folder Path = ', LevelsFldr, '\n')

    #RatesFldr  = str(sys.argv[3])
    RatesFldr  = WORKSPACEFldr + '/CoarseAIR/run_N3_NASA/Test/N3/N2/Rates/T_10000_10000/'
    print('    Rates Folder Path = ', RatesFldr, '\n')

    #OutputFldr  = str(sys.argv[4])
    OutputFldr  = WORKSPACEFldr + '/CoarseAIR/run_N3_NASA/Postprocessing/'
    print('    Output Folder Path = ', OutputFldr, '\n')

    #BinsStart      = int(sys.argv[5])
    BinsStart      = 1
    print('    Start Bins = ', BinsStart, '\n')    

    #BinsEnd      = int(sys.argv[6])
    BinsEnd      = 9390
    print('    End Bins = ', BinsEnd, '\n')

    #T0_VecTemp = numpy.array([float(sys.argv[7])])
    T0_VecTemp = 10000
    print('    Temperature = ', T0_VecTemp, '\n')



    fNoRunPath = OutputFldr + "/NoRunLevels.out"
    fNoRun     = open(fNoRunPath, "a+")
    fNoRun.write("# No Run Levels:\n")
    fErrorPath = OutputFldr + "/ErrorLevels.out"
    fError     = open(fErrorPath, "a+")
    fError.write("# Error Levels:\n")

    for iBins in range(BinsStart-1, BinsEnd):


        PathToBinRate = RatesFldr + 'Bin' + str(iBins+1) + '.dat'
        exists = os.path.isfile(PathToBinRate)
        if exists:

            ToBinRate = open(PathToBinRate)
            Length    = len(ToBinRate.readlines())

            if (Length < 6):

                print('BIN ', iBins+1, ' DOES NOT HAVE RATES!\n')
                fError.write("%i\r\n" % (iBins+1))

        else:
            print('BIN ', iBins+1, ' DOES NOT HAVE FILE!\n')
            fNoRun.write("%i\r\n" % (iBins+1))
        