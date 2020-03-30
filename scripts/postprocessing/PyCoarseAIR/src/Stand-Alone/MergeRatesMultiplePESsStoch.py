import sys
import os
import numpy
import pandas
import scipy
from scipy.io import FortranFile

######################################################################################################################################
### RUNNING
######################################################################################################################################

#WORKSPACEFldr = str(sys.argv[1])
WORKSPACEFldr = '/home/venturi/WORKSPACE/'
print('\n    WORKSPACEFldr Folder Path = ', WORKSPACEFldr, '\n')

RatesFldr = [' ', ' ']
RatesFldr[0] = WORKSPACEFldr + '/CoarseAIR/run_O3_8PESs_CG/Test/O3/O2/Rates/'
print('    PES 1, Rates Folder Path = ', RatesFldr[0], '\n')
RatesFldr[1]  = WORKSPACEFldr + '/CoarseAIR/run_O3_PES9_BNN/Test/O3/O2/Rates/'
print('    PES 2, Rates Folder Path = ', RatesFldr[1], '\n')

PESsWeights = numpy.array([22/27, 5/27])
print('    PESs Weights = ', PESsWeights, '\n')

#OutputFldr  = str(sys.argv[3])
OutputFldr = WORKSPACEFldr + '/CoarseAIR/TempRates/'
print('    Output Folder Path = ', OutputFldr, '\n')

#BinsStart      = int(sys.argv[4])
BinsStart  = 1
print('    Start Bins = ', BinsStart, '\n')    

#BinsEnd      = int(sys.argv[5])
BinsEnd    = 200
print('    End Bins = ', BinsEnd, '\n')

#NLevels      = int(sys.argv[6])
NLevels    = 200
print('    NLevels = ', NLevels, '\n')    

#TVec      = float(sys.argv[7])
TVec       = numpy.array([5000])
print('    Temperature Vec = ', TVec, '\n')
NT         = TVec.shape[0]

NPEStoch   = 50
print('    Nb of BNN Samples = ', NPEStoch, '\n')


fNoRunPath = OutputFldr + "/NoRunLevels.out"
fNoRun     = open(fNoRunPath, "a+")
fNoRun.write("# No Run Levels:\n")
fErrorPath = OutputFldr + "/ErrorLevels.out"
fError     = open(fErrorPath, "a+")
fError.write("# Error Levels:\n")


for jPES in range(NPEStoch):
    RatesMat = numpy.zeros(((BinsEnd-BinsStart+1), NLevels*3+1, TVec.shape[0]))
    iPES=0
    iT=0
    for TTra in TVec:
        RateFldr = RatesFldr[iPES] + '/T_' + str(TTra) + '_' + str(TTra) + '/'
        for iBins in range(BinsStart-1, BinsEnd):
            PathToBinRate = RateFldr  + 'Bin' + str(iBins+1) + '.dat'
            exists = os.path.isfile(PathToBinRate)
            if exists:
                ToBinRate = open(PathToBinRate)
                Length    = len(ToBinRate.readlines())
                if (Length < 6):
                    print('T = ', TTra, '; BIN ', iBins+1, ' DOES NOT HAVE RATES!\n')
                    fError.write("%i\r\n" % (iBins+1))
                else:
                    DataMat   = pandas.read_csv(PathToBinRate, header=None, skiprows=5, sep='\s+')
                    iProc     = DataMat[1].apply(pandas.to_numeric, errors='coerce')
                    iProc     = iProc.values
                    print(PathToBinRate)
                    #print('iProc     = ', iProc)
                    #print('RatesTemp = ', DataMat[2].apply(pandas.to_numeric, errors='coerce') )
                    for iP in range(iProc.shape[0]):
                        #RatesMat[iBins-(BinsStart-1),iProc[iP]-1,iT] = RatesMat[iBins-(BinsStart-1),iProc[iP]-1,iT] + float(DataMat[2][iP].replace('D', 'E')) * PESsWeights[iPES]
                        RatesMat[iBins-(BinsStart-1),iProc[iP]-1,iT] = RatesMat[iBins-(BinsStart-1),iProc[iP]-1,iT] + float(DataMat[2][iP]) * PESsWeights[iPES]
                ToBinRate.close()
            else:
                print('T = ', TTra, '; BIN ', iBins+1, ' DOES NOT HAVE FILE!\n')
                fNoRun.write("%i\r\n" % (iBins+1))
        iT=iT+1
        iPES=1
        iT=0
        for TTra in TVec:
            RateFldr = RatesFldr[iPES] + '/T_' + str(TTra) + '_' + str(TTra) + '/'
            for iBins in range(BinsStart-1, BinsEnd):
                PathToBinRate = RateFldr  + 'Bin' + str(iBins+1) + '.dat.' + str(jPES+1)
                exists = os.path.isfile(PathToBinRate)
                if exists:
                    ToBinRate = open(PathToBinRate)
                    Length    = len(ToBinRate.readlines())
                    if (Length < 6):
                        print('T = ', TTra, '; BIN ', iBins+1, ' DOES NOT HAVE RATES!\n')
                        fError.write("%i\r\n" % (iBins+1))
                    else:
                        DataMat   = pandas.read_csv(PathToBinRate, header=None, skiprows=5, sep='\s+')
                        iProc     = DataMat[1].apply(pandas.to_numeric, errors='coerce')
                        iProc     = iProc.values
                        print(PathToBinRate)
                        #print('iProc     = ', iProc)
                        #print('RatesTemp = ', DataMat[2].apply(pandas.to_numeric, errors='coerce') )
                        for iP in range(iProc.shape[0]):
                            #RatesMat[iBins-(BinsStart-1),iProc[iP]-1,iT] = RatesMat[iBins-(BinsStart-1),iProc[iP]-1,iT] + float(DataMat[2][iP].replace('D', 'E')) * PESsWeights[iPES]
                            RatesMat[iBins-(BinsStart-1),iProc[iP]-1,iT] = RatesMat[iBins-(BinsStart-1),iProc[iP]-1,iT] + float(DataMat[2][iP]) * PESsWeights[iPES]
                    ToBinRate.close()
                else:
                    print('T = ', TTra, '; BIN ', iBins+1, ' DOES NOT HAVE FILE!\n')
                    fNoRun.write("%i\r\n" % (iBins+1))
            iT=iT+1
    #
    #
    # for iBins in range(BinsStart-1, BinsEnd):
    #     for kBins in range(2*NLevels+1, 3*NLevels+1):
    #         for iT in range(NT):
    #             if (RatesMat[iBins-(BinsStart-1),kBins,iT] > 0.0):
    #                 jBins = kBins - NLevels
    #                 RatesMat[iBins-(BinsStart-1),jBins,iT] = RatesMat[iBins-(BinsStart-1),jBins,iT] + RatesMat[iBins-(BinsStart-1),kBins,iT]
    #                 RatesMat[iBins-(BinsStart-1),kBins,iT] = 0.0
    #
    #
    for TTra in TVec:
        RateFldr = OutputFldr + '/T_' + str(TTra) + '_' + str(TTra) + '/'
        if not os.path.exists(RateFldr):
            os.makedirs(RateFldr)
        for iBins in range(BinsStart-1, BinsEnd):
            PathToBinRate = RateFldr  + '/Bin' + str(iBins+1) + '.dat.' + str(jPES+1)
            fBinRates = open(PathToBinRate, "w")
            fBinRates.write("# \n")
            fBinRates.write("# \n")
            fBinRates.write("# \n")
            fBinRates.write("# \n")
            fBinRates.write("# \n")
            for jBins in range(NLevels*3+1):
                if (numpy.sum(RatesMat[iBins-(BinsStart-1),jBins,0]) > 0.0):
                    fBinRates.write("                  - %20d    %.10e    %.10e\n" % (jBins+1, RatesMat[iBins-(BinsStart-1),jBins,0], 0.0))
            fBinRates.close()
















