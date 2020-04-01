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
WORKSPACEFldr = '/home/aracca/WORKSPACE/'
print('\n    WORKSPACEFldr Folder Path = ', WORKSPACEFldr, '\n')

#RatesFldr  = str(sys.argv[2])
RatesFldr  = WORKSPACEFldr + '/CoarseAIR/run_O3_ALL/Test/O3/O2/Rates/'
print('    Rates Folder Path = ', RatesFldr, '\n')

#OutputFldr  = str(sys.argv[3])
OutputFldr = WORKSPACEFldr + '/CoarseAIR/TempRates/'
print('    Output Folder Path = ', OutputFldr, '\n')

#BinsStart      = int(sys.argv[4])
BinsStart  = 1
print('    Start Bins = ', BinsStart, '\n')    

#BinsEnd      = int(sys.argv[5])
BinsEnd    = 10
print('    End Bins = ', BinsEnd, '\n')

#NLevels      = int(sys.argv[6])
NLevels    = 6115
print('    NLevels = ', NLevels, '\n')    

#TVec      = float(sys.argv[7])
TVec       = numpy.array([1500, 3000, 8000, 10000, 12000, 14000])
print('    Temperature Vec = ', TVec, '\n')
NT         = TVec.shape[0]


fNoRunPath = OutputFldr + "/NoRunLevels.out"
fNoRun     = open(fNoRunPath, "a+")
fNoRun.write("# No Run Levels:\n")
fErrorPath = OutputFldr + "/ErrorLevels.out"
fError     = open(fErrorPath, "a+")
fError.write("# Error Levels:\n")


RatesMat = numpy.zeros(((BinsEnd-BinsStart+1), NLevels*3+1, TVec.shape[0]))


iT=0
for TTra in TVec:
    RateFldr = RatesFldr + '/T_' + str(TTra) + '_' + str(TTra) + '/'
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
                #print('RatesTemp = ', RatesTemp)
                for iP in range(iProc.shape[0]):
                    RatesMat[iBins-(BinsStart-1),iProc[iP]-1,iT] = float(DataMat[2][iP].replace('D', 'E'))
                    #RatesMat[iBins-(BinsStart-1),iProc[iP]-1,iT] = RatesTemp[iP]
            ToBinRate.close()
        else:
            print('T = ', TTra, '; BIN ', iBins+1, ' DOES NOT HAVE FILE!\n')
            fNoRun.write("%i\r\n" % (iBins+1))
    iT=iT+1


for iBins in range(BinsStart-1, BinsEnd):
    for kBins in range(2*NLevels+1, 3*NLevels+1):
        for iT in range(NT):
            if (RatesMat[iBins-(BinsStart-1),kBins,iT] > 0.0):
                jBins = kBins - NLevels
                RatesMat[iBins-(BinsStart-1),jBins,iT] = RatesMat[iBins-(BinsStart-1),jBins,iT] + RatesMat[iBins-(BinsStart-1),kBins,iT]
                RatesMat[iBins-(BinsStart-1),kBins,iT] = 0.0


NProcTot = 0
NDissTot = 0
NInelTot = 0
NExchTot = 0
for iBins in range(BinsStart-1, BinsEnd):    
    for kBins in range(RatesMat.shape[1]):
        if (numpy.sum(RatesMat[iBins-(BinsStart-1),kBins,:]) > 0.0):
            NProcTot = NProcTot + 1
            if (kBins == 0):
                NDissTot = NDissTot + 1
            elif (kBins <= NLevels):
                NInelTot = NInelTot + 1
            else:
                NExchTot = NExchTot + 1

print('Nb of Total        Processes: ', NProcTot)
print('Nb of Dissociation Processes: ', NDissTot)
print('Nb of Dissociation Processes: ', NInelTot)
print('Nb of Exchange     Processes: ', NExchTot)
DissRatesMat = numpy.zeros((NDissTot, NT))
InelRatesMat = numpy.zeros((NInelTot, NT))
ExchRatesMat = numpy.zeros((NExchTot, NT))
DissIMat     = numpy.zeros(NDissTot, dtype=numpy.int32)
InelIMat     = numpy.zeros(NInelTot, dtype=numpy.int32)
ExchIMat     = numpy.zeros(NExchTot, dtype=numpy.int32)
InelJMat     = numpy.zeros(NInelTot, dtype=numpy.int32)
ExchJMat     = numpy.zeros(NExchTot, dtype=numpy.int32)


NDissTot = 0
NInelTot = 0
NExchTot = 0
for iBins in range(BinsStart-1, BinsEnd):    
    for kBins in range(RatesMat.shape[1]):
        if (numpy.sum(RatesMat[iBins-(BinsStart-1),kBins,:]) > 0.0):
            NProcTot  = NProcTot + 1
            TempRates = RatesMat[iBins-(BinsStart-1),kBins,:]
            #lineRates = FormatStr.format(TempRates[0], TempRates[1], TempRates[2], TempRates[3], TempRates[4], TempRates[5])
            if (kBins == 0):
                DissIMat[NDissTot]       = iBins + 1
                DissRatesMat[NDissTot,:] = TempRates
                NDissTot                 = NDissTot + 1
                #lineProc = FormatStrProcDiss.format(iBins + 1)
                #fDiss.write(lineRates    + '\n')
                #fDissProc.write(lineProc + '\n')
            elif (kBins <= NLevels):
                InelIMat[NInelTot]       = iBins + 1
                InelJMat[NInelTot]       = kBins
                InelRatesMat[NInelTot,:] = TempRates
                NInelTot                 = NInelTot + 1
                #lineProc = FormatStrProc.format(iBins + 1, kBins)
                #fInel.write(lineRates    + '\n')
                #fInelProc.write(lineProc + '\n')
            else:
                ExchIMat[NExchTot]       = iBins + 1
                ExchJMat[NExchTot]       = kBins - NLevels
                ExchRatesMat[NExchTot,:] = TempRates
                NExchTot = NExchTot + 1
                #lineProc = FormatStrProc.format(iBins + 1, kBins)
                #fExch.write(lineRates    + '\n')
                #fExchProc.write(lineProc + '\n')

InelProc = numpy.transpose( numpy.stack( (InelIMat, InelJMat) ) )
ExchProc = numpy.transpose( numpy.stack( (ExchIMat, ExchJMat) ) )


################################################ WRITING RATES in CSV FILE ##############
HeaderStrProc     = '# iBin, jBin'
HeaderStrProcDiss = '# iBin'
HeaderStr         = '# @T=' + str(TVec[0])
for i in range(1, RatesMat.shape[2]):
    HeaderStr = HeaderStr + ', @T=' + str(TVec[i])

fDiss     = open(OutputFldr + '/DissRates.csv','w')
fDiss.write(HeaderStr + '\n')
fInel     = open(OutputFldr + '/InelRates.csv','w')
fInel.write(HeaderStr + '\n')
fExch     = open(OutputFldr + '/ExchRates.csv','w')
fExch.write(HeaderStr + '\n')
fDissProc = open(OutputFldr + '/DissProc.csv','w')
fDissProc.write(HeaderStrProcDiss + '\n')
fInelProc = open(OutputFldr + '/InelProc.csv','w')
fInelProc.write(HeaderStrProc + '\n')
fExchProc = open(OutputFldr + '/ExchProc.csv','w')
fExchProc.write(HeaderStrProc + '\n')

numpy.savetxt(fDiss,     DissRatesMat, delimiter=",", fmt='%.10e')
numpy.savetxt(fDissProc, DissIMat,     delimiter=",", fmt='%d')
numpy.savetxt(fInel,     InelRatesMat, delimiter=",", fmt='%.10e')
numpy.savetxt(fInelProc, InelProc,     delimiter=",", fmt='%d')
numpy.savetxt(fExch,     ExchRatesMat, delimiter=",", fmt='%.10e')
numpy.savetxt(fExchProc, ExchProc,     delimiter=",", fmt='%d')

fDiss.close()
fInel.close()
fExch.close()
fDissProc.close()
fInelProc.close()
fExchProc.close()
#########################################################################################


################################## WRITING RATES in BIN FORMAT for FORTRAN ##############
fInfo = open(OutputFldr + '/UnfInfo.txt','w')
fInfo.write('Nb of Dissociation Processes = {:d}\n'.format(NDissTot) )
fInfo.write('Nb of Inelastic    Processes = {:d}\n'.format(NInelTot) )
fInfo.write('Nb of Exchange     Processes = {:d}\n'.format(NExchTot) )
fInfo.write('Nb of Temperatures           = {:d}\n'.format(NT) )
fInfo.write('Temperatures:\n')
fInfo.write(numpy.array2string(TVec))
fInfo.close()

fDiss = FortranFile(OutputFldr + '/DissRates.unf','w')
fInel = FortranFile(OutputFldr + '/InelRates.unf','w')
fExch = FortranFile(OutputFldr + '/ExchRates.unf','w')

fDiss.write_record(DissIMat)
fInel.write_record(InelProc)
fExch.write_record(ExchProc)
fDiss.write_record(DissRatesMat)
fInel.write_record(InelRatesMat)
fExch.write_record(ExchRatesMat)

fDiss.close()
fInel.close()
fExch.close()
# fDiss = FortranFile(OutputFldr + '/DissRates.unf','r')
# fInel = FortranFile(OutputFldr + '/InelRates.unf','r')
# fExch = FortranFile(OutputFldr + '/ExchRates.unf','r')

# DissIMatNEW = fDiss.read_ints(numpy.int32)
# InelProc = fInel.read_ints(numpy.int32).reshape((NInelTot,2))
# InelIMatNEW = InelProc[:,0]
# IneljMatNEW = InelProc[:,1]
# ExchProc = fExch.read_ints(numpy.int32).reshape((NExchTot,2))
# ExchIMatNEW = ExchProc[:,0]
# ExchjMatNEW = ExchProc[:,1]
# DissRatesMatNEW = fDiss.read_reals(float).reshape((NDissTot,NT))
# InelRatesMatNEW = fInel.read_reals(float).reshape((NInelTot,NT))
# ExchRatesMatNEW = fExch.read_reals(float).reshape((NExchTot,NT))

# print(np.array_equal(DissIMatNEW,DissIMat))
# print(np.array_equal(InelIMatNEW,InelIMat))
# print(np.array_equal(IneljMatNEW,IneljMat))
# print(np.array_equal(ExchIMatNEW,ExchIMat))
# print(np.array_equal(ExchJMatNEW,ExchJMat))
# print(np.array_equal(DissRatesMatNEW,DissRatesMat))
# print(np.array_equal(InelRatesMatNEW,InelRatesMat))
# print(np.array_equal(ExchRatesMatNEW,ExchRatesMat))

# fDiss.close()
# fInel.close()
# fExch.close()
#########################################################################################
