import sys
import os
import numpy
import pandas
import scipy
from scipy.io import FortranFile
import matplotlib.pyplot as plt
import pickle


######################################################################################################################################
### RUNNING
######################################################################################################################################

#WORKSPACEFldr = str(sys.argv[1])
WORKSPACEFldr = '/home/aracca/WORKSPACE/'
print('\n    WORKSPACEFldr Folder Path = ', WORKSPACEFldr, '\n')

RatesFldr = [' ', ' ']
RatesFldr = WORKSPACEFldr + '/CoarseAIR/run_O3_ALL/Test/O3/O2/Rates/'
print('    Rates Folder Path = ', RatesFldr, '\n')

QNFile = WORKSPACEFldr + '/CoarseAIR/run_O3_ALL/Test/O3/O2/O2_6115/qnsFirst.dat'
print('    Path to Level-to-QNs file = ', QNFile, '\n')

#OutputFldr  = str(sys.argv[3])
OutputFldr = WORKSPACEFldr + '/CoarseAIR/TempRates/'
print('    Output Folder Path = ', OutputFldr, '\n')

#BinsStart      = int(sys.argv[4])
BinsStart  = 1
print('    Start Bins = ', BinsStart, '\n')    

#BinsEnd      = int(sys.argv[5])
BinsEnd    = 6115
print('    End Bins = ', BinsEnd, '\n')

#NLevels      = int(sys.argv[6])
NLevels    = 6115
print('    NLevels = ', NLevels, '\n')    

#TVec      = float(sys.argv[7])
TVec       = numpy.array([10000])
print('    Temperature Vec = ', TVec, '\n')
NT         = TVec.shape[0]


fNoRunPath = OutputFldr + "/NoRunLevels.out"
fNoRun     = open(fNoRunPath, "a+")
fNoRun.write("# No Run Levels:\n")
fErrorPath = OutputFldr + "/ErrorLevels.out"
fError     = open(fErrorPath, "a+")
fError.write("# Error Levels:\n")


RatesMat    = numpy.zeros(((BinsEnd-BinsStart+1), NLevels*3+1, TVec.shape[0]))

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
                #print('RatesTemp = ', DataMat[2].apply(pandas.to_numeric, errors='coerce') )
                for iP in range(iProc.shape[0]):
                    RatesMat[iBins-(BinsStart-1),iProc[iP]-1,iT] = float(DataMat[2][iP].replace('D', 'E')) 
                    #RatesMat[iBins-(BinsStart-1),iProc[iP]-1,iT] = float(DataMat[2][iP])
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


f = open(RatesFldr + 'RatesMat.pckl', 'wb')
pickle.dump(RatesMat, f)
f.close()


# f = open(RatesFldr + 'RatesMat.pckl', 'rb')
# RatesMat = pickle.load(f)
# f.close()


# for TTra in TVec:
#     RateFldr = OutputFldr + '/TEMP_T_' + str(TTra) + '_' + str(TTra) + '/'
#     if not os.path.exists(RateFldr):
#         os.makedirs(RateFldr)
#     for iBins in range(0, 2):
#         PathToBinRate = RateFldr  + '/Bin' + str(iBins+1) + '.dat'
#         fBinRates = open(PathToBinRate, "w")
#         fBinRates.write("# \n")
#         fBinRates.write("# \n")
#         fBinRates.write("# \n")
#         fBinRates.write("# \n")
#         fBinRates.write("# \n")
#         for jBins in range(NLevels*3+1):
#             if (numpy.sum(RatesMat[iBins-(BinsStart-1),jBins,0]) > 0.0):
#                 fBinRates.write("                  - %20d    %.10e    %.10e\n" % (jBins+1, RatesMat[iBins-(BinsStart-1),jBins,0], 0.0))
#         fBinRates.close()



DataMat   = pandas.read_csv(QNFile, header=None, skiprows=1, sep='\s+')
DataMat   = DataMat.apply(pandas.to_numeric, errors='coerce')
LevelToQN = DataMat.values[:,0:2]
LevelToQN = LevelToQN.astype(int)


DeltaJMax = 3
DeltaVMax = 2
Areaqn    = (2*DeltaJMax+1) * (2*DeltaVMax+1)
QNVec     = numpy.zeros((NLevels,Areaqn))
QNTot     = numpy.zeros(NLevels)
for jBins in range(NLevels):
    jLevel  = jBins
    Jjqn    = LevelToQN[jLevel,1]
    Jvqn    = LevelToQN[jLevel,0]
    ii      = 0 
    for kBins in range(NLevels):
        kLevel = kBins
        Kjqn   = LevelToQN[kLevel,1]
        Kvqn   = LevelToQN[kLevel,0]
        if ( (abs(Kjqn - Jjqn) <= DeltaJMax) and (abs(Kvqn - Jvqn) <= DeltaVMax) ):
            QNVec[jBins,ii] = kLevel
            ii              = ii + 1      
    QNTot[jBins] = ii
    print("iLevel = ", jLevel, "; vqn = ", Jvqn, "; jqn = ", Jjqn, "; QNTot[jBins] = ", QNTot[jBins], "; QNVec[jBins,:] = ", QNVec[jBins,:])
QNTot = QNTot.astype(int)
QNVec = QNVec.astype(int)

f = open(RatesFldr + 'QNVec.pckl', 'wb')
pickle.dump(QNVec, f)
f.close()

f = open(RatesFldr + 'QNTot.pckl', 'wb')
pickle.dump(QNTot, f)
f.close()

# f = open(RatesFldr + 'QNVec.pckl', 'rb')
# QNVec = pickle.load(f)
# f.close()

# f = open(RatesFldr + 'QNTot.pckl', 'rb')
# QNTot = pickle.load(f)
# f.close()


RatesMatNew = numpy.zeros(((BinsEnd-BinsStart+1), NLevels*3+1, TVec.shape[0]))

#for iBins in range(BinsStart-1, BinsEnd):
# for iBins in range(2):
#     print('iBins=',iBins)
#     for jBins in range(NLevels):
#         kBins                         = jBins
#         RatesMatNew[iBins,jBins+1,iT] = RatesMatNew[iBins,jBins+1,iT] + RatesMat[iBins-(BinsStart-1),kBins+1,iT]
#         RatesMatNew[iBins,jBins+1,iT] = RatesMatNew[iBins,jBins+1,iT] / 1
#     for jBins in range(NLevels, 2*NLevels):
#         kBins                         = jBins
#         RatesMatNew[iBins,jBins+1,iT] = RatesMatNew[iBins,jBins+1,iT] + RatesMat[iBins-(BinsStart-1),kBins+1,iT]
#         RatesMatNew[iBins,jBins+1,iT] = RatesMatNew[iBins,jBins+1,iT] / 1

for iBins in range(BinsStart-1, BinsEnd):
#for iBins in range(2):
    RatesMatNew[iBins,0,iT] = RatesMatNew[iBins,0,iT]
    print('iBins=',iBins)
    for jBins in range(NLevels):
        for ii in range(QNTot[jBins]):
            kBins                         = QNVec[jBins,ii]
            RatesMatNew[iBins,jBins+1,iT] = RatesMatNew[iBins,jBins+1,iT] + RatesMat[iBins-(BinsStart-1),kBins+1,iT]
        RatesMatNew[iBins,jBins+1,iT] = RatesMatNew[iBins,jBins+1,iT] / QNTot[jBins]           
    for jBins in range(NLevels, 2*NLevels):
        for ii in range(QNTot[jBins-NLevels]):
            kBins                         = QNVec[jBins-NLevels,ii]
            RatesMatNew[iBins,jBins+1,iT] = RatesMatNew[iBins,jBins+1,iT] + RatesMat[iBins-(BinsStart-1),kBins+1+NLevels,iT]
        RatesMatNew[iBins,jBins+1,iT] = RatesMatNew[iBins,jBins+1,iT] / QNTot[jBins-NLevels]


plt.figure()
ax = plt.gca()
ax.scatter(numpy.linspace(0, RatesMat.shape[1]-1, RatesMat.shape[1], dtype=int), RatesMat[0,:,0],    c='blue')
ax.scatter(numpy.linspace(0, RatesMat.shape[1]-1, RatesMat.shape[1], dtype=int), RatesMatNew[0,:,0], c='red')
ax.set_yscale('log')
ax.set_ylim([1.e-20,1.e-8])
#plt.semilogy(numpy.linspace(0, RatesMat.shape[1]-1, RatesMat.shape[1], dtype=int), RatesMat[0,:,0])
#plt.semilogy(numpy.linspace(0, RatesMat.shape[1]-1, RatesMat.shape[1], dtype=int), RatesMatNew[0,:,0])
plt.savefig("Bin1.png")


f = open(RatesFldr + 'RatesMatNew.pckl', 'wb')
pickle.dump(RatesMatNew, f)
f.close()


# f = open(RatesFldr + 'RatesMatNew.pckl', 'rb')
# RatesMat = pickle.load(f)
# f.close()


for TTra in TVec:
    RateFldr = OutputFldr + '/AVERAGED_T_' + str(TTra) + '_' + str(TTra) + '/'
    if not os.path.exists(RateFldr):
        os.makedirs(RateFldr)
    for iBins in range(BinsStart-1, BinsEnd):
        PathToBinRate = RateFldr  + '/Bin' + str(iBins+1) + '.dat'
        fBinRates = open(PathToBinRate, "w")
        fBinRates.write("# \n")
        fBinRates.write("# \n")
        fBinRates.write("# \n")
        fBinRates.write("# \n")
        fBinRates.write("# \n")
        for jBins in range(NLevels*3+1):
            if (numpy.sum(RatesMat[iBins-(BinsStart-1),jBins,0]) > 0.0):
                fBinRates.write("                  - %20d    %.10e    %.10e\n" % (jBins+1, RatesMatNew[iBins-(BinsStart-1),jBins,0], 0.0))
        fBinRates.close()
