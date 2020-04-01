import numpy 
import scipy.io
import sys
import pandas
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import matplotlib
matplotlib.use('TkAgg')
from   pylab import *
from   matplotlib import cm



######################################################################################################################################
### RUNNING
######################################################################################################################################
if __name__ == '__main__':


    NBinsTot    = numpy.array([9390, 9390, 9390])
    BinsVecPlot = numpy.array([1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000])
    BinsVecPlot = BinsVecPlot - 1


    #WORKSPACEFldr = str(sys.argv[1])
    #WORKSPACEFldr = '/Users/sventuri/WORKSPACE'
    WORKSPACEFldr = '/home/venturi/WORKSPACE'
    print('\n    WORKSPACEFldr Folder Path = ', WORKSPACEFldr, '\n')

    #Run1Fldr = str(sys.argv[2])
    Run1Fldr = WORKSPACEFldr + '/CoarseAIR/run_N3_NASA'
    print('    Folder 1 = ', Run1Fldr, '\n')

    #Run2Fldr = str(sys.argv[3])
    Run2Fldr = WORKSPACEFldr + '/CoarseAIR/run_N3_NN'
    print('    Folder 2 = ', Run2Fldr, '\n')

    #LevelsFldr = str(sys.argv[4])
    LevelsFldr = Run1Fldr + '/Test/N3/N2/'
    print('    Levels Folder Path = ', LevelsFldr, '\n')

    #Rates1Fldr  = str(sys.argv[5])
    Rates1Fldr  = Run1Fldr + '/Test/N3/N2/Rates/T_10000_10000/'
    print('    Rates 1 Folder Path = ', Rates1Fldr, '\n')

    #Rates2Fldr  = str(sys.argv[6])
    Rates2Fldr  = Run2Fldr + '/Test/N3/N2/Rates/T_10000_10000/'
    print('    Rates 2 Folder Path = ', Rates2Fldr, '\n')

    #OutputFldr  = str(sys.argv[7])
    OutputFldr  = WORKSPACEFldr + '/CoarseAIR/NASA_vs_NN/'
    print('    Output Folder Path = ', OutputFldr, '\n')

    #BinsStart      = int(sys.argv[8])
    BinsStart      = 1
    print('    Start Bins = ', BinsStart, '\n')    

    #BinsEnd      = int(sys.argv[9])
    BinsEnd      = 9390
    print('    End Bins = ', BinsEnd, '\n')

    #T0_VecTemp = numpy.array([float(sys.argv[10])])
    T0_VecTemp = 10000
    print('    Temperature = ', T0_VecTemp, '\n')



    PathToBinRate = LevelsFldr + 'levels_cut.inp'
    LevelsDData   = pandas.read_csv(PathToBinRate, header=None, skiprows=15, delim_whitespace=True)
    LevelsDData   = LevelsDData.apply(pandas.to_numeric, errors='coerce')
    LevelsData    = LevelsDData.values
    

    vqn  = LevelsData[:,0].astype(int)
    jqn  = LevelsData[:,1].astype(int)
    EeV  = LevelsData[:,2] * 27.2113839712790
    EGam = LevelsData[:,3] 
    rMin = LevelsData[:,4]       
    rMax = LevelsData[:,5]          
    VMin = LevelsData[:,6] * 27.2113839712790        
    VMax = LevelsData[:,7] * 27.2113839712790               
    Tau  = LevelsData[:,8]           
    rIn  = LevelsData[:,9]           
    rOut = LevelsData[:,10]



    RateMat1        = numpy.zeros((numpy.amax(NBinsTot),numpy.amax(NBinsTot),2))
    RateMatSD1      = numpy.zeros((numpy.amax(NBinsTot),numpy.amax(NBinsTot),2)) 
    DissRates1      = numpy.zeros(NBinsTot[0])
    ProcessesRates1 = numpy.zeros((NBinsTot[0],4))

    RateMat2        = numpy.zeros((numpy.amax(NBinsTot),numpy.amax(NBinsTot),2)) 
    RateMatSD2      = numpy.zeros((numpy.amax(NBinsTot),numpy.amax(NBinsTot),2)) 
    DissRates2      = numpy.zeros(NBinsTot[0])
    ProcessesRates2 = numpy.zeros((NBinsTot[0],4)) 


    for iBins in range(BinsStart-1, BinsEnd):

        #if (iBins) in BinsVecPlot:
            #fig, ((ax1), (ax2)) = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(18,10))
            #fig = gcf()
            #plt.set_cmap('jet')
            #iBinStr = 'Bin ' + str(iBins+1)
            #fig.suptitle(iBinStr, fontsize=14)
            ##fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, sharex=True, sharey=True)
   

        print(iBins+1)

        PathToBinRate1 = Rates1Fldr + 'Bin' + str(iBins+1) + '.dat'
        exists = os.path.isfile(PathToBinRate1)
        if exists:

            ToBinRate = open(PathToBinRate1)
            Length    = len(ToBinRate.readlines())

            if (Length > 5):

                BinDData = pandas.read_csv(PathToBinRate1, header=None, skiprows=5, delim_whitespace=True)
                Length   = int(BinDData.size/4)

                for iiProc in range(Length):

                    iProc = BinDData[1][iiProc]

                    if (iProc == 1):

                        RateTemp                 = float(BinDData[2][iiProc].replace('D', 'E'))
                        DissRates1[iBins]        = RateTemp
                        ProcessesRates1[iBins,0] = RateTemp

                    elif (iProc <= NBinsTot[0] + 1):

                        iP   = 0
                        Proc = int(BinDData[1][iiProc]) - 2

                        RateTemp                       = float(BinDData[2][iiProc].replace('D', 'E'))
                        RateMat1[iBins, Proc, iP]   = RateTemp
                        RateMatSD1[iBins, Proc, iP] = float(BinDData[3][iiProc].replace('D', 'E'))
                        
                    elif (iProc <= NBinsTot[0] + NBinsTot[1] + 1):

                       iP   = 1
                       Proc = int(BinDData[1][iiProc]) - NBinsTot[0] - 2

                       RateTemp                       = float(BinDData[2][iiProc].replace('D', 'E'))
                       RateMat1[iBins, Proc, iP]   = RateTemp
                       RateMatSD1[iBins, Proc, iP] = float(BinDData[3][iiProc].replace('D', 'E'))

                    else:

                        iP   = 1
                        #iP   = 2
                        Proc = int(BinDData[1][iiProc]) - NBinsTot[0] - NBinsTot[1] - 2

                        RateTemp                       = RateMat1[iBins, Proc, iP, 0]   + float(BinDData[2][iiProc].replace('D', 'E'))
                        RateMat1[iBins, Proc, iP]   = RateTemp
                        RateMatSD1[iBins, Proc, iP] = RateMatSD1[iBins, Proc, iP, 0] + float(BinDData[3][iiProc].replace('D', 'E'))

                ProcessesRates1[iBins, 1] = numpy.sum(RateMat1[iBins, :, 0])
                ProcessesRates1[iBins, 2] = numpy.sum(RateMat1[iBins, :, 1])
                #ProcessesRates[iBins, 3] = numpy.sum(RateMat[iBins, :, 2, 0])


        PathToBinRate2 = Rates2Fldr + 'Bin' + str(iBins+1) + '.dat'
        exists = os.path.isfile(PathToBinRate2)
        if exists:

            ToBinRate2 = open(PathToBinRate2)
            Length    = len(ToBinRate2.readlines())

            if (Length > 5):

                BinDData = pandas.read_csv(PathToBinRate2, header=None, skiprows=5, delim_whitespace=True)
                Length   = int(BinDData.size/4)

                for iiProc in range(Length):

                    iProc = BinDData[1][iiProc]

                    if (iProc == 1):

                        RateTemp                 = float(BinDData[2][iiProc].replace('D', 'E'))
                        DissRates2[iBins]        = RateTemp
                        ProcessesRates2[iBins,0] = RateTemp

                    elif (iProc <= NBinsTot[0] + 1):

                        iP   = 0
                        Proc = int(BinDData[1][iiProc]) - 2

                        RateTemp                    = float(BinDData[2][iiProc].replace('D', 'E'))
                        RateMat2[iBins, Proc, iP]   = RateTemp
                        RateMatSD2[iBins, Proc, iP] = float(BinDData[3][iiProc].replace('D', 'E'))
                        
                    elif (iProc <= NBinsTot[0] + NBinsTot[1] + 1):

                       iP   = 1
                       Proc = int(BinDData[1][iiProc]) - NBinsTot[0] - 2

                       RateTemp                    = float(BinDData[2][iiProc].replace('D', 'E'))
                       RateMat2[iBins, Proc, iP]   = RateTemp
                       RateMatSD2[iBins, Proc, iP] = RateMatSD2[iBins, Proc, iP, 0] + float(BinDData[3][iiProc].replace('D', 'E'))

                    else:

                        iP   = 1
                        #iP   = 2
                        Proc = int(BinDData[1][iiProc]) - NBinsTot[0] - NBinsTot[1] - 2

                        RateTemp                    = RateMat2[iBins, Proc, iP, 0]   + float(BinDData[2][iiProc].replace('D', 'E'))
                        RateMat2[iBins, Proc, iP]   = RateTemp
                        RateMatSD2[iBins, Proc, iP] = RateMatSD2[iBins, Proc, iP, 0] + float(BinDData[3][iiProc].replace('D', 'E'))

                ProcessesRates2[iBins, 1] = numpy.sum(RateMat2[iBins, :, 0])
                ProcessesRates2[iBins, 2] = numpy.sum(RateMat2[iBins, :, 1])
                #ProcessesRates[iBins, 3] = numpy.sum(RateMat[iBins, :, 2, 0])


        # if (iBins) in BinsVecPlot:

        #     size1  = (numpy.log10(RateMat1[iBins, :, 0, 0]) + 15.0) * 10.0
        #     color1 = numpy.log10(RateMat1[iBins, :, 0, 0])
        #     size2  = (numpy.log10(RateMat2[iBins, :, 0, 0]) + 15.0) * 10.0
        #     color2 = numpy.log10(RateMat2[iBins, :, 0, 0])

        #     ax1.plot(jqn,VMin,'ko',markersize=1)
        #     ax1.plot(jqn,VMax,'ko',markersize=1)
        #     SC1 = ax1.scatter(jqn, EeV, s=size1, c=color1, alpha=0.8, vmin=-13, vmax=-8)
        #     ax1.set_xlabel(r'J', fontsize=15)
        #     ax1.set_ylabel(r'E [eV]', fontsize=15)
        #     ax1.set_title('NASA')
        #     ax1.grid(True)
        #     ax1.set_xlim([numpy.amin(jqn), numpy.amax(jqn)])
        #     ax1.set_ylim([numpy.amin(EeV), numpy.amax(EeV)])
        #     CB1 = fig.colorbar(SC1, shrink=0.8, extend='both', ax=ax1)
        #     CB1.set_label('$\log_10(K_{i,j})$')

        #     ax2.plot(jqn,VMin,'ko',markersize=1)
        #     ax2.plot(jqn,VMax,'ko',markersize=1)
        #     SC2 = ax2.scatter(jqn, EeV, s=size2, c=color2, alpha=0.8, vmin=-13, vmax=-8)
        #     ax2.set_xlabel(r'J', fontsize=15)
        #     ax2.set_ylabel(r'E [eV]', fontsize=15)
        #     ax2.set_title('NN')
        #     ax2.grid(True)
        #     ax2.set_xlim([numpy.amin(jqn), numpy.amax(jqn)])
        #     ax2.set_ylim([numpy.amin(EeV), numpy.amax(EeV)])
        #     CB2 = fig.colorbar(SC2, shrink=0.8, extend='both', ax=ax2)
        #     CB2.set_label('$\log_10(K_{i,j})$')

        #     # ax3.set_xlabel(r'J', fontsize=15)
        #     # ax3.set_ylabel(r'E [eV]', fontsize=15)
        #     # ax3.set_title('Volume and percent change')
        #     # ax3.grid(True)    
        #     # ax3.set_xlim([numpy.amin(jqn), numpy.amax(jqn)])
        #     # ax3.set_ylim([numpy.amin(EeV), numpy.amax(EeV)])          
             
        #     fig.tight_layout() 
        #     manager = plt.get_current_fig_manager()
        #     manager.full_screen_toggle()  
        #     plt.show() 
        #     FigPath = OutputFldr + '/Bin_' + str(iBins+1) + 'Inel.png'
        #     fig.savefig(FigPath, dpi=1000)
        #     #plt.close()
        
    scipy.io.savemat(OutputFldr + '/Rates_NASA_NN.mat', dict(RatesMatrixNASA=RateMat1, RatesMatrixNN=RateMat2, DissNASA=DissRates1, DissNN=DissRates2) )

    # ProcVec   = numpy.arange(NBinsTot[0])

    # colorDiss = numpy.log10(numpy.maximum(ProcessesRates[:,0], 1.e-100))
    # sizeDiss  = (15.0 + numpy.log10(numpy.maximum(ProcessesRates[:,0], 1.e-100))) * 5.0

    # colorExch = numpy.log10(numpy.maximum(ProcessesRates[:,2], 1.e-100))
    # sizeExch  = (15.0 + numpy.log10(numpy.maximum(ProcessesRates[:,2], 1.e-100))) * 5.0

    # fig, ((ax1), (ax2)) = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(18,10))
    # fig = gcf()
    # plt.set_cmap('jet')
    # iBinStr = '$N_2$'
    # fig.suptitle(iBinStr, fontsize=14)

    # ax1.plot(jqn,VMin,'ko',markersize=1)
    # ax1.plot(jqn,VMax,'ko',markersize=1)
    # SC1 = ax1.scatter(jqn[ProcVec], EeV[ProcVec], s=sizeDiss, c=colorDiss, alpha=0.8, vmin=-13, vmax=-8)
    # ax1.set_xlabel(r'J', fontsize=15)
    # ax1.set_ylabel(r'E [eV]', fontsize=15)
    # ax1.set_title('Dissociation Rates')
    # ax1.grid(True)
    # ax1.set_xlim([numpy.amin(jqn), numpy.amax(jqn)])
    # ax1.set_ylim([numpy.amin(EeV), numpy.amax(EeV)])
    # CB1 = fig.colorbar(SC1, shrink=0.8, extend='both', ax=ax1)
    # CB1.set_label('$\log_10(K_{i,j})$')

    # ax2.plot(jqn,VMin,'ko',markersize=1)
    # ax2.plot(jqn,VMax,'ko',markersize=1)
    # SC2 = ax2.scatter(jqn[ProcVec], EeV[ProcVec], s=sizeExch, c=colorExch, alpha=0.8, vmin=-13, vmax=-8)
    # ax2.set_xlabel(r'J', fontsize=15)
    # ax2.set_ylabel(r'E [eV]', fontsize=15)
    # ax2.set_title('Exchange Rates')
    # ax2.grid(True)
    # ax2.set_xlim([numpy.amin(jqn), numpy.amax(jqn)])
    # ax2.set_ylim([numpy.amin(EeV), numpy.amax(EeV)])
    # CB2 = fig.colorbar(SC2, shrink=0.8, extend='both', ax=ax2)
    # CB2.set_label('$\log_10(K_{i,j})$')

    # fig.tight_layout() 
    # manager = plt.get_current_fig_manager()
    # manager.full_screen_toggle()  
    # #plt.show()   
    # FigPath = OutputFldr + '/Processes.png'
    # fig.savefig(FigPath, dpi=1000)
    # #plt.close()
