import os
import sys
import shutil
import fileinput
import subprocess
import numpy  as np
import pandas as pd
from matplotlib import pyplot as plt



def ComputeProbabilities_InverseRates(Iter, p2, Rates, MinRate):
    Iter  += 1
    p1     = 1. / np.maximum(Rates, MinRate)
    p2    += p1 / np.sum(p1) 
    p      = p2 / (Iter)
    return Iter, p2, p



def Plotting_Probabilities(FigFldr, LevelsData, Idx, LevelsVecTot, p):
    fig  = plt.figure(figsize=(16,8))
    ec   = plt.scatter(LevelsData['jqn'], LevelsData['eint'], c=p, cmap='Greens_r')
    plt.scatter(LevelsData['jqn'].to_numpy(int)[LevelsVecTot], LevelsData['eint'].to_numpy(np.float64)[LevelsVecTot], c='k')
    if (Idx > 0):
        plt.scatter(LevelsData['jqn'].to_numpy(int)[LevelsVecTot[-1]], LevelsData['eint'].to_numpy(np.float64)[LevelsVecTot[-1]], c='r')
    plt.xlabel('J')
    plt.ylabel('E [eV]')
    cbar = plt.colorbar(ec)
    #plt.show()
    plt.savefig(FigFldr+'/SampledLevels_'+str(Idx)+'.png', dpi=900)



def Reading_Rates(RatesFile, NLevels):
    Rates     = np.zeros(NLevels, dtype=np.float64)
    RatesData = pd.read_csv(RatesFile, header=None, skiprows=5)
    RatesData = RatesData[RatesData[0] > 2]
    RatesData = RatesData[RatesData[0] < NLevels +3]
    Rates[RatesData[0].to_numpy(int)-3] = RatesData[1].to_numpy(np.float64)
    return Rates



if __name__ == "__main__": 

    WORKSPACE_PATH = os.environ['WORKSPACE_PATH']
    CoarseAIRFldr     = WORKSPACE_PATH + '/CoarseAIR/coarseair/'
    #plt.style.use(CoarseAIRFldr+'/scripts/postprocessing/presentation.mplstyle')
    plt.style.use(CoarseAIRFldr+'/scripts/postprocessing/paper.mplstyle')



    ########################################################################################################
    ### INPUT VARIABLES

    ### Machine-Specific Input
    CoarseAIRFile  = CoarseAIRFldr  + '/scripts/launching/CoarseAIR_ActiveLearning.sh'                      
    DtbPath        = CoarseAIRFldr  + '/dtb/'        

    ### System-Specific Input
    InputFrmtFile  = CoarseAIRFldr  + '/input/O3/UMN/input/CoarseAIR_ActiveLearning.inp'                     ### NOTE: Do not forget to change quantities (e.g., no of trajectories) from
                                                                                                             ###       /CoarseAIR/CoarseAIR/scripts/PipeLine_CoarseAIR/CoarseAIR_Templates/O3_UMN/input/CoarseAIR.inp!
    System         = 'O3'
    Molecule       = 'O2'
    MinRatesVec    = [5.e-13,  1.e-13]
    TExtrVec       = [2500.0, 20000.0]

    ### Run-Specific Input
    RunFldr        = WORKSPACE_PATH + '/CoarseAIR/ActiveCoarseAIR_O3_UMN/'                                 
    NProc          = 20
    TVec           = np.array([10000.0])                                                                    
    NSamples       = 45


    for T in TVec:
        print('\n[Active-CoarseAIR] Temperature: T = ', T)



        ########################################################################################################
        ### INITIALIZING RUNS

        ### Creating Run Directory
        try:
            os.makedirs(RunFldr)
        except OSError as e:
            pass
        RunFldrT       = RunFldr + '/Test_T' + str(int(T)) + 'K/'
        try:
            os.makedirs(RunFldrT)
        except OSError as e:
            pass
        FigFldr        = RunFldr + '/Figs_T' + str(int(T)) + 'K/'
        try:
            os.makedirs(FigFldr)
        except OSError as e:
            pass


        ### Creating Input File 
        InputFldr      = RunFldr + '/Input_T' + str(int(T)) + 'K/'
        try:
            os.makedirs(InputFldr)
        except OSError as e:
            pass

        InputFile      = InputFldr + '/CoarseAIR.inp'
        shutil.copy2(InputFrmtFile, InputFile)


        ### Copying Executable
        ExecFile       = RunFldr + '/CoarseAIR_T' + str(int(T)) + 'K.sh'
        shutil.copy2(CoarseAIRFile, ExecFile)

        ### Modifying Executable
        with fileinput.FileInput(ExecFile, inplace=True) as file:
            for line in file:
                line = line.replace('***RUN_DIR***',              RunFldr)
                line = line.replace('***COARSEAIR_INPUT_DIR***',  InputFldr)
                line = line.replace('***COARSEAIR_OUTPUT_DIR***', RunFldrT)
                print(line, end='')
        print('\n[Active-CoarseAIR]   Copied and Modified CoarseAIR Executable')




        ########################################################################################################
        ### PREPROCESSING LEVELS

        ### Modifying Input File for Preprocessing Levels
        with fileinput.FileInput(InputFile, inplace=True) as file:
            for line in file:
                line = line.replace('Preprocessing Energy Levels List? = no', 'Preprocessing Energy Levels List? = yes')
                line = line.replace('***T***',        str(T))
                line = line.replace('***DtbPath***', DtbPath)
                print(line, end='')
        print('\n[Active-CoarseAIR]   Copied and Modified CoarseAIR Input File for Preprocessing Levels')


        ### Initializing Levels
        print('[Active-CoarseAIR]   Executing CoarseAIR for Preprocessing Levels')
        TempFile = RunFldrT+'/PreprocLevels.log'
        os.system('source ~/.zprofile; COARSEAIR_UPDATE; COARSEAIR_release >> '+RunFldrT+'/COARSEAIR_release.log; bash '+ExecFile+' '+str(1)+' >> '+TempFile)
        print('[Active-CoarseAIR]   Preprocessing Levels Executed. Open '+TempFile+' to Check Correct Execution')


        ### Finding Out Number of Levels
        NLevelsFile = RunFldrT + '/' + System + '/' + Molecule + '/NLevels.inp'
        NLevels     = pd.read_csv(NLevelsFile, header=None).to_numpy()[0,0]
        print('[Active-CoarseAIR]   Found ', NLevels, ' Levels')


        ### Reading Energy Levels
        LevelsFile         = RunFldrT + '/' + System + '/' + Molecule + '/levels_cut.inp'
        LevelsData         = pd.read_csv(LevelsFile, header=None, skiprows=15, delim_whitespace=True)
        LevelsData.columns = ['vqn','jqn','eint','egam','rmin','rmax','vmin','vmax','tau','ri','ro']
        LevelsData['eint'] = (LevelsData['eint'].to_numpy() - LevelsData['vmax'].min()) * 27.211399
        LevelsData['vmin'] = (LevelsData['vmin'].to_numpy() - LevelsData['vmax'].min()) * 27.211399
        LevelsData['vmax'] = (LevelsData['vmax'].to_numpy() - LevelsData['vmax'].min()) * 27.211399
        print('[Active-CoarseAIR]   Found the following Level Properies: ')
        print('[Active-CoarseAIR]   ', LevelsData.head() , '\n')    




        ########################################################################################################
        ### RUNNING FIRST BATCH OF TRAJECTORIES

        ### Modifying File for Running Trajectories
        with fileinput.FileInput(InputFile, inplace=True) as file:
            for line in file:
                line = line.replace('Preprocessing Energy Levels List? = yes', 'Preprocessing Energy Levels List? = no')
                line = line.replace('Running Trajectories? = no', 'Running Trajectories? = yes')
                line = line.replace('Postprocessing Trajectories? = no', 'Postprocessing Trajectories? = yes')
                print(line, end='')
        print('\n[Active-CoarseAIR]   Copied and Modified CoarseAIR Input File for Running Trajectories')


        ### Looking for the First Three Levels
        Level1    = LevelsData["eint"].idxmin() + 1
        Level2    = LevelsData["ro"].idxmax()   + 1
        Level3    = LevelsData["eint"].idxmax() + 1
        LevelsVec = np.array([Level1,Level2,Level3], dtype=int)
        print('[Active-CoarseAIR]   The First Batch of Trajectories will have Initial States:', LevelsVec)


        ### Writing Levels to Be Run
        fProcessesAll = open(InputFldr+'/ProcessesToRunList_ALL.inp', 'w')
        np.savetxt(fProcessesAll,                       LevelsVec[:,np.newaxis], fmt='%i')
        np.savetxt(InputFldr+'/ProcessesToRunList.inp', LevelsVec[:,np.newaxis], fmt='%i')
        print('[Active-CoarseAIR]   Saving them in '+InputFldr+'/ProcessesToRunList.inp')


        ### Running Trajectories for the First 3 Levels
        print('[Active-CoarseAIR]   Executing CoarseAIR for Running the First Batch of Trajectories')
        TempFile = RunFldrT+'/ComputeTrajs_1.log'
        os.system('source ~/.zprofile; COARSEAIR_UPDATE; COARSEAIR_release >> '+RunFldrT+'/COARSEAIR_release.log; bash '+ExecFile+' '+str(NProc)+' >> '+TempFile)
        print('[Active-CoarseAIR]   Done Executing CoarseAIR for Running the First Batch of Trajectories. Open '+TempFile+' to Check Correct Execution')




        ########################################################################################################
        ### RUNNING REMAINING TRAJECTORIES
        MinRate       = MinRatesVec[0] + (T - TExtrVec[0])*(MinRatesVec[1] - MinRatesVec[0])/(TExtrVec[1] - TExtrVec[0])

        ### Reading Rates and Computing Sampling Probabilities
        print('\n[Active-CoarseAIR]   At this Temperature, the Rates Cutting Value for Computing Probabilities is ', MinRate)
        print('[Active-CoarseAIR]   Reading Rates for Finding Initial Sampling Probabilies')
        p2            = np.zeros(NLevels, dtype=np.float64)
        Iter          = 0 

        LevelsVec     = pd.read_csv(InputFldr+'/ProcessesToRunList.inp', header=None).to_numpy(int)[:,0]
        LevelsVecTot  = []

        for iLevel in LevelsVec:
            LevelsVecTot.append(iLevel-1)

            RatesFile   = RunFldrT + '/' + System + '/Rates/T_' + str(int(T)) + '_' + str(int(T)) + '/i' + str(iLevel) + '.csv'
            Rates       = Reading_Rates(RatesFile, NLevels)
            print('[Active-CoarseAIR]   Rates for iLevel ', iLevel, ': ', Rates)
            
            Iter, p2, p = ComputeProbabilities_InverseRates(Iter, p2, Rates, MinRate)
            pMin        = np.amin(p)
            pMax        = np.amax(p)
            print('[Active-CoarseAIR]   Updated Probabilities = ', p, '; Min Probabiliy = ', pMin, '; Max Probabiliy = ', pMax)


        ### Plotting Levels and Probabilities
        Plotting_Probabilities(FigFldr, LevelsData, 0, LevelsVecTot, p)
        print('[Active-CoarseAIR]   Probabilities Plotted in ', FigFldr)


        iFig = 1
        for iSample in range(3, NSamples):
            iLevel = np.random.choice(np.arange(NLevels)+1, 1, p=p)[0]
            LevelsVecTot.append(iLevel-1)
            print('\n[Active-CoarseAIR]   Sampled Level No ', iSample+1, ':', iLevel)


            ### Writing Level to Be Run
            np.savetxt(fProcessesAll,                       [iLevel], fmt='%i')
            np.savetxt(InputFldr+'/ProcessesToRunList.inp', [iLevel], fmt='%i')
            print('[Active-CoarseAIR]   Writing New Level in '+InputFldr+'/ProcessesToRunList.inp')


            ### Running Trajectories for the First 3 Levels
            print('[Active-CoarseAIR]   Executing CoarseAIR for Running the First Batch of Trajectories')
            TempFile = RunFldrT+'/ComputeTrajs_'+str(iSample-1)+'.log'
            os.system('source ~/.zprofile; COARSEAIR_UPDATE; COARSEAIR_release >> '+RunFldrT+'/COARSEAIR_release.log; bash '+ExecFile+' '+str(NProc)+' >> '+TempFile)
            print('[Active-CoarseAIR]   Done Executing CoarseAIR for Running the First Batch of Trajectories. Open '+TempFile+' to Check Correct Execution')


            ### Reading Rates from File 
            RatesFile   = RunFldrT + '/' + System + '/Rates/T_' + str(int(T)) + '_' + str(int(T)) + '/i' + str(iLevel) + '.csv'
            Rates       = Reading_Rates(RatesFile, NLevels)
            print('[Active-CoarseAIR]   Rates for iLevel ', iLevel, ': ', Rates)


            ### Computing Probabilities 
            Iter, p2, p = ComputeProbabilities_InverseRates(Iter, p2, Rates, MinRate)
            pMin        = np.amin(p)
            pMax        = np.amax(p)
            print('[Active-CoarseAIR]   Updated Probabilities = ', p, '; Min Probabiliy = ', pMin, '; Max Probabiliy = ', pMax)


            ### Plotting Levels and Probabilities
            Plotting_Probabilities(FigFldr, LevelsData, iFig, LevelsVecTot, p)
            print('[Active-CoarseAIR]   Probabilities Plotted in ', FigFldr)
            iFig += 1


        fProcessesAll.close()
