# - Python --
##==============================================================================================================
# 
# Coarse-Grained method for Quasi-Classical Trajectories (CG-QCT) 
# 
# Copyright (C) 2018 Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
#
# Based on "VVTC" (Vectorized Variable stepsize Trajectory Code) by David Schwenke (NASA Ames Research Center). 
# 
# This program is free software; you can redistribute it and/or modify it under the terms of the 
# Version 2.1 GNU Lesser General Public License as published by the Free Software Foundation. 
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without e=ven the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU Lesser General Public License for more details. 
# 
# You should have received a copy of the GNU Lesser General Public License along with this library; 
# if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
# 
#---------------------------------------------------------------------------------------------------------------
##==============================================================================================================

import numpy

from VariablesDict   import Input, Flags, Options, PhysParams, PlotParams
from InputFile       import input_load
from Read            import read_input
from SystemVars_     import SystemVars, initialize_SystemVars, compute_SystemVars



def Run_KONIG(Input, Flags, Options, PhysParams, PlotParams)
    
    ############################################################################################################################################
    ## SPECIFING & LOADING VARIABLES
    Input, Flags, Options      = input_load(Input, Flags, Options)
    Input                      = read_input(Input)
    SystemVars                 = initialize_SystemVars(Input, SystemVars)
    SystemVars                 = compute_SystemVars(Input, SystemVars)
    ############################################################################################################################################


    ############################################################################################################################################
    if (Flags.ReUpload):
        Rates                 = numpy.zeros(SystemVars.TotDim[-1], Input.NBins[0], Input.NTint)
        #RatesSigma            = numpy.zeros(SystemVars.TotDim[-1], Input.NBins[0], Input.NTint)
        RatesMatrix           = numpy.zeros(Input.NBins[0], numpy.amax(SystemVars.yDim), 3, Input.NTint) 
        #RatesSigmaMatrix      = numpy.zeros(Input.NBins[0], numpy.amax(SystemVars.yDim), 3, Input.NTint) 
        DissRates             = numpy.zeros(Input.NBins[0], Input.NTint)
        #DissRatesSigma        = numpy.zeros(Input.NBins[0], Input.NTint)
        ProcessesRates        = numpy.zeros(Input.NBins[0], 4, Input.NTint)
        #ProcessesRatesSigma   = numpy.zeros(Input.NBins[0], 4, Input.NTint)
        ProcessesRatesOverall = []

    
    for iT in range(Input.NTint):
    
  
        if (Flag.ReUpload): 

    
            ########################################################################################################################################
            ## READING RATES
            if (Flags.PlotProcProb or Flags.PlotKs or Flags.PlotKsMoleFracs or Flags.PlotMovie or Flags.WriteNetworkCSV or Flags.ComputeMasterEq):
                RatesMatrix, DissRates, ProcessesRates = ReadRates(iT, RatesMatrix, DissRates, ProcessesRates, StartBin, FinalBin)

            
            ## READING AND COMPUTING QUANTITIES
            t, MolFracs, Ttrans, rho, P, nd, E = ReadBox(iT)


            if (PlotTemperaturesFlg):
                Tint = ReadInternalT(iT)

            
            QBins, EeV = ReadPartFunEnergy(iT)
            if (Flags.PlotProcProb or Flags.PlotKs or Flags.PlotKsMoleFracs or Flags.PlotEnergies or Flags.PlotMovie or Flags.PlotStepsFlg or Flags.WriteNetworkCSV or Flags.ComputeMasterEq):

                if (Flags.ReUpload):
                    if (Options.ReUpload == 1):
                        PopVec, PopOverQ, Pop, NSteps  = ReadPop(iT, QBins)
                        niRatio, ProcessesRatesOverall = ComputeProcessesRatesOverall(iT, Pop, ProcessesRates)

                        filenamePop = Input.PathToKONIGOutputFldr + '/T_' + str(Input.TtraVec[iT]) + '/output/' + '/Pop'
                        save(filenamePop, 'PopVec', 'PopOverQ', 'Pop', 'NSteps', 'niRatio', 'ProcessesRatesOverall', '-v7.3')

                    elif (Options.ReUpload == 2):
                        #clear PopVec PopOverQ Pop NSteps niRatio ProcessesRatesOverall
                        filenamePop = Input.PathToKONIGOutputFldr + '/T_' + str(Input.TtraVec[iT]) + '/output/' + '/Pop'

                        load(filenamePop, 'PopVec', 'PopOverQ', 'Pop', 'NSteps', 'niRatio', 'ProcessesRatesOverall')
                        #FileInelDissEigs = strcat('./FileInelDissEigs')
                        #load(FileInelDissEigs, 'A', 'D', 'DInv', 'R', 'L')


            #if (Flags.PlotProcProb or Flags.PlotEnergies or Flags.PlotMovie or Flags.PlotStepsFlg):
            NLevels, Levelvqn, Leveljqn, LevelEh, LevelEeV, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelg, LevelToBin, LevelQ, vToLevel, DeltaEintDiss, DeltaEintDepth, rIn, rOut, rMin, rMax, VMin, VMax, Tau, Egam,  dVIn, ddVIn, dVOut, ddVOut, dVdJIn, dVdJOut = ReadLevelInfo(iT, EeV)
            for iMol in range(Input.NMolecules):
                EeV[:,iMol] = EeV[:,iMol] + LevelEeV[0,iMol]
          
          
            ########################################################################################################################################
            ## COMPUTING QUANTITIES
            if (Flags.PlotProcProb or Flags.PlotEnergies or Flags.PlotMovie or Flags.PlotStepsFlg or Flags.WriteNetworkCSV or Flags.ComputeMasterEq):
                StpInstants = ComputeStepInstants(t)        

            
            ########################################################################################################################################
            ## READING / COMPUTING LEVELS-To-BINs MAPPING
            if (LevToBinFlg):

                if (LevToBinOption == 1):
                    LevToBin, qnToBin = ReadLevToBin(NLevels, Levelvqn, Leveljqn, LevelEeV, DeltaEintDiss, StpInstants, GroupNb)
                    BinEnergy    = 0.0 * EeV
                    BinEnergyMax = 0.0 * EeV

                elif (LevToBinOption == 2):
                    LevToBin, BinEnergy, BinEnergyMax, qnToBin = ComputeLevToBin(NLevels, Levelvqn, Leveljqn, LevelEeV, DeltaEintDiss, DeltaEintDepth)

            else:
                LevToBin     = 0.0 * LevelEeV + 1
                BinEnergy    = 0.0 * EeV
                BinEnergyMax = 0.0 * EeV
            ########################################################################################################################################
            

            if (Flags.PlotProcProb or Flags.PlotEnergies or Flags.PlotMovie or Flags.PlotStepsFlg):
                if (Flags.PlotMovie or Flags.PlotStepsFlg or Flags.PlotEnergies):
                    Steps = ComputeMovieSteps(iT, t)


                if (Flags.PlotEnergies):
                    eInt, eRot, eVib = ComputeEnergies(NLevels, LevelQ, Levelg, Pop, QBins, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelvqn, LevelToBin, Steps)
            
            else:
                vEeVVib = []
            ########################################################################################################################################

                
            Flags.ReUpload = False
           

        ##########################################################################################################################################
        ## PLOTTING

        # Plotting MOLE FRACTIONS
        if (Flags.PlotMolFracs):
            iFigure = PlotMolFracs(iT, iFigure, t, MolFracs)

        # Plotting TEMPERATURES
        if (PlotTemperaturesFlg):
            iFigure = PlotTemperatures(iT, iFigure, t, Tint, Ttrans)

        # Plotting QSS RATES
        if (Flags.PlotKs):
            iFigure = PlotKs(iT, iFigure, t, ProcessesRatesOverall)

        # QSS RATES PLUS MOLE FRACTION
        if (Flags.PlotKsMoleFracs):
            iFigure = PlotKsMoleFracs(iT, iFigure, t, ProcessesRatesOverall, MolFracs)

        # Plotting MOVIES
        if (Flags.PlotMovie):
            iFigure = PlotMovie(iT, iFigure, t, MolFracs, ProcessesRatesOverall, NLevels, LevelEeV, DeltaEintDiss, LevelQ, Levelg, Levelvqn, Leveljqn, LevelToBin, Pop, QBins, Steps, LevToBin, ProcessesRates)

        # Plotting STEPS
        if (Flags.PlotStepsFlg):
            iFigure = Flags.PlotSteps(iT, iFigure,t, MolFracs, ProcessesRatesOverall, NLevels, LevelEeV, LevelQ, Levelg, LevelToBin, Pop, QBins, StpInstants, Levelvqn, Leveljqn, LevToBin, DeltaEintDiss, ProcessesRates)

        # Plotting Energies
        if (Flags.PlotEnergies):
            iFigure = PlotEnergies(iT, iFigure, t, eInt, eRot, eVib, MolFracs, P, Steps)

        # Plotting Processes Probabilities
        if (Flags.PlotProcProb):
            iFigure = PlotProcProb(iT, iFigure, t, ProcessesRates, ProcessesRatesOverall, niRatio, LevelEeV, DeltaEintDiss, vEeVVib, Steps, StpInstants, Leveljqn, Levelvqn)

        # Reconstructing Q.B. Levels from updated Saha Eq.
        if (Flags.PlotNewSaha):
            iFigure = PlotNewSaha(iT, iFigure, iFigureStart, StpInstants, NLevels, Levelg, LevelEeV, DeltaEintDiss, LevToBin, BinEnergy, BinEnergyMax, t, MolFracs, nd, VMax, ProcessesRates, RatesMatrix, Tint)

        if (Flags.PlotTaus):
            iFigure = PlotTaus(iFigure) 

        if (Flags.WriteNetworkCSV):
            WriteNetworkCSV(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, Levelg, RatesMatrix, ProcessesRates, t, MolFracs, nd, Pop, StpInstants, LevelEeVVib0, LevelEeVRot)

        if (Flags.WriteRatesToParaviewCSV):
            WriteRatesToParaviewCSV(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, rOut, Levelg, RatesMatrix, ProcessesRates, t, MolFracs, nd, LevelEeVVib0, LevelEeVRot)

        if (Flags.ComputeMasterEq):
            Flags.ComputeMasterEquation(iT, NLevels, LevelEeV, LevelEeV0, Levelvqn, Leveljqn, rIn, rOut, Levelg, RatesMatrix, ProcessesRates, t, MolFracs, nd, LevelEeVVib0, LevelEeVRot)



######################################################################################################################################
### RUNNING
######################################################################################################################################
if __name__ == '__main__':

    if not os.path.exists(NNInput.PathToOutputFldrFldr):
        os.makedirs(NNInput.PathToOutputFldrFldr)

    Run_KONIG(Input, Flags, Options, PhysParams, PlotParams)