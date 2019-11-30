% -- MATLAB --
%%==============================================================================================================
% 
% Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
% 
% Copyright (C) 2018 Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
%
% Based on "VVTC" (Vectorized Variable stepsize Trajectory Code) by David Schwenke (NASA Ames Research Center). 
% 
% This program is free software; you can redistribute it and/or modify it under the terms of the 
% Version 2.1 GNU Lesser General Public License as published by the Free Software Foundation. 
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without e=ven the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU Lesser General Public License for more details. 
% 
% You should have received a copy of the GNU Lesser General Public License along with this library; 
% if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
% 
%---------------------------------------------------------------------------------------------------------------
%%==============================================================================================================

close all
clc
global ReUpload
ReUpload = 1
if ReUpload == 1
  clear all
  global ReUpload
  ReUpload = 1
end

global PathToOutput System T0_Vec NBins KinMthd NTtra NTint TInt_Vec NMolecules StartBin FinalBin NBinnedMol BinnedMolName NPESs

global NAtoms AtomsName MoleculesName DegeneracyFactor ColPartToComp BinnedMolToComp NComp CompNames CompColor AtomColor AtomSize AllMoleculesName ...
       PairColor AtomMass ComponentMass ColorVec ComponentDeg MoleculeMu MoleculedDissEn
     
global Plnck UKb Ue KeV AvN AMUToKg EhToeV DSWtoKg ATMToPa

global SystemPath RatesPath DatabasePath RunKONIGDir OutputPath OutputPathSaveFigs FigDirPath linS linST AxisFontSz AxisFontNm AxisLabelSz AxisLabelNm ...
       LegendFontSz LegendFontNm XLimPlot YLimPlot 

global PlotMolFracsFlg PlotTemperaturesFlg PlotKsFlg PlotKsMoleFracsFlg PlotMovieFlg PlotStepsFlg PlotEnergiesFlg PlotQSSDistrFlg PlotProcProbFlg PlotNewSahaFlg ...
       PlotTausFlg LevToBinFlg WriteNetworkCSVFlg WriteRatesToParaviewCSVFlg ComputeMasterEquationFlg PlotMolFracsStochFlg

global NStepsMovie TimeMinMovie TimeMaxMovie NSubplots NInstantsPerSubplot tInstants MinEvPlot MaxEvPlot StepsOverlappingSteps vqnColor SubplotsFlg ...
       QNSpace RxLxIdx ReadAllRatesFlg SaveRatesFlg CompOI PlotPairUser BinnedFlg MovieFlg WriteRatesFlag WriteSrcTermsFlag WriteInelasticFlag WriteAllFlag ...
       WriteDissFlag

     
global xDim yDim TotDim PlotPair ProcToLevIP TempChar TempChar2 Pair_Name Pair_To_Molecule Pair_to_Atoms iPInternal iPExternal SaveFigs

global RCVec BCVec GCVec KCVec OCVec PCVec WCVec JCVec YCVec CCVec MCVec


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SPECIFING & LOADING INPUT 
PathToOutput = '../Test/'
ReadInput()
System  = 'O3'
T0_Vec  = [10000];
NBins   = [ 6115];
KinMthd = ['STS'];
UpdateAllInput() 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHOOSING POSTPROCESSING FUNCTIONS
PlotMolFracsFlg            = 1
PlotTemperaturesFlg        = 0
PlotKsFlg                  = 1
PlotKsMoleFracsFlg         = 0
PlotMovieFlg               = 0
PlotStepsFlg               = 1
PlotEnergiesFlg            = 1
PlotQSSDistrFlg            = 0
PlotProcProbFlg            = 0
PlotNewSahaFlg             = 0
PlotTausFlg                = 0
LevToBinFlg                = 0
WriteNetworkCSVFlg         = 0
WriteRatesToParaviewCSVFlg = 0
ComputeMasterEquationFlg   = 0
PlotMolFracsStochFlg       = 1
  NStepsMovie           = 150;
  TimeMinMovie(1:3)     = 1.d-14;
  TimeMaxMovie(1)       = 1.e-5; %2.3E-07;
  TimeMaxMovie(2)       = 1.e-5; %2.3E-07;
  TimeMaxMovie(3)       = 1.e-5; %2.3E-07;
  StartBin              = 1
  FinalBin              = 6115
  NSubplots             = 3
  NInstantsPerSubplot   = 3
  tInstants             = [1.e-11, 1.e-10, 1.e-9, 1.e-8, 1.e-7, 1.e-6, 1.e-5];  
  MinEvPlot             = [-20.d0, -1.d2];
  MaxEvPlot             = [1.d2,  1.d2];
  StepsOverlappingSteps = 1
  vqnColor              = 0
  SubplotsFlg           = 0
  QNSpace               = 0
  ReadAllRatesFlg       = -1
  SaveRatesFlg          = 1
  XLimPlot              = [1.d-12, 1.d-1]
  YLimPlot              = [ 0.d0,  1.d0]
  CompOI                = 3
  PlotPairUser          = [0, 1, 1]
  BinnedFlg             = 5
  MovieFlg              = 0
  ColPartToComp         = 1
  iFigureStart          = 2001  
  WriteRatesFlag        = 1
  WriteSrcTermsFlag     = 0
  WriteDissFlag         = 0
  WriteInelasticFlag    = 1
  WriteAllFlag          = 0
  NPESs                 = 50  
  DissCorrectionFactor  = 1.0 %16.0/3.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UPDATE PARAMETERS & PATHS
[iFigure]   = UpdateParameters();

RatesPath   = strcat(PathToOutput,System(1,:),'/',MoleculesName(1,:),'/Rates/');                                                             % For Levels / Bins Rates
SystemPath  = strcat(PathToOutput,System(1,:));   

iFigure    = 3001;
SaveFigs   = 0
FigDirPath = strcat('./CGQCT-',System,'-Figures/');
if SaveFigs > 0
  [status,msg,msgID] = mkdir(FigDirPath)
end   

RunKONIGDir = 'RunKonig_normalCond';
OutputPath  = strcat(PathToOutput,RunKONIGDir,'/output/'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTING SYSTEM MAIN QUANTITIES
ComputeSystemQuantities() 
if ReUpload > 0
  Rates               = zeros(TotDim(end),NBins(1),NTint);
%   RatesSigma          = zeros(TotDim(end),NBins(1),NTint);
  RatesMatrix         = zeros(NBins(1),max(yDim),3,NTint); 
%   RatesSigmaMatrix    = zeros(NBins(1),max(yDim),3,NTint); 
  DissRates           = zeros(NBins(1),NTint);
%   DissRatesSigma      = zeros(NBins(1),NTint);
  ProcessesRates      = zeros(NBins(1),4,NTint); 
%   ProcessesRatesSigma = zeros(NBins(1),4,NTint); 
  ProcessesRatesOverall = [];
end

    
for iT = 1:length(T0_Vec)
    
  
  if ReUpload > 0 

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% READING RATES
    if PlotProcProbFlg == 1 || PlotKsFlg == 1 || PlotKsMoleFracsFlg == 1 || PlotMovieFlg == 1 || WriteNetworkCSVFlg == 1 || ComputeMasterEquationFlg == 1
      [RatesMatrix, DissRates, ProcessesRates] = ReadRates(iT, RatesMatrix, DissRates, ProcessesRates, StartBin, FinalBin, DissCorrectionFactor);
    end
 
     if PlotMolFracsStochFlg == 1 
%       NPESs = 100
%       if ReadAllRatesFlg == -1
%         filenamePES = strcat('./StochPES')
%         load(filenamePES, 'StochPESR', 'StochPESAng', 'StochPESEeV')
%       else
%         [StochPESR, StochPESAng, StochPESEeV] = ReadStochPES();
%       end 
%       
%       RatesMatrixStoch         = zeros(NBins(1),max(yDim),3,NPESs); 
%       %RatesSigmaMatrixStoch    = zeros(NBins(1),max(yDim),3,NPESs); 
%       DissRatesStoch           = zeros(NBins(1),NPESs);
%       %DissRatesSigmaStoch      = zeros(NBins(1),NPESs);
%       ProcessesRatesStoch      = zeros(NBins(1),4,NPESs); 
%       %ProcessesRatesSigmaStoch = zeros(NBins(1),4,NPESs); 
%       [RatesMatrixStoch, DissRatesStoch, ProcessesRatesStoch] = ReadRatesStoch(iT, RatesMatrixStoch, DissRatesStoch, ProcessesRatesStoch, StartBin, FinalBin);
%       
      [tStoch, MolFracsStoch, TtransStoch, rhoStoch, PStoch, ndStoch, EStoch] = ReadBoxStoch(iT); 
    else
      [t, MolFracs, Ttrans, rho, P, nd, E] = ReadBox(iT); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    %% READING AND COMPUTING QUANTITIES
    [t, MolFracs, Ttrans, rho, P, nd, E] = ReadBox(iT); 


    if PlotTemperaturesFlg == 1
      [Tint] = ReadInternalT(iT);
    end

    
    [QBins, EeV] = ReadPartFunEnergy(iT); 
    if PlotProcProbFlg == 1 || PlotKsFlg == 1 || PlotKsMoleFracsFlg == 1 || PlotEnergiesFlg == 1 || PlotMovieFlg == 1 || PlotStepsFlg == 1 || WriteNetworkCSVFlg == 1 || ComputeMasterEquationFlg == 1
      if ReUpload == 1
        [PopVec, PopOverQ, Pop, NSteps] = ReadPop(iT, QBins);
        [niRatio, ProcessesRatesOverall] = ComputeProcessesRatesOverall(iT, Pop, ProcessesRates);
        filenamePop = strcat(OutputPath,'/T_',num2str(T0_Vec(iT)),'/output/','/Pop')
        save(filenamePop, 'PopVec', 'PopOverQ', 'Pop', 'NSteps', 'niRatio', 'ProcessesRatesOverall', '-v7.3')
      elseif ReUpload == 2
        clear PopVec PopOverQ Pop NSteps niRatio ProcessesRatesOverall
        filenamePop = strcat(OutputPath,'/T_',num2str(T0_Vec(iT)),'/output/','/Pop')
        load(filenamePop, 'PopVec', 'PopOverQ', 'Pop', 'NSteps', 'niRatio', 'ProcessesRatesOverall')
        %FileInelDissEigs = strcat('./FileInelDissEigs')
        %load(FileInelDissEigs, 'A', 'D', 'DInv', 'R', 'L')
      end
    end
    

    %if PlotProcProbFlg == 1 || PlotEnergiesFlg == 1 || PlotMovieFlg == 1 || PlotStepsFlg == 1
      [NLevels, Levelvqn, Leveljqn, LevelEh, LevelEeV, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelg, LevelToBin, LevelQ, vToLevel, DeltaEintDiss, DeltaEintDepth, rIn, rOut, rMin, rMax, VMin, VMax, Tau, Egam,  dVIn, ddVIn, dVOut, ddVOut, dVdJIn, dVdJOut] = ReadLevelInfo(iT, EeV);
      for iMol = 1:NMolecules
        EeV(:,iMol) = EeV(:,iMol) + LevelEeV(1,iMol); %- MoleculedDissEn(iMol);
      end
    %end
  
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% COMPUTING QUANTITIES
    if PlotProcProbFlg == 1 || PlotEnergiesFlg == 1 || PlotQSSDistrFlg == 1 || PlotMovieFlg == 1 || PlotStepsFlg == 1 || WriteNetworkCSVFlg == 1 || ComputeMasterEquationFlg == 1
      [StpInstants] = ComputeStepInstants(t);        
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% READING / COMPUTING LEVELS-To-BINs MAPPING
    if LevToBinFlg == 1
      [LevToBin, qnToBin] = ReadLevToBin(NLevels, Levelvqn, Leveljqn, LevelEeV, DeltaEintDiss, StpInstants, GroupNb);
      BinEnergy    = 0.d0 .* EeV;
      BinEnergyMax = 0.d0 .* EeV;
    elseif LevToBinFlg == 2
      [LevToBin, BinEnergy, BinEnergyMax, qnToBin] = ComputeLevToBin(NLevels, Levelvqn, Leveljqn, LevelEeV, DeltaEintDiss, DeltaEintDepth);
    else
      LevToBin     = 0.d0 .* LevelEeV + 1;
      BinEnergy    = 0.d0 .* EeV;
      BinEnergyMax = 0.d0 .* EeV;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    if PlotProcProbFlg == 1 || PlotEnergiesFlg == 1 || PlotMovieFlg == 1 || PlotStepsFlg == 1 || PlotQSSDistrFlg == 1
      if PlotMovieFlg == 1 || PlotStepsFlg == 1 || PlotEnergiesFlg == 1  
        [Steps] = ComputeMovieSteps(iT, t);  
      end

      if PlotEnergiesFlg == 1
        [eInt, eRot, eVib] = ComputeEnergies(NLevels, LevelQ, Levelg, Pop, QBins, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelvqn, LevelToBin, Steps);
      end
      
      if PlotQSSDistrFlg == 1
         [qssInt, qssRot, qssVib] = ComputeQssDistr(NLevels, LevelQ, Levelg, Pop, QBins, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelvqn, LevelToBin, Steps); 
      end
    else
      vEeVVib = [];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    ReUpload = 0
  end
   

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% PLOTTING

  % Plotting MOLE FRACTIONS
  if PlotMolFracsFlg == 1
    [iFigure] = PlotMolFracs(iT, iFigure, t, MolFracs);     
  end


  % Plotting Stochastic MOLE FRACTIONS
  if PlotMolFracsStochFlg == 1
    [iFigure] = PlotMolFracsStoch(iT, iFigure, tStoch, MolFracsStoch);     
  end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              

  % Plotting TEMPERATURES
  if PlotTemperaturesFlg == 1
    [iFigure] = PlotTemperatures(iT, iFigure, t, Tint, Ttrans);
  end 


  % Plotting QSS RATES
  if PlotKsFlg == 1 
    [iFigure] = PlotKs(iT, iFigure, t, ProcessesRatesOverall);  
  end


  % QSS RATES PLUS MOLE FRACTION
  if PlotKsMoleFracsFlg == 1 
    [iFigure] = PlotKsMoleFracs(iT, iFigure, t, ProcessesRatesOverall, MolFracs);                    
  end

  
  % Plotting MOVIES
  if PlotMovieFlg == 1
    iFigure = PlotMovie(iT, iFigure, t, MolFracs, ProcessesRatesOverall, NLevels, LevelEeV, DeltaEintDiss, LevelQ, Levelg, Levelvqn, Leveljqn, LevelToBin, Pop, QBins, Steps, LevToBin, ProcessesRates);  
  end
  

  % Plotting STEPS
  if PlotStepsFlg == 1
    [iFigure] = PlotSteps(iT, iFigure, t, MolFracs, ProcessesRatesOverall, NLevels, LevelEeV, LevelEeV0, LevelQ, Levelg, LevelToBin, Pop, QBins, StpInstants, Levelvqn, Leveljqn, LevToBin, DeltaEintDiss, ProcessesRates);
  end
  
  
  % Plotting Energies
  if PlotEnergiesFlg == 1
    [iFigure] = PlotEnergies(iT, iFigure, t, eInt, eRot, eVib, MolFracs, P, Steps);
  end
  
  
  % Plotting Processes Probabilities
  if PlotProcProbFlg == 1
    [iFigure] = PlotProcProb(iT, iFigure, t, ProcessesRates, ProcessesRatesOverall, niRatio, LevelEeV, DeltaEintDiss, vEeVVib, Steps, StpInstants, Leveljqn, Levelvqn);
  end
  
  
  % Reconstructing Q.B. Levels from updated Saha Eq.
  if PlotNewSahaFlg == 1
    [iFigure] = PlotNewSaha(iT, iFigure, iFigureStart, StpInstants, NLevels, Levelg, LevelEeV, DeltaEintDiss, LevToBin, BinEnergy, BinEnergyMax, t, MolFracs, nd, VMax, ProcessesRates, RatesMatrix, Tint)
  end 
  
  
  if PlotTausFlg == 1
    [iFigure] = PlotTaus(iFigure);
  end
  
  
  if WriteNetworkCSVFlg == 1
    WriteNetworkCSV(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, Levelg, RatesMatrix, ProcessesRates, t, MolFracs, nd, Pop, StpInstants, LevelEeVVib0, LevelEeVRot)
  end
  
  
  if WriteRatesToParaviewCSVFlg == 1
    WriteRatesToParaviewCSV(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, rOut, Levelg, RatesMatrix, ProcessesRates, t, MolFracs, nd, LevelEeVVib0, LevelEeVRot)
  end
  
  if ComputeMasterEquationFlg == 1
   ComputeMasterEquation(iT, NLevels, LevelEeV, LevelEeV0, Levelvqn, Leveljqn, rIn, rOut, Levelg, RatesMatrix, ProcessesRates, t, MolFracs, nd, LevelEeVVib0, LevelEeVRot)
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
end