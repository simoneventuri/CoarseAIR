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
ReUpload = 2
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
       PlotTausFlg LevToBinFlg WriteNetworkCSVFlg WriteRatesToParaviewCSVFlg ComputeMasterEquationFlg PlotMolFracsStochFlg ...
       PlotPColorFlg PlotSamplesFlg PlotMomentsFlg iPESStart iPESEnd PlotPostCut TimePostCutVec WriteStochQoIsFlg
     
global NStepsMovie TimeMinMovie TimeMaxMovie NSubplots NInstantsPerSubplot tInstants MinEvPlot MaxEvPlot StepsOverlappingSteps vqnColor SubplotsFlg ...
       QNSpace RxLxIdx ReadAllRatesFlg SaveRatesFlg CompOI PlotPairUser BinnedFlg MovieFlg WriteRatesFlag WriteSrcTermsFlag WriteInelasticFlag WriteAllFlag ...
       WriteDissFlag MergePairsFlg tQSS
     
global PlotKsStochFlg PlotEnergiesStochFlg ComputePESToMolFracCorrFlg ComputePESToQSSRatesCorrFlg ComputePESToTausCorrFlg 
global NHL Network_Folder AbscissaConverter PreLogShift BondOrderFun PIPFun AnglesVecGrid RMinGrid RMaxGrid NPointsGrid 
     
global xDim yDim TotDim PlotPair ProcToLevIP TempChar TempChar2 Pair_Name Pair_To_Molecule Pair_to_Atoms iPInternal iPExternal SaveFigs

global RCVec BCVec GCVec KCVec OCVec PCVec WCVec JCVec YCVec CCVec MCVec


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SPECIFING & LOADING INPUT 
PathToOutput = '../Test/'
ReadInput()
System  = 'O3'
T0_Vec  = [10000];
NBins   = [ 200];
KinMthd = ['STS'];
UpdateAllInput() 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHOOSING POSTPROCESSING FUNCTIONS
PlotMolFracsFlg            = 0
PlotTemperaturesFlg        = 0
PlotKsFlg                  = 0
PlotKsMoleFracsFlg         = 0
PlotMovieFlg               = 0
PlotStepsFlg               = 0
PlotEnergiesFlg            = 0
PlotQSSDistrFlg            = 0
PlotProcProbFlg            = 0
PlotNewSahaFlg             = 0
PlotTausFlg                = 0
LevToBinFlg                = 0
WriteNetworkCSVFlg         = 0
WriteRatesToParaviewCSVFlg = 0
ComputeMasterEquationFlg   = 0
PlotMolFracsStochFlg       = 0
  NPESTot         = 50 
  iPESStart       = 1
  iPESEnd         = 50
  PlotPColorFlg   = 0
  PlotSamplesFlg  = 1
  PlotMomentsFlg  = 0
  PlotPostCut     = 1
  TimePostCutVec  = [1.e-6]
PlotKsStochFlg              = 0
  tQSS            = 50.0;
  %[50.0, 2.e-4, 1.5e-5, 4.e-6]
PlotEnergiesStochFlg        = 1

ComputePESToMolFracCorrFlg  = 0
  MoleFracOI          = 0.5
  TopPercentage_Mol   = 5;
  NbSampledPoints_Mol = 100;
  DistThreshold_Mol   = 0.35;
ComputePESToQSSRatesCorrFlg = 0
  TopPercentage_QSS   = [0,3,0,3];
  NbSampledPoints_QSS = [0,100,0,100];
  DistThreshold_QSS   = [0.0,0.3,0.0,0.3];
ComputePESToTausCorrFlg     = 0
  StrExch             = 'WExch';
  TopPercentage_Tau   = [3,0];
  NbSampledPoints_Tau = [100,0];
  DistThreshold_Tau   = [0.3,0.0];
  ReadPESFlg        = true
  NHL               = [6,20,10,1];
  PreLogShift       = -3.5
  Network_Folder    = '/home/aracca/WORKSPACE/CoarseAIR/coarseair/dtb/O3/PESs/BNN/PES9_AbInitio_20_10_LHS50/'
  AbscissaConverter = 1.0
  BondOrderFun      = 'MorseFun'
  PIPFun            = 'Simone'
  AnglesVecGrid     = sort([[35:5:175],[106.75:10:126.75]])
  RMinGrid          = 1.5
  RMaxGrid          = 8.0
  NPointsGrid       = 150

WriteStochQoIsFlg           = 0
  NStepsMovie           = 800
  TimeMinMovie(1:3)     = 1.d-12;
  TimeMaxMovie(1)       = 5.e-5; %[1.e2, 5.e-4, 5e-4, 5.e-5;
  TimeMaxMovie(2)       = 5.e-5; %2.3E-07;
  TimeMaxMovie(3)       = 5.e-5; %2.3E-07;
  StartBin              = 1
  FinalBin              = 200;
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
  MergePairsFlg         = 1
  iComp                 = 2
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

RunKONIGDir = 'RunKonig_InelExch';
OutputPath  = strcat(PathToOutput,RunKONIGDir,'/output/'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTING SYSTEM MAIN QUANTITIES
ComputeSystemQuantities() 
    
for iT = 1:length(T0_Vec)
    
  
  if ReUpload > 0 

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% READING RATES
    if PlotKsStochFlg == 1 || PlotEnergiesStochFlg == 1
      NPESs = NPESTot
      
      if ReUpload >= 1
        RatesMatrixStoch         = zeros(NBins(1),max(yDim),3,NPESs); 
        RatesMatrixStochMmnts    = zeros(NBins(1),max(yDim),3,NPESs); 
        DissRatesStoch           = zeros(NBins(1),NPESs);
        DissRatesStochMmnts      = zeros(NBins(1),NPESs);
        ProcessesRatesStoch      = zeros(NBins(1),4,NPESs); 
        ProcessesRatesStochMmnts = zeros(NBins(1),4,NPESs); 
      end
      
      [RatesMatrixStoch, DissRatesStoch, ProcessesRatesStoch, RatesMatrixStochMmnts, DissRatesStochMmnts, ProcessesRatesStochMmnts] = ReadRatesStoch(iT, RatesMatrixStoch, DissRatesStoch, ProcessesRatesStoch, RatesMatrixStochMmnts, DissRatesStochMmnts, ProcessesRatesStochMmnts, StartBin, FinalBin);
    else
      
      if ReUpload >= 1
        RatesMatrix         = zeros(NBins(1),max(yDim),3,NTint); 
        %RatesSigmaMatrix    = zeros(NBins(1),max(yDim),3,NTint); 
        DissRates           = zeros(NBins(1),NTint);
        %DissRatesSigma      = zeros(NBins(1),NTint);
        ProcessesRates      = zeros(NBins(1),4,NTint); 
        %ProcessesRatesSigma = zeros(NBins(1),4,NTint); 
        ProcessesRatesOverall = [];
      end
      
      if PlotProcProbFlg == 1 || PlotKsFlg == 1 || PlotKsMoleFracsFlg == 1 || PlotMovieFlg == 1 || WriteNetworkCSVFlg == 1 || ComputeMasterEquationFlg == 1
        [RatesMatrix, DissRates, ProcessesRates] = ReadRates(iT, RatesMatrix, DissRates, ProcessesRates, StartBin, FinalBin, DissCorrectionFactor);
      end
      
    end
    
    if ComputePESToMolFracCorrFlg == 1 || ComputePESToQSSRatesCorrFlg == 1 || ComputePESToTausCorrFlg == 1
      [RGrid, EGrid, EGridMean, EGridSD] = ComputeStochPES(ReadPESFlg);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% READ BOX
    if PlotMolFracsStochFlg == 1 || PlotEnergiesStochFlg == 1 || PlotKsStochFlg == 1
      NPESs = NPESTot       
      [tStoch, MolFracsStoch, TtransStoch, rhoStoch, PStoch, ndStoch, EStoch] = ReadBoxStoch(iT); 
    else
      [t, MolFracs, Ttrans, rho, P, nd, E] = ReadBox(iT); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
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
    if PlotKsStochFlg == 1 || PlotEnergiesStochFlg == 1
      if ReUpload == 1
        [PopVecStoch, PopOverQStoch, PopStoch, NStepsStoch] = ReadPopStoch(iT, QBins);    
        [niRatioStoch, ProcessesRatesOverallStoch] = ComputeProcessesRatesOverallStoch(iT, PopStoch, ProcessesRatesStoch);
        filenamePop = strcat(OutputPath,'/T_',num2str(T0_Vec(iT)),'/output/','/Pop')
        save(filenamePop, 'PopVecStoch', 'PopOverQStoch', 'PopStoch', 'NStepsStoch', 'niRatioStoch', 'ProcessesRatesOverallStoch', '-v7.3')
      elseif  ReUpload == 2
        clear PopVec PopOverQ Pop NSteps niRatio ProcessesRatesOverall
        filenamePop = strcat(OutputPath,'/T_',num2str(T0_Vec(iT)),'/output/','/Pop')
        load(filenamePop, 'PopVecStoch', 'PopOverQStoch', 'PopStoch', 'NStepsStoch', 'niRatioStoch', 'ProcessesRatesOverallStoch')
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
    if PlotEnergiesStochFlg == 1
      [StpInstants] = ComputeStepInstants(tStoch);  
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
    

    if PlotProcProbFlg == 1 || PlotEnergiesFlg == 1 || PlotMovieFlg == 1 || PlotStepsFlg == 1 || PlotQSSDistrFlg == 1 || PlotEnergiesStochFlg == 1  
      if PlotMovieFlg == 1 || PlotStepsFlg == 1 || PlotEnergiesFlg == 1
        [Steps] = ComputeMovieSteps(iT, t);  
      elseif  PlotEnergiesStochFlg == 1
        [Steps] = ComputeMovieSteps(iT, tStoch); 
      end

      if PlotEnergiesFlg == 1
        [eInt, eRot, eVib] = ComputeEnergies(NLevels, LevelQ, Levelg, Pop, QBins, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelvqn, LevelToBin, Steps);
      end
      
       if PlotEnergiesStochFlg == 1  
        [eIntStoch, eRotStoch, eVibStoch] = ComputeEnergiesStoch(NLevels, LevelQ, Levelg, PopStoch, QBins, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelvqn, LevelToBin, Steps);
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
  if PlotKsStochFlg == 1 || ComputePESToQSSRatesCorrFlg == 1
     [iFigure, KDiss, PDF_KDiss, KDissQSS, PDF_KDissQSS, KExch, PDF_KExch, KExchQSS, PDF_KExchQSS] = PlotKsStoch(iT, iFigure, tStoch, ProcessesRatesOverallStoch); 
     FileName = strcat('StochRates_', num2str(T0_Vec(1)), 'K' );
     save(FileName, 'KDiss', 'PDF_KDiss', 'KDissQSS', 'PDF_KDissQSS', 'KExch', 'PDF_KExch', 'KExchQSS', 'PDF_KExchQSS', '-v7.3')
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
  
  if PlotEnergiesStochFlg == 1
    [iFigure, TauPressure, tauIntStoch, tauIntPStoch, PDF_tauIntP, tauVibStoch, tauVibPStoch, PDF_tauVibP, tauRotStoch, tauRotPStoch, PDF_tauRotP] = PlotEnergiesStoch(iT, iFigure, tStoch, eIntStoch, eRotStoch, eVibStoch, MolFracsStoch, PStoch, Steps)
    if strcmp(StrExch, 'WExch')
      tauIntPStoch_WExch = tauIntPStoch;
      PDF_tauIntP_WExch  = PDF_tauIntP;
      tauVibPStoch_WExch = tauVibPStoch;
      PDF_tauVibP_WExch  = PDF_tauVibP;
      tauRotPStoch_WExch = tauRotPStoch;
      PDF_tauRotP_WExch  = PDF_tauRotP;
      FileName = strcat('StochTausWExch_', num2str(T0_Vec(1)), 'K' );
      save(FileName, 'tauIntPStoch_WExch', 'PDF_tauIntP_WExch', 'tauVibPStoch_WExch', 'PDF_tauVibP_WExch', 'tauRotPStoch_WExch', 'PDF_tauRotP_WExch', 'TauPressure', '-v7.3');  
    else
      tauIntPStoch_WOExch = tauIntPStoch;
      PDF_tauIntP_WOExch  = PDF_tauIntP;
      tauVibPStoch_WOExch = tauVibPStoch;
      PDF_tauVibP_WOExch  = PDF_tauVibP;
      tauRotPStoch_WOExch = tauRotPStoch;
      PDF_tauRotP_WOExch  = PDF_tauRotP;
      FileName = strcat('StochTausWOExch_', num2str(T0_Vec(1)), 'K' );
      save(FileName, 'tauIntPStoch_WOExch', 'PDF_tauIntP_WOExch', 'tauVibPStoch_WOExch', 'PDF_tauVibP_WOExch', 'tauRotPStoch_WOExch', 'PDF_tauRotP_WOExch', 'TauPressure', '-v7.3');  
    end
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
  
  
  if ComputePESToMolFracCorrFlg == 1
    [iFigure] = ComputePESToMoleFracCorr(iFigure, iComp, MoleFracOI, tStoch, MolFracsStoch, RGrid , EGrid, EGridMean, EGridSD, TopPercentage_Mol, NbSampledPoints_Mol, DistThreshold_Mol);
  end
  
  if ComputePESToQSSRatesCorrFlg == 1
    [iFigure] = ComputePESToQSSRatesCorr(iFigure, KDiss, KDissQSS, KExch, KExchQSS, RGrid, EGrid, EGridMean, EGridSD, TopPercentage_QSS, NbSampledPoints_QSS, DistThreshold_QSS);
  end
  
  if ComputePESToTausCorrFlg == 1
    [iFigure] = ComputePESToTausCorr(iFigure, tauIntPStoch, tauVibPStoch, tauRotPStoch, RGrid, EGrid, EGridMean, EGridSD, TopPercentage_Tau, NbSampledPoints_Tau, DistThreshold_Tau, StrExch);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  if WriteStochQoIsFlg == 1
    
    TVec = [2500, 5000, 7500, 10000]
    
    iT=1;
    for TT = TVec
      FileName = strcat('StochRates_', num2str(TT), 'K' );
      load(FileName, 'KDiss', 'PDF_KDiss', 'KDissQSS', 'PDF_KDissQSS', 'KExch', 'PDF_KExch', 'KExchQSS', 'PDF_KExchQSS');
      
      KDiss_All(iT,:)    = KDiss(:);
      KDissQSS_All(iT,:) = KDissQSS(:);
      KExch_All(iT,:)    = KExch(:);
      KExchQSS_All(iT,:) = KExchQSS(:);
      
      iT=iT+1;
    end
    
    figure
    semilogy(10000./TVec,KDiss_All,'k')
    hold on
    semilogy(10000./TVec,KDissQSS_All,'r')
    
    figure
    semilogy(10000./TVec,KExch_All,'k')
    hold on
    semilogy(10000./TVec,KExchQSS_All,'r')
    
    
    iT=1;
    for TT = TVec
      FileName = strcat('StochTausWExch_', num2str(TT), 'K' );
      load(FileName, 'tauIntPStoch_WExch',  'PDF_tauIntP_WExch',  'tauVibPStoch_WExch',  'PDF_tauVibP_WExch',  'tauRotPStoch_WExch',  'PDF_tauRotP_WExch',  'TauPressure');
    
      tauRotP_All(iT,:) = tauRotPStoch_WExch(:);
      tauVibP_All(iT,:) = tauVibPStoch_WExch(:);
      
      iT=iT+1;
    end
    
    figure
    semilogy(TVec,tauRotP_All,'k')
    hold on
    semilogy(TVec,tauVibP_All,'r')
      
    
    for TT = TVec
      FileName = strcat('StochRates_', num2str(TT), 'K' );
      load(FileName, 'KDiss', 'PDF_KDiss', 'KDissQSS', 'PDF_KDissQSS', 'KExch', 'PDF_KExch', 'KExchQSS', 'PDF_KExchQSS');
      FileName = strcat('StochTausWExch_', num2str(TT), 'K' );
      load(FileName, 'tauIntPStoch_WExch',  'PDF_tauIntP_WExch',  'tauVibPStoch_WExch',  'PDF_tauVibP_WExch',  'tauRotPStoch_WExch',  'PDF_tauRotP_WExch',  'TauPressure');
    
      WriteStochQoIs(TT, TauPressure, KDiss, PDF_KDiss, KDissQSS, PDF_KDissQSS, KExch, PDF_KExch, KExchQSS, PDF_KExchQSS, ...
                     tauIntPStoch_WExch,  PDF_tauIntP_WExch,  tauVibPStoch_WExch,  PDF_tauVibP_WExch,  tauRotPStoch_WExch,  PDF_tauRotP_WExch);
    end
    
    filename = './VargaQoIs.csv';
    delimiter = ',';
    startRow = 2;
    formatSpec = '%s%s%s%s%s%s%s%s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
      raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));
    for col=[1,2,3,4,5,6,7,8]
      % Converts text in the input cell array to numbers. Replaced non-numeric
      % text with NaN.
      rawData = dataArray{col};
      for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
          result = regexp(rawData(row), regexstr, 'names');
          numbers = result.numbers;

          % Detected commas in non-thousand locations.
          invalidThousandsSeparator = false;
          if numbers.contains(',')
            thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
            if isempty(regexp(numbers, thousandsRegExp, 'once'))
              numbers = NaN;
              invalidThousandsSeparator = true;
            end
          end
          % Convert numeric text to numbers.
          if ~invalidThousandsSeparator
            numbers = textscan(char(strrep(numbers, ',', '')), '%f');
            numericData(row, col) = numbers{1};
            raw{row, col} = numbers{1};
          end
        catch
          raw{row, col} = rawData{row};
        end
      end
    end
    TVec         = cell2mat(raw(:, 1));
    KDissQSS     = cell2mat(raw(:, 2));
    KDissThermal = cell2mat(raw(:, 3));
    KExchQSS     = cell2mat(raw(:, 4));
    KExchThermal = cell2mat(raw(:, 5));
    TauVib = cell2mat(raw(:, 6));
    TauRot = cell2mat(raw(:, 7));
    TauInt = cell2mat(raw(:, 8));
    clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp;

    figure
    semilogy(10000./TVec, KDissQSS)
    hold on
    semilogy(10000./TVec, KDissThermal)

    figure
    semilogy(10000./TVec, KExchQSS)
    hold on
    semilogy(10000./TVec, KExchThermal)

    figure
    semilogy(TVec, TauVib)
    hold on
    semilogy(TVec, TauRot)

  end
  
end