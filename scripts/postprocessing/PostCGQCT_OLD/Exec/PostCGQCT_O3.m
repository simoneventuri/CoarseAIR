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
% untitled
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
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
ReUpload = 1               % =0 for not reading data; =1 ror reading .dat data; =2 for reading .m data 
if ReUpload == 1
  clear all
  ReUpload = 1
end

global PathToOutput System T0_Vec NBins KinMthd NTtra NTint TInt_Vec NMolecules StartBin FinalBin NBinnedMol BinnedMolName NPESs

global NAtoms AtomsName MoleculesName DegeneracyFactor ColPartToComp BinnedMolToComp NComp CompNames CompColor AtomColor AtomSize AllMoleculesName ...
       PairColor AtomMass ComponentMass ColorVec ComponentDeg NComp RxLxIdx
     
global Plnck UKb Ue KeV AvN AMUToKg EhToeV DSWtoKg
global SystemPath RatesPath DatabasePath OutputPath SaveFigs FigDirPath linS linST AxisFontSz AxisFontNm AxisLabelSz AxisLabelNm ...
       LegendFontSz LegendFontNm XLimPlot YLimPlot 

global PlotDiatPotFlg PlotPESFlg PerturbEnergyFlg PlotTrajectoryFlg AnalyzeTrajectoriesFlg CheckbSensitivityFlg ProduceMatFlg PlotLevelRatesFlg ...
       PlotProcessesRatesFlg PlotStSRatesQNFlg PlotRatesMatEeVFlg PlotRatesMovieFlg PlotRatesMovieCompFlg WriteRatesToParaviewCSVFlg LevToBinFlg ComputeEigFlg ...
       CheckEmptyRatesFlg WriteRatesToNetworkFlg ComputePESToRatesCorrFlg PlotStochProcessesRatesFlg

global SigmaOn MergePairsFlg PlotPairUser iTraj iNode iProc tMin tMax vqnColor Angles MinPES MaxPES iLevelToPerturb EnergyUM EeVSpace iLevelsVec ...
       SubplotsFlg OverlapFlg NBinsNew NBinsNewQB NBinsNewBLow MinEnQB NBinsNewBBack jqnBack ReadAllRatesFlg SaveRatesFlg tInstants PairsToSum ...
       OutputFileName ForceFlg ParaViewFlg ConsiderDegFlg WriteExchangeFlg WriteOrigInvFlg WritePajekFlg WriteDirUndFlg
   
global xDim yDim TotDim PlotPair ProcToLevIP TempChar TempChar2 Pair_Name Pair_To_Molecule Pair_to_Atoms iPInternal iPExternal

global RCVec BCVec GCVec KCVec OCVec PCVec WCVec JCVec YCVec CCVec MCVec

     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SPECIFING & LOADING INPUT 
PathToOutput = '../Test'
ReadInput()
System  = 'O3'
T0_Vec  = [2500];
NBins   = [ 6115];
KinMthd = [ 'STS'];
UpdateAllInput() 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHOOSING POSTPROCESSING FUNCTIONS
PlotDiatPotFlg             = 0
PlotPESFlg                 = 0
PerturbEnergyFlg           = 0
PlotTrajectoryFlg          = 0
AnalyzeTrajectoriesFlg     = 0
CheckbSensitivityFlg       = 0
CheckEmptyRatesFlg         = 0
ProduceMatFlg              = 1
ComputeEigFlg              = 0
PlotLevelRatesFlg          = 0
PlotProcessesRatesFlg      = 1
PlotStSRatesQNFlg          = 0
PlotRatesMatEeVFlg         = 0
PlotRatesMovieFlg          = 0
PlotRatesMovieCompFlg      = 0
WriteRatesToParaviewCSVFlg = 0
WriteRatesToNetworkFlg     = 0
LevToBinFlg                = 0
ComputePESToRatesCorrFlg   = 0
PlotStochProcessesRatesFlg = 0
  SigmaOn         = ' no'
  MergePairsFlg   = 1
  iBins           = 1
  StartBin        = 1
  FinalBin        = 6115
  PlotPairUser    = [1, 1, 1]
  iTraj           = 1
  iNode           = 1
  iProc           = 1
  tMin            = 0.d0
  tMax            = 200000.d0
  vqnColor        = 0
  Angles          = [120]
  MinPES          = [1.5, 1.5]
  MaxPES          = [  5,   5]
  iLevelToPerturb = [1000]
  EnergyUM        = 'eV'
  EeVSpace        = 1
  iLevelsVec      = [1, 1500, 3000, 4500, 6000, 7500, 9000, 10500, 13000];
  SubplotsFlg     = 1
  OverlapFlg      = 0
  NBinsNew        = [   50,   49]
  NBinsNewQB      = [   10,    4] 
  NBinsNewBLow    = [    3,    3]
  MinEnQB         = [ -1.5, -2.0]
  NBinsNewBBack   = [    0,    1]
  jqnBack         = [  500,  100]
  ReadAllRatesFlg = 0
  tInstants       = [1.e-14, 1.e-12, 1.e-10, 1.e-9, 5.e-9, 1.e-8, 3.e-6, 1.e-3];
  PairsToSum      = [2];
  OutputFileName  = 'DiatomicPotential_NASA.csv';
  ForceFlg        = 0
  ParaViewFlg     = 0
  
  ConsiderDegFlg   = 0      % = 1 for Considering the Master Eq. in n; = 0 for Considering the Master Eq. in n/g
  WriteExchangeFlg = 1      % = 1 for Considering Exchange
  WriteOrigInvFlg  = [1,0]  % = [1,*] for writing Kij; = [*,1] for writing 1/Kij
  WritePajekFlg    = 1      % = 0 for writing a .csv list of edges; = 1 for writing a Pajek .net file 
  WriteDirUndFlg   = [1,0] 
  
  DissCorrectionFactor  = 1.0 %16.0/3.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UPDATE PARAMETERS & PATHS
[iFigure] = UpdateParameters();
  
RatesPath   = strcat('../Test/',System(1,:),'/',MoleculesName(1,:),'/Rates/');                                                             % For Levels / Bins Rates
SystemPath  = strcat('../Test/',System(1,:));                                                                                              % For qnsEnBin.dat (Levels Info)

iFigure    = 2001;
SaveFigs   = 0
FigDirPath = strcat('./CGQCT-',System,'-Figures/');
if SaveFigs > 0
  [status,msg,msgID] = mkdir(FigDirPath)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTING SYSTEM MAIN QUANTITIES
ComputeSystemQuantities() 
if ReUpload >= 1
  RatesMatrix         = zeros(NBins(1),max(yDim),3,NTint); 
  %RatesSigmaMatrix    = zeros(NBins(1),max(yDim),3,NTint); 
  DissRates           = zeros(NBins(1),NTint);
  %DissRatesSigma      = zeros(NBins(1),NTint);
  ProcessesRates      = zeros(NBins(1),4,NTint); 
  %ProcessesRatesSigma = zeros(NBins(1),4,NTint); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
for iT = 1:length(T0_Vec)
  
  
  if PlotTrajectoryFlg == 1
    [iFigure] = PlotTrajectory(iFigure, iT)
  end 
  
  
  if ReUpload > 0
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% READING AND COMPUTING QUANTITIES
    %if 
      [QBins, EeV] = ReadPartFunEnergy(iT);
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% READING AND COMPUTING QUANTITIES
    %if 
      [NLevels, Levelvqn, Leveljqn, LevelEh, LevelEeV, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelg, LevelToBin, LevelQ, vToLevel, DeltaEintDiss, DeltaEintDepth, rIn, rOut, rMin, rMax, VMin, VMax, Tau, Egam, dVIn, ddVIn, dVOut, ddVOut, dVdJIn, dVdJOut] = ReadLevelInfo(iT, EeV);
      if ReUpload < 2
        for iMol=size(LevelEeV,2)
          EeV(:,iMol) = EeV(:,iMol) + LevelEeV(1,iMol);
        end
      end
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% READING RATES
    if CheckEmptyRatesFlg == 1
     CheckEmptyRates(iT, StartBin, FinalBin);
    end 
    
    if ComputeEigFlg == 1 || ProduceMatFlg == 1 || PlotLevelRatesFlg == 1 || PlotRatesMovieCompFlg == 1 || PlotRatesMovieFlg == 1 || PlotRatesMatEeVFlg == 1 || PlotStSRatesQNFlg == 1 || PlotProcessesRatesFlg || WriteRatesToParaviewCSVFlg == 1 || WriteRatesToNetworkFlg == 1
       [RatesMatrix, DissRates, ProcessesRates, KRec, KRecOverg] = ReadRates(iT, RatesMatrix, DissRates, ProcessesRates, StartBin, FinalBin, DissCorrectionFactor, NLevels, Levelg, LevelEeV);
    end
    
    if ComputePESToRatesCorrFlg == 1 || PlotStochProcessesRatesFlg
      NPESs = 100
      if ReadAllRatesFlg == -1
        filenamePES = strcat('./StochPES')
        load(filenamePES, 'StochPESR', 'StochPESAng', 'StochPESEeV')
      else
        [StochPESR, StochPESAng, StochPESEeV] = ReadStochPES();
      end 
      
      RatesMatrixStoch         = zeros(NBins(1),max(yDim),3,NPESs); 
      %RatesSigmaMatrixStoch    = zeros(NBins(1),max(yDim),3,NPESs); 
      DissRatesStoch           = zeros(NBins(1),NPESs);
      %DissRatesSigmaStoch      = zeros(NBins(1),NPESs);
      ProcessesRatesStoch      = zeros(NBins(1),4,NPESs); 
      %ProcessesRatesSigmaStoch = zeros(NBins(1),4,NPESs); 
      [RatesMatrixStoch, DissRatesStoch, ProcessesRatesStoch] = ReadRatesStoch(iT, RatesMatrixStoch, DissRatesStoch, ProcessesRatesStoch, StartBin, FinalBin);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% READING / COMPUTING LEVELS-To-BINs MAPPING
    if LevToBinFlg == 1
      [LevToBin, qnToBin] = ReadLevToBin(NLevels, Levelvqn, Leveljqn, LevelEeV, DeltaEintDiss);
    elseif LevToBinFlg == 2
      [LevToBin, BinEnergy, BinEnergyMax, qnToBin] = ComputeLevToBin(NLevels, Levelvqn, Leveljqn, LevelEeV, DeltaEintDiss, DeltaEintDepth);
    else
      LevToBin = 0.d0 .* LevelEeV + 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    ReUpload = 0
  end
  
  
  if AnalyzeTrajectoriesFlg == 1
    [iFigure, xxFin, EintFin, RbondFin, AngMomFin, jqnFin, vqnFin, iTypeFin] = AnalyzeTrajectories(iFigure, iT, iBins, Levelvqn, Leveljqn, LevelEeV)
  end
  
  
  if CheckbSensitivityFlg == 1
    [iFigure, FinState, Cross, CrossCDF, CrossVar] = CheckbSensitivity(iFigure, iT, iBins, Levelvqn, Leveljqn)
  end 
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% WRITING CSV FILES
  if WriteRatesToParaviewCSVFlg == 1
    WriteRatesToParaviewCSV(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, rOut, Levelg, RatesMatrix, ProcessesRates)
  end
  
  if WriteRatesToNetworkFlg == 1
    %WriteRatesToNetwork(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, Levelg, RatesMatrix, ProcessesRates, LevelEeVVib0, LevelEeVRot, DeltaEintDiss)
    %WriteRatesToNetwork_Infomap(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, Levelg, RatesMatrix, ProcessesRates, LevelEeVVib0, LevelEeVRot, DeltaEintDiss)
    WriteRatesAsNetwork_ForAmal(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, Levelg, RatesMatrix, ProcessesRates, LevelEeVVib0, LevelEeVRot, DeltaEintDiss)
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING DIATOMIC POTENTIALS & PES
if PlotDiatPotFlg > 0
  [iFigure] = PlotDiatPot(iFigure, NLevels, Leveljqn, Levelvqn, LevelEh, LevelEeV, rMin, rMax, VMin, VMax, rIn, rOut, LevToBin, DeltaEintDiss, LevelEeVVib0, LevelEeVRot)
end 

if PlotPESFlg == 1
  [iFigure] = PlotPES(iFigure)
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERTURBING DIATOMIC POTENTIAL
if PerturbEnergyFlg == 1
  [iFigure] = PerturbEnergy(iFigure, Levelvqn, Leveljqn, LevelEeV);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTING EIGEN-FUNCTIONS
if ComputeEigFlg > 0
  [iFigure, EigenVals] = ComputeEig(iFigure, iT, RatesMatrix, NLevels, Levelg, LevelEeV);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING QUANTITIES

if PlotLevelRatesFlg == 1
  %iFigure = PlotLevelRates(iFigure, NLevels, StartBin, FinalBin, RatesMatrix, RatesSigmaMatrix, LevelEeV, Levelvqn)
  iFigure = PlotLevelRates(iFigure, NLevels, StartBin, FinalBin, RatesMatrix, LevelEeV, Levelvqn)
end

% Plotting Dissociation Rates
if PlotProcessesRatesFlg == 1
  %iFigure = PlotProcessesRates(iFigure, ProcessesRates, ProcessesRatesSigma, LevelEeV, EeV, LevToBin);  
  iFigure = PlotProcessesRates(iFigure, NLevels, ProcessesRates, Levelg, LevelEeV, EeV, DeltaEintDiss, Levelvqn, Leveljqn,LevToBin, rOut, rIn, Egam, QBins)
end

% Plotting Rates Matrix [Ev(iLevel) -> Ev(jLevel)]
if PlotRatesMatEeVFlg == 1
  [iFigure] = PlotRatesMatEeV(iFigure, RatesMatrix, ProcessesRates, Levelg, LevelEeV, EeV)
end

% Plotting Rates [iLevel -> (jLevel(jqn),jLevel(vqn))]
if PlotStSRatesQNFlg == 1
  [iFigure] = PlotStSRatesQN(iFigure, LevelEeV, Leveljqn, Levelvqn, RatesMatrix)
end

for iP=1:3
  if PlotPairUser(iP) ~= 0

    % Plotting 
    if PlotRatesMovieFlg == 1
      [iFigure] = PlotRatesMovie(iFigure, iP, EeV, RatesMatrix)
    end

    % Plotting 
    if PlotRatesMovieCompFlg == 1
      [iFigure] = PlotRatesMovieComp(iFigure, iP, EeV, RatesMatrix)
    end

  end 
end

% Plotting Stochastic Dissociation Rates
if PlotStochProcessesRatesFlg == 1
  iFigure = PlotStochProcessesRates(iFigure, NLevels, ProcessesRatesStoch, Levelg, LevelEeV, EeV, DeltaEintDiss, Levelvqn, Leveljqn, LevToBin, rOut, rIn, Egam, Tau, LevelQ, dVIn, ddVIn, dVOut, ddVOut, dVdJIn, dVdJOut, QBins)
end

if ComputePESToRatesCorrFlg == 1
  IniBin  = 13
  FinBin  = 15
  FinProc = 1
  [iFigure, CorrVec] = ComputePESToRatesCorr(iFigure, StochPESR, StochPESAng, StochPESEeV, RatesMatrixStoch, IniBin, FinBin, FinProc);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%