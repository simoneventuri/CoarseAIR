clear all
close all
clc


global PathToOutput System T0_Vec NBins KinMthd NT  tra NTint TInt_Vec NMolecules StartBin FinalBin NBinnedMol BinnedMolName

global NAtoms AtomsName MoleculesName DegeneracyFactor ColPartToComp BinnedMolToComp NComp CompNames CompColor AtomColor AtomSize AllMoleculesName ...
       PairColor AtomMass ComponentMass ColorVec ComponentDeg
     
global Plnck UKb Ue KeV AvN AMUToKg EhToeV DSWtoKg
global SystemPath RatesPath DatabasePath OutputPath SaveFigs FigDirPath linS linST AxisFontSz AxisFontNm LegendFontSz LegendFontNm

global PlotDiatPotFlg PlotPESFlg PerturbEnergyFlg PlotTrajectoryFlg AnalyzeTrajectoriesFlg CheckbSensitivityFlg ProduceMatFlg PlotLevelRatesFlg ...
       PlotProcessesRatesFlg PlotStSRatesQNFlg PlotRatesMatEeVFlg PlotRatesMovieFlg PlotRatesMovieCompFlg WriteCSVFilesFlg LevToBinFlg

global SigmaOn MergePairsFlg PlotPairUser iTraj iNode iProc tMin tMax vqnColor Angles MinPES MaxPES iLevelToPerturb EnergyUM EeVSpace iLevelsVec ....
       SubplotsFlg OverlapFlg NBinsNew NBinsNewQB NBinsNewBLow MinEnQB NBinsNewBBack jqnBack ReadAllRatesFlg SaveRatesFlg
     
global xDim yDim TotDim PlotPair ProcToLevIP TempChar TempChar2 Pair_Name Pair_To_Molecule Pair_to_Atoms iPInternal iPExternal


     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SPECIFING & LOADING INPUT 
System  = 'CO2'
T0_Vec  = [5000];
MoleculesName = [  'CO';  'O2'];
NBins         = [ 13521;  6078];
KinMthd       = [ 'STS'; 'STS'];
NMolecules    = length(NBins);
PairToMol = [1, 1, 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PathToMapping(1) = {'~/WORKSPACE/CoarseAIR/coarseair/dtb/CO2/CO/LevelsMapping_CO.dat'};
PathToMapping(2) = {'~/WORKSPACE/CoarseAIR/coarseair/dtb/CO2/O2/LevelsMapping_O2_New_to_NASA.dat'};

LevMapping      = zeros(max(NBins),NMolecules);
for iMol = 1:NMolecules
  delimiter = '';
  formatSpec = '%f%[^\n\r]';
  fileID = fopen(char(PathToMapping(iMol)),'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
  fclose(fileID);
  LevMapping(1:NBins(iMol),iMol) = dataArray{:, 1};
  clearvars filename delimiter formatSpec fileID dataArray ans;
end
ProcessToJ      = zeros(sum(NBins(PairToMol(:)))+1,1);
ProcessToJ(1)   = 1;
Start           = 1;
StartBis        = 1;
for iP = 1:3
  ProcessToJ(Start+1:Start+NBins(PairToMol(iP)),1) = LevMapping(1:NBins(PairToMol(iP)),PairToMol(iP)) + StartBis;
  Start    = Start    + NBins(PairToMol(iP));
  StartBis = StartBis + max(LevMapping(:,PairToMol(iP))); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iPES = 2:3
  
  for iT = 1:length(T0_Vec)
    
    DirName              = strcat('~/WORKSPACE/CoarseAIR/RESULTS/CO2/StS_PES',num2str(iPES),'/Test/',System(1,:),'/',MoleculesName(1,:),'/Rates/T_',num2str(T0_Vec(iT)),'_',num2str(T0_Vec(iT)),'_New');
    [status, msg, msgID] = mkdir(DirName)

    for iBin = 1:NBins(1)

      filenameRatesNew = strcat(DirName,'/Bin',num2str(iBin),'.dat');
      fileIDNew        = fopen(filenameRatesNew,'w');

      filenameRates = strcat('~/WORKSPACE/CoarseAIR/RESULTS/CO2/StS_PES',num2str(iPES),'/Test/',System(1,:),'/',MoleculesName(1,:),'/Rates/T_',num2str(T0_Vec(iT)),'_',num2str(T0_Vec(iT)),'/Bin',num2str(iBin),'.dat');
      startRow = 6;
      formatSpec = '%30s%10f%20f%20f%[^\n\r]';
      fileID = fopen(filenameRates,'r');
      dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
      dataArray{1} = strtrim(dataArray{1});
      fclose(fileID);
      String  = dataArray{:, 1};
      Process = dataArray{:, 2};
      Rate    = dataArray{:, 3};
      RateSD  = dataArray{:, 4};
      clearvars startRow formatSpec fileID dataArray ans;
      endRow = 5;
      formatSpec = '%s%[^\n\r]';
      fileID = fopen(filenameRates,'r');
      dataArray = textscan(fileID, formatSpec, endRow, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
      dataArray{1} = strtrim(dataArray{1});
      fclose(fileID);
      HeaderLines = dataArray{:, 1};
      clearvars filename endRow formatSpec fileID dataArray ans;
      startRow = 2;
      endRow = 2;
      formatSpec = '%84s%[^\n\r]';
      fileID = fopen(filenameRates,'r');
      dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
      dataArray{1} = strtrim(dataArray{1});
      fclose(fileID);
      HeaderLine2 = dataArray{:, 1};
      clearvars filename startRow endRow formatSpec fileID dataArray ans;
      
      formatSpec = '%-104s\r\n';
      fprintf(fileIDNew,formatSpec,HeaderLines{1});
      formatSpec = '%-84s%18d\r\n';
      fprintf(fileIDNew,formatSpec,HeaderLine2{1},sum(max(LevMapping(:,PairToMol(:))))+1);
      for iLine = 3:5
        formatSpec = '%-104s\r\n';
        fprintf(fileIDNew,formatSpec,HeaderLines{iLine});
      end

      for iLine = 1:length(Rate)
        Process(iLine) = ProcessToJ(Process(iLine));
        formatSpec     = '%-19s %20d %19.10E %19.10E\n';
        fprintf(fileIDNew,formatSpec,String{iLine},Process(iLine),Rate(iLine),RateSD(iLine));
      end
      
      fclose(fileIDNew);
      
    end

  end

end