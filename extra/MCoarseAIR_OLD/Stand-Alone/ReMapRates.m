clear all
close all
clc


global PathToOutput System T0_Vec NBins KinMthd NTtra NTint TInt_Vec NMolecules StartBin FinalBin NBinnedMol BinnedMolName

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
T0_Vec  = [10000];
MoleculesName = [  'O2';  'CO'];
NBins         = [  6112; 13521];
KinMthd       = [ 'STS'; 'STS'];
NMolecules    = length(NBins);
PairToMol = [1, 2, 2];
NPESs     = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PathToMapping(2) = {'~/WORKSPACE/CoarseAIR/coarseair/dtb/CO2/CO/LevelsMapping_CO.dat'};
PathToMapping(1) = {'~/WORKSPACE/CoarseAIR/coarseair/dtb/CO2/O2/LevelsMapping_O2_venturi_to_new.dat'};

LevMapping = zeros(max(NBins),NMolecules);
for iMol = 1:NMolecules
  delimiter = '';
  formatSpec = '%f%[^\n\r]';
  fileID = fopen(char(PathToMapping(iMol)),'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
  fclose(fileID);
  LevMapping(1:NBins(iMol),iMol) = dataArray{:, 1};
  clearvars filename delimiter formatSpec fileID dataArray ans;
end

ProcessToJ         = zeros(sum(NBins(PairToMol(:)))+1,2);
ProcessToJ(1,:)    = [1,1];
ProcessOldToNew(1) = 1;
StartOldToNew      = 1;
StartToJ           = 1;
for iP = 1:3
  
  ProcessOldToNew(StartToJ+1:StartToJ+NBins(PairToMol(iP)),1) = LevMapping(1:NBins(PairToMol(iP)),PairToMol(iP)) + StartOldToNew;
  StartOldToNew                                               = StartOldToNew + max(LevMapping(:,PairToMol(iP)));
  
  ProcessToJ(StartToJ+1:StartToJ+NBins(PairToMol(iP)),1)      = LevMapping(1:NBins(PairToMol(iP)),PairToMol(iP));
  ProcessToJ(StartToJ+1:StartToJ+NBins(PairToMol(iP)),2)      = PairToMol(iP)+1;
  StartToJ                                                    = StartToJ + NBins(PairToMol(iP));
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Rates      = zeros(max(max(LevMapping)),1+NMolecules);
RatesSD    = zeros(max(max(LevMapping)),1+NMolecules);
RatesPES   = zeros(max(max(LevMapping)),1+NMolecules);
RatesSDPES = zeros(max(max(LevMapping)),1+NMolecules);
      
for iT = 1:length(T0_Vec)
  
%   for iPES = 1:NPESs
%     filenameRatesPESNew(iPES,:) = strcat('~/WORKSPACE/CoarseAIR/RESULTS/CO2/StS_PES',num2str(iPES),'/Test/',System(1,:),'/',MoleculesName(1,:),'/Rates/T_',num2str(T0_Vec(iT)),'_',num2str(T0_Vec(iT)),'/Rates.csv');   
%     fileIDPES(iPES)             = fopen(filenameRatesPESNew(iPES,:),'w');
%   end
%   filenameRatesNew = strcat('~/WORKSPACE/CoarseAIR/RESULTS/CO2/Rates_',num2str(T0_Vec(iT)),'K.csv');  
%   fileIDNew        = fopen(filenameRatesNew,'w');
%   formatSpec       = '# InitialBin, FinalMolecule, FinalBin, Process, Rate, RateSD\n';
%   fprintf(fileIDNew,formatSpec);
  
  for iBin = 1427:NBins(1)
    
    Rates   = 0.d0 * Rates;
    RatesSD = 0.d0 * RatesSD;

    for iPES = 1:NPESs
      
      RatesPES   = 0.d0 * RatesPES;
      RatesSDPES = 0.d0 * RatesSDPES;

      if LevMapping(iBin,PairToMol(1)) ~= 0

        filenameBinRatesNew   = strcat('~/WORKSPACE/CoarseAIR/RESULTS/O2C/StS_Merged_',num2str(T0_Vec(iT)),'/Test/',System(1,:),'/',MoleculesName(1,:),'/Rates_New/T_',num2str(T0_Vec(iT)),'_',num2str(T0_Vec(iT)),'/Bin',num2str(LevMapping(iBin,PairToMol(1))),'.dat');
        fileIDBinRatesNew     = fopen(filenameBinRatesNew,'w');
        formatSpecBinRatesNew = '                            %12d%20.10E%20.10E\n';

        %filenameRates = strcat('~/WORKSPACE/CoarseAIR/RESULTS/CO2/StS_PES',num2str(iPES),'/Test/',System(1,:),'/',MoleculesName(1,:),'/Rates/T_',num2str(T0_Vec(iT)),'_',num2str(T0_Vec(iT)),'/Bin',num2str(iBin),'.dat');
        filenameBinRates = strcat('~/WORKSPACE/CoarseAIR/RESULTS/O2C/StS_Merged_',num2str(T0_Vec(iT)),'/Test/',System(1,:),'/',MoleculesName(1,:),'/Rates/T_',num2str(T0_Vec(iT)),'_',num2str(T0_Vec(iT)),'/Bin',num2str(iBin),'.dat');
        startRow       = 6;
        formatSpec     = '%*28s%12f%20f%20f%[^\n\r]';
        fileIDBinRates = fopen(filenameBinRates,'r');
        dataArray      = textscan(fileIDBinRates, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileIDBinRates);
        TempProcesses  = dataArray{:, 1};
        TempRates      = dataArray{:, 2};
        TempRatesSD    = dataArray{:, 3};
        clearvars startRow formatSpec fileIDBinRates dataArray ans;


        fileIDBinRates = fopen(filenameBinRates);
        for i=1:5
          tline = fgetl(fileIDBinRates);
          if i==2
            tline = [tline(1:97),num2str(StartOldToNew,'               %20d'),'\n'];
          else
            tline = strcat(tline,'\n');
          end
          fprintf(fileIDBinRatesNew,tline);
        end
        fclose(fileIDBinRates);
        for iLine = 1:length(TempProcesses)
          if ProcessToJ(TempProcesses(iLine),1) ~= 0
            fprintf(fileIDBinRatesNew,formatSpecBinRatesNew,ProcessOldToNew(TempProcesses(iLine),1),TempRates(iLine),TempRatesSD(iLine));
          end
        end
        fclose(fileIDBinRatesNew);

  %       for iLine = 1:length(TempRates)
  %         if ProcessToJ(TempProcesses(iLine),1) ~= 0
  %           RatesPES(ProcessToJ(TempProcesses(iLine),1),ProcessToJ(TempProcesses(iLine),2))   = RatesPES(ProcessToJ(TempProcesses(iLine),1),ProcessToJ(TempProcesses(iLine),2))   + TempRates(iLine);
  %           RatesSDPES(ProcessToJ(TempProcesses(iLine),1),ProcessToJ(TempProcesses(iLine),2)) = RatesSDPES(ProcessToJ(TempProcesses(iLine),1),ProcessToJ(TempProcesses(iLine),2)) + TempRatesSD(iLine).^2;
  %         end
  %       end
  %       
  %       Rates   = Rates   + RatesPES;
  %       RatesSD = RatesSD + RatesSDPES;


  %       if Rates(1,1) > 0.d0
  %         formatSpec = '%5d, %1d, %5d, %13E, %13E\n';
  %         fprintf(fileIDNew,formatSpec,iBin,0,0,Rates(1,1),RatesSD(1,1));
  %       end
  %       for iMol = 1:2
  %         for jBin = 1:NBins(iMol)
  %           if Rates(jBin,iMol+1) > 0.d0
  %             formatSpec = '%5d, %1d, %5d, %13E, %13E\n';
  %             fprintf(fileIDNew,formatSpec,iBin,iMol,jBin,Rates(jBin,iMol+1),RatesSD(jBin,iMol+1));
  %           end
  %         end
  %       end
  
      end
      
    end
    
    Rates   =      Rates   ./ NPESs;
    RatesSD = sqrt(RatesSD ./ NPESs);
       
  end

end