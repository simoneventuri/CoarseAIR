%% The Function reads the rates (K(i,j)) and computes the Processes (Dissociation, Exchange 1,2 and 3) overall rates (K(i,Process))
%
%  Input Arguments:  - iT:            Index for the current Translational Temperature
%                    - iLevelsStart:  (Optional) Initial Level/Bin to read. If not given, it is assumed equal to be equal to 1
%                    - iLevelsEnd:    (Optional) Final Level/Bin to read. If not given, it is assumed equal to be equal to NBins
%
%  Input Global Var: - T0_Vec:        Vector of Translational Temperatures (e.g.: [10000])
%                    - RatesPath:     The path to the output folder (e.g.: ../Test/N3/N2/Rates)
%

%function [RatesSigma, RatesMatrix, RatesSigmaMatrix, DissRates, DissRatesSigma, ProcessesRates, ProcessesRatesSigma] = ReadRates(iT, RatesSigma, RatesMatrix, RatesSigmaMatrix, DissRates, DissRatesSigma, ProcessesRates, ProcessesRatesSigma, iBinsStart, iBinsEnd)    
function [RatesMatrix, DissRates, ProcessesRates] = ReadRatesAll(iT, RatesMatrix, DissRates, ProcessesRates, iLevelsStart, iLevelsEnd)    

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
  % without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  % See the GNU Lesser General Public License for more details. 
  % 
  % You should have received a copy of the GNU Lesser General Public License along with this library; 
  % if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
  % 
  %---------------------------------------------------------------------------------------------------------------
  %%==============================================================================================================

  global T0_Vec RatesPath NTint TotDim

  
  filename = strcat(RatesPath,'/T_',num2str(T0_Vec(iT)),'_',num2str(T0_Vec(iT)),'/Rates.csv')
  
  NTsFromFile = 1;
  
  iTToCol(1) = 1;

  delimiter = ',';
  startRow = 6;
  formatSpec = '%f%f%f%*s%*s%*s%*s%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  SystemPESEcutMoleculeNoProcesses = dataArray{:, 1};
  iLevels  = dataArray{:, 1};
  jP       = dataArray{:, 2};
  jLevels  = dataArray{:, 3};
  clearvars delimiter startRow formatSpec fileID dataArray ans;

%   TotDimTemp = [1, TotDim+1];
%   jProcess(:) = TotDimTemp(jP(:)+1)' + jLevels(:) - 1;



  delimiter = ',';
  startRow = 6;
  formatSpec = '%*s%*s%*s';
  for i=1:2*NTsFromFile
    formatSpec = strcat(formatSpec,'%f');
  end
  formatSpec = strcat(formatSpec,'%[^\n\r]');
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  AllRatesTemp = [dataArray{1:end-1}];
  clearvars filename delimiter startRow formatSpec fileID dataArray ans;


  
  for iLine = 1:size(AllRatesTemp,1)
    
    if iLevels(iLine) >= iLevelsStart && iLevels(iLine) <= iLevelsEnd

      if jP(iLine) ~= 0

        RatesMatrix(iLevels(iLine),jLevels(iLine),jP(iLine),iT)      = AllRatesTemp(iLine,iTToCol(iT));
        %RatesSigmaMatrix(iLevels(iLine),jLevels(iLine),jP(iLine),iT) = AllRatesTemp(iLine,iTToCol(iT)+1);  

      else

        ProcessesRates(iLevels(iLine),1,iT)      = AllRatesTemp(iLine,iTToCol(iT));
        %ProcessesRatesSigma(iLevels(iLine),1,iT) = AllRatesTemp(iLine,iTToCol(iT)+1);  
        DissRates(iLevels(iLine),iT)             = AllRatesTemp(iLine,iTToCol(iT));  
        %DissRatesSigma(iLevels(iLine),iT)        = AllRatesTemp(iLine,iTToCol(iT)+1);  

      end

    end
    
  end

  for iP=1:3
    for iLevelss = 1:size(RatesMatrix,1)
      ProcessesRates(iLevelss,iP+1,iT)      = sum(RatesMatrix(iLevelss,:,iP,iT));
      %ProcessesRatesSigma(iLevelss,iP+1,iT) = sqrt(sum(RatesSigmaMatrix(iLevelss,:,iP,iT).^2));
    end
  end
  
  

end