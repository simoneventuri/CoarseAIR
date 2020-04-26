%% The Function reads the rates (K(i,j)) and computes the Processes (Dissociation, Exchange 1,2 and 3) overall rates (K(i,Process))
%
%  Input Arguments:  - iLevelsStart:  (Optional) Initial Level/Bin to read. If not given, it is assumed equal to be equal to 1
%                    - iLevelsEnd:    (Optional) Final Level/Bin to read. If not given, it is assumed equal to be equal to NBins
%
%  Input Global Var: - T0_Vec:        Vector of Translational Temperatures (e.g.: [10000])
%                    - RatesPath:     The path to the output folder (e.g.: ../Test/N3/N2/Rates)
%

%function [RatesSigma, RatesMatrix, RatesSigmaMatrix, DissRates, DissRatesSigma, ProcessesRates, ProcessesRatesSigma] = ReadRates(iT, RatesSigma, RatesMatrix, RatesSigmaMatrix, DissRates, DissRatesSigma, ProcessesRates, ProcessesRatesSigma, iBinsStart, iBinsEnd)    
function [RatesMatrix, DissRates, ProcessesRates] = ReadRatesAllTint(RatesMatrix, DissRates, ProcessesRates, iLevelsStart, iLevelsEnd)    

  %%==============================================================================================================
  % 
  % Coarse-Grained method for Quasi-Classical Trajectories (CG-QCT) 
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

  filename = strcat(RatesPath,'/AllRates.csv')

  
  startRow = 4;
  endRow = 4;
  formatSpec = '%*1s%10s%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
  for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
  end
  numericData = NaN(size(dataArray{1},1),size(dataArray,2));
  rawData = dataArray{1};
  for row=1:size(rawData, 1);
    % Create a regular expression to detect and remove non-numeric prefixes and
    % suffixes.
    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
    try
      result = regexp(rawData{row}, regexstr, 'names');
      numbers = result.numbers;

      % Detected commas in non-thousand locations.
      invalidThousandsSeparator = false;
      if any(numbers==',');
        thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
        if isempty(regexp(numbers, thousandsRegExp, 'once'));
          numbers = NaN;
          invalidThousandsSeparator = true;
        end
      end
      % Convert numeric text to numbers.
      if ~invalidThousandsSeparator;
        numbers = textscan(strrep(numbers, ',', ''), '%f');
        numericData(row, 1) = numbers{1};
        raw{row, 1} = numbers{1};
      end
    catch me
    end
  end
  NTsFromFile = cell2mat(raw);
  clearvars startRow endRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me;

  
  
  startRow = 4;
  endRow = 4;
  formatSpec = '%*29s%20s%20s%20s%s%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
  for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
  end
  numericData = NaN(size(dataArray{1},1),size(dataArray,2));
  for col=[1,2,3,4]
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
      % Create a regular expression to detect and remove non-numeric prefixes and
      % suffixes.
      regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
      try
        result = regexp(rawData{row}, regexstr, 'names');
        numbers = result.numbers;

        % Detected commas in non-thousand locations.
        invalidThousandsSeparator = false;
        if any(numbers==',');
          thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
          if isempty(regexp(numbers, thousandsRegExp, 'once'));
            numbers = NaN;
            invalidThousandsSeparator = true;
          end
        end
        % Convert numeric text to numbers.
        if ~invalidThousandsSeparator;
          numbers = textscan(strrep(numbers, ',', ''), '%f');
          numericData(row, col) = numbers{1};
          raw{row, col} = numbers{1};
        end
      catch me
      end
    end
  end
  I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),1); % Find columns with non-numeric cells
  raw(:,I) = [];
  columnIndices = cumsum(~I);
  TsFromFile = cell2mat(raw);
  clearvars startRow endRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me I columnIndices;

  iTToCol = zeros(NTint,1);
  for iT = 1:NTint
    for jT = 1:length(TsFromFile)
      if T0_Vec(iT) == TsFromFile(jT)
        iTToCol(iT) = jT*2 - 1;
      end
    end
  end


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

        RatesMatrix(iLevels(iLine),jLevels(iLine),jP(iLine),:)      = AllRatesTemp(iLine,iTToCol(:));
        %RatesSigmaMatrix(iLevels(iLine),jLevels(iLine),jP(iLine),:) = AllRatesTemp(iLine,iTToCol(:)+1);  

      else

        ProcessesRates(iLevels(iLine),1,:)      = AllRatesTemp(iLine,iTToCol(:));
        %ProcessesRatesSigma(iLevels(iLine),1,:) = AllRatesTemp(iLine,iTToCol(:)+1);  
        DissRates(iLevels(iLine),:)             = AllRatesTemp(iLine,iTToCol(:));  
        %DissRatesSigma(iLevels(iLine),:)        = AllRatesTemp(iLine,iTToCol(:)+1);  

      end

    end
    
  end

  for iP=1:3
    for iT = 1:NTint      
      for iLevelss = 1:size(RatesMatrix,1)
        ProcessesRates(iLevelss,iP+1,iT)      = sum(RatesMatrix(iLevelss,:,iP,iT));
        %ProcessesRatesSigma(iLevelss,iP+1,iT) = sqrt(sum(RatesSigmaMatrix(iLevelss,:,iP,iT).^2));
      end
    end
  end
  
  

end