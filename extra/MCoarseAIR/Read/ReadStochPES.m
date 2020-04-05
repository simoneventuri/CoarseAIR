%% The Function reads the values at the Stochastic PES surface at grid points
%

function [StochPESR, StochPESAng, StochPESEeV] = ReadStochPES()    

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

  global NBins T0_Vec RatesPath ProcToLevIP ReadAllRatesFlg ProduceMatFlg
  
  
  filename = strcat('../Test/PlotPES/GridForStochPESInfo.dat')
  delimiter = ' ';
  formatSpec = '%s%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
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
  I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),2); % Find rows with non-numeric cells
  raw(I,:) = [];
  Temp = cell2mat(raw(:, 1));
  
  NPESSamples = Temp(1);
  RMin        = Temp(2);
  RMax        = Temp(3);
  NPoints     = Temp(4);
  hGrid       = Temp(5);
  NAng        = Temp(6);
  StochPESAng = Temp(7:end);
  clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me I J K Temp;

  
  StochPESR = [];
  for iA=1:NAng
    ii=1;
    for i=1:NPoints
      for j=1:i
        Rp(1) = RMin + hGrid * (i-1);
        Rp(3) = RMin + hGrid * (j-1);
        Rp(2) = sqrt( Rp(1)^2 + Rp(3)^2 - 2.d0*Rp(1)*Rp(3) * cos(StochPESAng(iA)/180.d0*pi) );
        StochPESR(iA,ii,:) = [Rp(1), Rp(2), Rp(3)];
        ii = ii+1;
      end
    end
  end
  NTotPoints = ii-1
  
  
  StochPESEeV = zeros(NAng, NTotPoints, NPESSamples );
  for iA=1:NAng
 
    filename = strcat('../Test/PlotPES/StochPESFromGrid.csv.', num2str(StochPESAng(iA)))
    delimiter = ',';
    startRow = 4;
    formatSpec = '%q%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
      raw(1:length(dataArray{col}),col) = dataArray{col};
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));
    rawData = dataArray{1};
    for row=1:size(rawData, 1);
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
    I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),2); % Find rows with non-numeric cells
    raw(I,:) = [];
    Temp = cell2mat(raw(:, 1));
    
    ii=1;
    jj=1;
    for i=1:NTotPoints
      %Temp(ii:ii+NPESSamples-1)
      StochPESEeV(iA,i,1:NPESSamples) = Temp(ii:ii+NPESSamples-1);
      ii = ii+NPESSamples+1;
    end
    clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me I J K Temp;
    
  end
  
  if ProduceMatFlg == 1
      
    filenameRates = strcat('./StochPES')
    save(filenameRates,'StochPESR', 'StochPESAng', 'StochPESEeV','-v7.3');

  end
  
  
end