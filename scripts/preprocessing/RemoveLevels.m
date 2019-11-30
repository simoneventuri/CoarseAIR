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
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU Lesser General Public License for more details. 
% 
% You should have received a copy of the GNU Lesser General Public License along with this library; 
% if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
% 
%---------------------------------------------------------------------------------------------------------------
%%==============================================================================================================

close all
clear all 
clc


FileOrig = '/Users/sventuri/WORKSPACE/CoarseAIR/coarseair/dtb/CO2/O2/O2_levels_venturi_new.dat';
FileNew  = '/Users/sventuri/WORKSPACE/CoarseAIR/coarseair/dtb/CO2/O2/O2_levels_venturi_Final.dat';


startRow = 16;
formatSpec = '%6s%5s%15s%15s%15s%15s%15s%15s%15s%15s%15s%[^\n\r]';
fileID = fopen(FileOrig,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
  raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,4,5,6,7,8,9,10,11]
  % Converts text in the input cell array to numbers. Replaced non-numeric
  % text with NaN.
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
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells
vqn  = cell2mat(raw(:, 1));
jqn  = cell2mat(raw(:, 2));
Var1 = cell2mat(raw(:, 3));
Var2 = cell2mat(raw(:, 4));
Var3 = cell2mat(raw(:, 5));
Var4 = cell2mat(raw(:, 6));
Var5 = cell2mat(raw(:, 7));
Var6 = cell2mat(raw(:, 8));
Var7 = cell2mat(raw(:, 9));
Var8 = cell2mat(raw(:, 10));
Var9 = cell2mat(raw(:, 11));
clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;
NLevels = length(Var1)


filename = '/Users/sventuri/Desktop/ToRemove.dat';
delimiter = '';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
kLevels = dataArray{:, 1};
clearvars filename delimiter formatSpec fileID dataArray ans;


fileID = fopen(FileNew,'w');
lLevels = 1;
for iLevels = 1:NLevels
  flagg = 1;
  for jLines = 1:length(kLevels)
    if iLevels == kLevels(jLines)
      flagg = 0;
    end
  end
  if flagg == 1
    fprintf(fileID,'%6d %4d %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E\n',vqn(lLevels), jqn(lLevels), Var1(lLevels), Var2(lLevels), Var3(lLevels), Var4(lLevels), Var5(lLevels), Var6(lLevels), Var7(lLevels), Var8(lLevels), Var9(lLevels) );
    lLevels = lLevels+1;
  end
end
fclose(fileID);
