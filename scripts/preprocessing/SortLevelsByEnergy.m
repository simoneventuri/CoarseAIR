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


FileOrig = '/home/venturi/WORKSPACE/output.dat';
FileNew  = '/home/venturi/WORKSPACE/CoarseAIR/coarseair/dtb/Molecules/N2/LeRoy/MyLeroy_FromRobyn_Cut.inp';


startRow = 16;
formatSpec = '%6s%5s%15s%15s%15s%15s%15s%15s%15s%15s%s%[^\n\r]';
fileID = fopen(FileOrig,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,4,5,6,7,8,9,10,11]
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
R = cellfun(@(x) (~isnumeric(x) && ~islogical(x)) || isnan(x),raw); % Find non-numeric cells
raw(R) = {1.0E-99}; % Replace non-numeric cells
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
clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;NLevels = length(Var1)

aa=(Var2  > 0.d0) + (Var2 < 1.d-100);
Var2(aa==2) = 1.d-103;  

%[En, ToNew] = sort(Var1);
ToNew       = [1:length(vqn)];

fileID = fopen(FileNew,'w');
for iLevels = 1:NLevels
  if Var2(ToNew(iLevels)) > 0.d0
    if Var2(ToNew(iLevels)) <= 1.d-200
      Var2(ToNew(iLevels)) = 0.d0;
    elseif Var2(ToNew(iLevels)) <= 1.d-100
      Var2(ToNew(iLevels)) = 1.0000000E-99;
    end
  end
  fprintf(fileID,'%6d %4d %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E\n',vqn(ToNew(iLevels)), jqn(ToNew(iLevels)), Var1(ToNew(iLevels)), Var2(ToNew(iLevels)), Var3(ToNew(iLevels)), Var4(ToNew(iLevels)), Var5(ToNew(iLevels)), Var6(ToNew(iLevels)), Var7(ToNew(iLevels)), Var8(ToNew(iLevels)), Var9(ToNew(iLevels)) );
end
fclose(fileID);
