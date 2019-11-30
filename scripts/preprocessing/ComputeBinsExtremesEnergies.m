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
EH_To_EeV = 27.21138397127;


System     = 'CO2'
Molecule   = 'CO'
LevelsFile = 'CO_levels_venturi.dat'
NBins      =  60
dtbPath    = '/Users/sventuri/WORKSPACE/QCT/qct/dtb/'

NLevelsBased = 1
EnergyBased  = 0


filename = strcat(dtbPath,System,'/',Molecule,'/',LevelsFile);
startRow = 16;
formatSpec = '%*11s%15s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
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
EH = cell2mat(raw);
clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me I J K;

EeV=EH.*EH_To_EeV;
min(EeV)
max(EeV)
NLevels=length(EeV)
NLevelsBin=floor(NLevels/NBins)
[EeVSorted,Crrsp] = sort(EeV,1,'ascend');

if NLevelsBased == 1
    
    FileName = strcat('./E_low_',Molecule,'.dat');
    fileID = fopen(FileName,'w');
    fprintf(fileID,'%18s\n','# Bin          Energy [eV]');
    iLevels=1;
    for iBins=1:NBins
        iLevelsInBin=1;
        BinEeV(iBins)=EeVSorted(iLevels);
        A = [iBins; BinEeV(iBins)];
        fprintf(fileID,'%5i %20.10f\n',A);
        while iLevelsInBin <= NLevelsBin
            iLevels=iLevels+1;
            iLevelsInBin=iLevelsInBin+1;
        end
    end
    fclose(fileID);
    
elseif EnergyBased == 1
    
    FileName = strcat('./E_low_',Molecule,'.dat');
    fileID = fopen(FileName,'w');
    fprintf(fileID,'%18s\n','# Bin          Energy [eV]');
    iLevels=1;
    DeltaE = (max(EeVSorted)-min(EeVSorted)) / (NBins+1);
    for iBins=1:NBins
        BinEeV(iBins) = min(EeVSorted) + (iBins-1)*DeltaE;
        A = [iBins; BinEeV(iBins)];
        fprintf(fileID,'%5i %20.10f\n',A);
    end
    fclose(fileID);
    
end
