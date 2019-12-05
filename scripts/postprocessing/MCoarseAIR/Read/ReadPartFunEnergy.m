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

function [QBins, EeV] = ReadPartFunEnergy(iT)     

    %% (METHOD DEPENDENT)

    global T0_Vec BinnedMolName NBinnedMol NBins DatabasePath OutputPath SystemPath KinMthd

    for iBinnedMol=1:NBinnedMol

        % Reading Binned Molecules' Bins Partion Functions and Bins Energies
        filename = strcat(SystemPath,'/',BinnedMolName(iBinnedMol,:),'/',BinnedMolName(iBinnedMol,:),'_',num2str(NBins(iBinnedMol)),'/T',num2str(T0_Vec(iT)),'.dat')
        startRow = 2;
%         formatSpec = '%14f%14f%f%[^\n\r]';
%         fileID = fopen(filename,'r');
%         dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        formatSpec = '%20f%20f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileID);
        TempVec=dataArray{:, 1};
        Temp=length(TempVec);
        QBinsRatio(1:Temp,iBinnedMol) = TempVec;
        QBins(1:Temp,iBinnedMol)      = dataArray{:, 2};
        EeV(1:Temp,iBinnedMol)        = dataArray{:, 3};
        clearvars filename startRow formatSpec fileID dataArray ans;
        
%         filename = strcat(DatabasePath,'/thermo/',BinnedMolName(iBinnedMol,:),'_',num2str(T0_Vec(iT)))
%         delimiter = ' ';
%         formatSpec = '%s%s%[^\n\r]';
%         fileID = fopen(filename,'r');
%         dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
%         fclose(fileID);
%         raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
%         for col=1:length(dataArray)-1
%             raw(1:length(dataArray{col}),col) = dataArray{col};
%         end
%         numericData = NaN(size(dataArray{1},1),size(dataArray,2));
%         for col=[1,2]
%             % Converts text in the input cell array to numbers. Replaced non-numeric
%             % text with NaN.
%             rawData = dataArray{col};
%             for row=1:size(rawData, 1);
%                 % Create a regular expression to detect and remove non-numeric prefixes and
%                 % suffixes.
%                 regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
%                 try
%                     result = regexp(rawData{row}, regexstr, 'names');
%                     numbers = result.numbers;
%     
%                     % Detected commas in non-thousand locations.
%                     invalidThousandsSeparator = false;
%                     if any(numbers==',');
%                         thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
%                         if isempty(regexp(numbers, thousandsRegExp, 'once'));
%                             numbers = NaN;
%                             invalidThousandsSeparator = true;
%                         end
%                     end
%                     % Convert numeric text to numbers.
%                     if ~invalidThousandsSeparator;
%                         numbers = textscan(strrep(numbers, ',', ''), '%f');
%                         numericData(row, col) = numbers{1};
%                         raw{row, col} = numbers{1};
%                     end
%                 catch me
%                 end
%             end
%         end
%         I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),2); % Find rows with non-numeric cells
%         raw(I,:) = [];
%         Matrix = cell2mat(raw);
%         clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me I J K;
%         temp=size(Matrix,1);
%         QBins(1:temp,iBinnedMol) = Matrix(:,1);
%         EeV(1:temp,iBinnedMol)   = Matrix(:,2);
%         clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me I J K;
%           
        
%         if KinMthd(iBinnedMol,:) == 'CGM'
%           filename = strcat(SystemPath,'/',BinnedMolName(iBinnedMol,:),'/',BinnedMolName(iBinnedMol,:),'_',num2str(NBins(iBinnedMol)),'/E_low.dat');
%           delimiter = ' ';
%           startRow = 2;
%           formatSpec = '%*s%f%*s%*s%[^\n\r]';
%           fileID = fopen(filename,'r');
%           dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%           fclose(fileID);
%           EeV(1:Temp,iBinnedMol) = dataArray{:, 1};
%           clearvars filename delimiter startRow formatSpec fileID dataArray ans;
%         end      
        
    end
    
end
