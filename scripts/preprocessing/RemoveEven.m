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


FileOrig = '/Users/sventuri/WORKSPACE/CoarseAIR/dtb/Molecules/O2/UMN/FromUMN_Sorted.inp';
FileNew  = '/Users/sventuri/WORKSPACE/CoarseAIR/dtb/Molecules/O2/UMN/FromUMN_OnlyOdd.inp';


filename = FileOrig;
startRow = 16;
formatSpec = '%6f%5f%15s%15C%15C%15C%15C%15C%15s%15s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray{3} = strtrim(dataArray{3});
dataArray{9} = strtrim(dataArray{9});
dataArray{10} = strtrim(dataArray{10});
dataArray{11} = strtrim(dataArray{11});
fclose(fileID);
vv = dataArray{:, 1};
jj = dataArray{:, 2};
x1 = dataArray{:, 3};
x2 = dataArray{:, 4};
x3 = dataArray{:, 5};
x4 = dataArray{:, 6};
x5 = dataArray{:, 7};
x6 = dataArray{:, 8};
x7 = dataArray{:, 9};
x8 = dataArray{:, 10};
x9 = dataArray{:, 11};
clearvars filename startRow formatSpec fileID dataArray ans;
NLevels = length(vv)



fileID = fopen(FileNew,'w');
for i=1:NLevels
    if mod(jj(i),2) ~= 0
        fprintf(fileID,'%6d %4d %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E\n', vv(i), jj(i), x1(i), x2(i), x3(i), x4(i), x5(i), x6(i), x7(i), x8(i), x9(i));
    end
end
fclose(fileID);