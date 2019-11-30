close all
clear all
clc


filenameOrig = '/Users/sventuri/WORKSPACE/CoarseAIR/coarseair/dtb/CO2/O2/O2_levels_venturi.dat';
filenameNew = '/Users/sventuri/WORKSPACE/CoarseAIR/coarseair/dtb/CO2/O2/O2_levels_venturi_new.dat';


delimiter = ' ';
startRow = 16;
formatSpec = '%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filenameOrig,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
vOrig = dataArray{:, 1};
jOrig = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
NLevelsOrig = length(vOrig);

delimiter = ' ';
startRow = 16;
formatSpec = '%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filenameNew,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
vNew = dataArray{:, 1};
jNew = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
NLevelsNew = length(vNew);

figure
plot(vOrig,jOrig,'ro')
hold on
plot(vNew,jNew,'go')


OrigToNew = 0 * vOrig;
for iLevels = 1:NLevelsOrig
  
  for jLevels = 1:NLevelsNew
    
    if vNew(jLevels) == vOrig(iLevels) && jNew(jLevels) == jOrig(iLevels)
      OrigToNew(iLevels) = jLevels;
    end
    
  end
  
end 

fileID = fopen('~/Desktop/LevelsMapping_O2.dat','w');
fprintf(fileID,'%d\n',OrigToNew);