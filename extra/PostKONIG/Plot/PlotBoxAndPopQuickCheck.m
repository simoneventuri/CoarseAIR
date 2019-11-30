%close all
clc 
clear all

NBins        = 9390
OutputFolder = '/Users/sventuri/Desktop/TempTemp/'
%FileN       = '_NASA_ORIG'
%FileN       = '_NASA'
FileN       = '_NASA_NEW'
%FileN       = '_NN_40000'
%FileN        = '_NN_5000'
%FileN        = '_NN_5000_Plus'
%FileN        = '_NN_5000_Minus'
%FileN        = '_NN_1000'
%FileN       = '_LEPS'
%FileN        = '_NoTri'
%FileN       = '_NASA_CG'
iStepsVec    = [1, 100, 200, 300, 400, 500, 600]

filename = strcat(OutputFolder,'/box',FileN,'.dat');
formatSpec = '%20f%20f%20f%20f%20f%20f%20f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
time = dataArray{:, 1};
XN = dataArray{:, 2};
XN2 = dataArray{:, 3};
T = dataArray{:, 4};
rho = dataArray{:, 5};
P = dataArray{:, 6};
nd = dataArray{:, 7};
E = dataArray{:, 8};
clearvars filename formatSpec fileID dataArray ans;
  
figure(1)
semilogx(time,XN)
hold on
xlim([1.e-8, 1.e-3])

figure(2)
semilogx(time,E)
hold on
xlim([1.e-8, 1.e-3])


iFigure = 2

filename = strcat(OutputFolder,'/pop_N2',FileN,'.dat');
delimiter = ' ';
formatSpec = '%f%f%*s%*s%[^\n\r]';

for iStepss = iStepsVec

  iFigure = iFigure + 1;
  figure(iFigure)
  endRow = 1;
  for iSteps = 1:iStepss
    startRow = endRow + 2
    endRow   = startRow + NBins - 1

    if iSteps==iStepss
      fileID = fopen(filename,'r');
      dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
      fclose(fileID);
      EeV = dataArray{:, 1};
      Pop = dataArray{:, 2};
      clearvars dataArray ans;

      semilogy(EeV,Pop,'o')
      hold on
      ylim([1.e5, 1.e22])

    end

  end 
  
end
