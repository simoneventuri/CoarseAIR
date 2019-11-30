clear all
close all
clc

figure
for iPES = 1:9
  filename = strcat('/Users/sventuri/WORKSPACE/CoarseAIR/run_O3/Test/PlotPES/PES_',num2str(iPES),'/PESDistVsAngles.csv')
  delimiter = ',';
  startRow = 2;
  formatSpec = '%f%f%f%f%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  r1 = dataArray{:, 1};
  r2 = dataArray{:, 2};
  r3 = dataArray{:, 3};
  Angle = dataArray{:, 4};
  E = dataArray{:, 5};
  clearvars filename delimiter startRow formatSpec fileID dataArray ans;
  
  plot(Angle,E)
  hold on
end