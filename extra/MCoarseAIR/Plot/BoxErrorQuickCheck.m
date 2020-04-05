close all
clc
clear all

filename = '/Users/sventuri/WORKSPACE/CG-QCT/run_N3_NASA_Orig/Test/RunKonig/output/T_10000/output/box.dat';
formatSpec = '%20f%20f%20f%20f%20f%20f%20f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
time_NASA = dataArray{:, 1};
XN_NASA   = dataArray{:, 2};
XN2_NASA  = dataArray{:, 3};
clearvars filename formatSpec fileID dataArray ans;


filename = '/Users/sventuri/WORKSPACE/CG-QCT/run_N3_CG_NASA_40/Test/RunKonig/output/T_10000/output/box.dat';
formatSpec = '%20f%20f%20f%20f%20f%20f%20f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
time_CG_NASA = dataArray{:, 1};
XN_CG_NASA   = dataArray{:, 2};
XN2_CG_NASA  = dataArray{:, 3};
clearvars filename formatSpec fileID dataArray ans;


filename = '/Users/sventuri/WORKSPACE/CG-QCT/run_N3_CG_NN_40/Test/RunKonig/output/T_10000/output/box.dat';
formatSpec = '%20f%20f%20f%20f%20f%20f%20f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
time_CG_NN = dataArray{:, 1};
XN_CG_NN   = dataArray{:, 2};
XN2_CG_NN  = dataArray{:, 3};
clearvars filename formatSpec fileID dataArray ans;


Error_CG = (XN_CG_NASA-XN_NASA)./XN_NASA;
figure
semilogx(time_NASA,Error_CG)


Error_CG_NN = (XN_CG_NN-XN_NASA)./XN_NASA;
figure
semilogx(time_NASA,Error_CG_NN)