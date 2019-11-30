close all
clc
clear all


filename = './BinIndex-150Bound_50Quasi.dat';
delimiter = '\t';
formatSpec = '%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
Id    = dataArray{:, 1};
Temp1 = dataArray{:, 2};
Temp2 = dataArray{:, 3};
Bin   = dataArray{:, 4};
clearvars filename delimiter formatSpec fileID dataArray ans;


FileName1 = strcat('./InelDiss.dat');
fileID1   = fopen(FileName1,'w');
for i = 1:length(Bin)
  fprintf(fileID1,'%i,%i\n', Id(i)+1, Bin(i)+1);
end
fclose(fileID1);