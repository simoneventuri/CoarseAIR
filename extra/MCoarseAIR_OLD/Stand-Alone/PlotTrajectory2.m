close all
clear all
clc


filename = '/Users/sventuri/WORKSPACE/CoarseAIR/run_CHN/Test/T_7000_7000/Bins_100_0/Node_1/Proc_1/XXEvo-54.out';
delimiter = ',';
startRow = 3;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
Times = dataArray{:, 1};
x1 = dataArray{:, 2};
y1 = dataArray{:, 3};
z1 = dataArray{:, 4};
x2 = dataArray{:, 5};
y2 = dataArray{:, 6};
z2 = dataArray{:, 7};
x3 = dataArray{:, 8};
y3 = dataArray{:, 9};
z3 = dataArray{:, 10};
vx1 = dataArray{:, 11};
vy1 = dataArray{:, 12};
vz1 = dataArray{:, 13};
vx2 = dataArray{:, 14};
vy2 = dataArray{:, 15};
vz2 = dataArray{:, 16};
vx3 = dataArray{:, 17};
vy3 = dataArray{:, 18};
vz3 = dataArray{:, 19};
iP = dataArray{:, 20};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;


figure
plot3(x1,y1,z1)
hold on
plot3(x2,y2,z2)
plot3(x3,y3,z3)
xlim([-20,20])
ylim([-20,20])
zlim([-20,20])


iP