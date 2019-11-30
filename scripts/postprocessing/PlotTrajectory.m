close all
clear all
clc

% filename = '/Users/sventuri/WORKSPACE/CoarseAIR/run_CO2_StS/Test/T_5000_5000/Bins_7000_0/Node_1/Proc_1/XXEvo-8.out';
% %filename = '/Users/sventuri/WORKSPACE/CoarseAIR/run_N3_NASA/Test/T_10000_10000/Bins_1_0/Node_1/Proc_1/XXEvo-1.out';
% delimiter = ',';
% startRow = 3;
% formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% Times = dataArray{:, 1};
% x1  = dataArray{:, 2};
% y1  = dataArray{:, 3};
% z1  = dataArray{:, 4};
% x2  = dataArray{:, 5};
% y2  = dataArray{:, 6};
% z2  = dataArray{:, 7};
% x3  = dataArray{:, 8};
% y3  = dataArray{:, 9};
% z3  = dataArray{:, 10};
% vx1 = dataArray{:, 11};
% vy1 = dataArray{:, 12};
% vz1 = dataArray{:, 13};
% vx2 = dataArray{:, 14};
% vy2 = dataArray{:, 15};
% vz2 = dataArray{:, 16};
% vx3 = dataArray{:, 17};
% vy3 = dataArray{:, 18};
% vz3 = dataArray{:, 19};
% iP  = dataArray{:, 20};
% clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%filename = '/Users/sventuri/WORKSPACE/CoarseAIR/run_CO2_StS/Test/T_5000_5000/Bins_7000_0/Node_1/Proc_1/XXEvo-8.out';
%filename = '/Users/sventuri/WORKSPACE/CoarseAIR/run_N3_NASA/Test/T_10000_10000/Bins_1_0/Node_1/Proc_1/XXEvo-1.out';
filename = '/Users/sventuri/Downloads/Proc_1/XXEvo-5.out';
delimiter = ',';
startRow = 3;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
Times = dataArray{:, 1};
x1  = dataArray{:, 2};
y1  = dataArray{:, 3};
z1  = dataArray{:, 4};
x2  = dataArray{:, 5};
y2  = dataArray{:, 6};
z2  = dataArray{:, 7};
x3  = dataArray{:, 8};
y3  = dataArray{:, 9};
z3  = dataArray{:, 10};
x4  = dataArray{:, 11};
y4  = dataArray{:, 12};
z4  = dataArray{:, 13};
vx1 = dataArray{:, 14};
vy1 = dataArray{:, 15};
vz1 = dataArray{:, 16};
vx2 = dataArray{:, 17};
vy2 = dataArray{:, 18};
vz2 = dataArray{:, 19};
vx3 = dataArray{:, 20};
vy3 = dataArray{:, 21};
vz3 = dataArray{:, 22};
vx4 = dataArray{:, 23};
vy4 = dataArray{:, 24};
vz4 = dataArray{:, 25};
iP  = dataArray{:, 26};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

figure
scatter3(x1,y1,z1)
hold on
scatter3(x2,y2,z2)
scatter3(x3,y3,z3)
scatter3(x4,y4,z4)

iP


% figure
% scatter3(x1(2:100),y1(2:100),z1(2:100))
% hold on
% scatter3(x2(2:100),y2(2:100),z2(2:100))
% scatter3(x3(2:100),y3(2:100),z3(2:100))

