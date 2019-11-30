close all
clc

filename = '/Users/sventuri/Desktop/SSSS/Inelastic/ForceAtlas2/G_20000.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
id         = dataArray{:, 1};
v          = dataArray{:, 2};
Longitude  = dataArray{:, 3};
Latitude   = dataArray{:, 4};
rIn        = dataArray{:, 5};
EeVVib     = dataArray{:, 6};
EeVRot     = dataArray{:, 7};
x1         = dataArray{:, 8};
y1         = dataArray{:, 9};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

x1      = x1 - x1(1);
y1      = y1 - y1(1);
alp1    = - atan(y1(end) ./ x1(end));
x1New   = x1 .* cos(alp1) - y1 .* sin(alp1);
y1New   = x1 .* sin(alp1) + y1 .* cos(alp1);
x1      = x1New;
y1      = y1New;
% len1    = sqrt(x1(end).^2 + y1(end).^2);
% x1      = x1 ./ len1;
% y1      = y1 ./ len1;



filename = '/Users/sventuri/Desktop/SSSS/InelExchange/ForceAtlas2/G_16000.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
x2         = dataArray{:, 8};
y2         = dataArray{:, 9};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

x2      = x2 - x2(1);
y2      = y2 - y2(1);
alp2    = -atan(y2(end) ./ x2(end));
x2New   = x2 .* cos(alp2) - y2 .* sin(alp2);
y2New   = x2 .* sin(alp2) + y2 .* cos(alp2);
x2      = - x2New;
y2      = y2New;
% len2    = sqrt(x2(end).^2 + y2(end).^2);
% x2      = x2 ./ len2;
% y2      = y2 ./ len2;


figure
scatter(x1,y1)
hold on
scatter(x2,y2)



FileName = strcat('/Users/sventuri/Desktop/SSSS/Inelastic/ForceAtlas2/G_20000_Bis.csv');
fileID = fopen(FileName,'w');
fprintf(fileID,'id,v,Longitude,Latitude,rIn,EeVVib,EeVRot,x,y\n');
for i = 1:max(id)
  fprintf(fileID,'%i,%i,%e,%e,%e,%e,%e,%e,%e\n', i, v(i), -Longitude(i), Latitude(i), rIn(i), EeVVib(i), EeVRot(i), y1(i), x1(i));
end
fclose(fileID);

FileName = strcat('/Users/sventuri/Desktop/SSSS/InelExchange/ForceAtlas2/G_16000_Bis.csv');
fileID = fopen(FileName,'w');
fprintf(fileID,'id,v,Longitude,Latitude,rIn,EeVVib,EeVRot,x,y\n');
for i = 1:max(id)
  fprintf(fileID,'%i,%i,%e,%e,%e,%e,%e,%e,%e\n', i, v(i), -Longitude(i), Latitude(i), rIn(i), EeVVib(i), EeVRot(i), y2(i), x2(i));
end
fclose(fileID);