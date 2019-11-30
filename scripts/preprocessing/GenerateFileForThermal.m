close all
clear all
clc

NLevels = 7010%[ 217; 6997; 480]

File   = strcat('./InelDiss.csv');
fileID = fopen(File,'w');
fprintf(fileID,'id,bin\n')
for iLevels=1:NLevels
  fprintf(fileID,'%i,%i\n', iLevels,1);
end
fclose(fileID);