% This Program converts the Energy Levels List from Minnesota's format to the CoarseAIR's (NASA's) one

clear all;
close all;
clc;


% Initially, 2 Files have been manually created, by separating info about J=0 Levels from the rest of data;

File1   = '/home/venturi/WORKSPACE/CoarseAIR/coarseair/dtb/Molecules/N2/UMN_ForN4/Original_FromRoss_1.inp';
File2   = '/home/venturi/WORKSPACE/CoarseAIR/coarseair/dtb/Molecules/N2/UMN_ForN4/Original_FromRoss_2.inp';

NewFile = '/home/venturi/WORKSPACE/CoarseAIR/coarseair/dtb/Molecules/N2/UMN_ForN4/FromRoss.inp';

fileID = fopen(File1,'r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
fclose(fileID);
matA = vec2mat(A,7);

fileID = fopen(File2,'r');
formatSpec = '%f';
B = fscanf(fileID,formatSpec);
fclose(fileID);
matB = vec2mat(B,7);
NLevels = size(matB,1)

fileID = fopen(NewFile,'w');
fmt = '%5d %5d %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E\n';
fprintf(fileID, '######################################################################################################################################################\n');
fprintf(fileID, '# jqn   : the rotational q.n. of the i''th quantum state\n');
fprintf(fileID, '# vqn   : the vibrational q.n. of the i''th quantum state\n');
fprintf(fileID, '# eint  : internal energy of i''th quantum state [Eh]\n');
fprintf(fileID, '# egam  : Half width of i''th quantum state\n');
fprintf(fileID, '# rmin  : the position of the potential minimum (included centrifugal potential) for i''th quantum state\n');
fprintf(fileID, '# vmin  : the value of the potential minimun (inc. cent. pot.)\n');
fprintf(fileID, '# vmax  : the value of the local potential maximum (inc. cent. pot.)\n');
fprintf(fileID, '# tau   : the vibrational period of the i''th quantum state\n');
fprintf(fileID, '# ri    : inner turning point\n');
fprintf(fileID, '# ro    : outter turning point\n');
fprintf(fileID, '# rmax  : location of maximum in centrifugal barrier\n');
fprintf(fileID, '######################################################################################################################################################\n');
fprintf(fileID, '#    vqn  jqn   eint            egam            rmin          rmax          vmin            vmax            tau           ri             ro\n');
fprintf(fileID, '######################################################################################################################################################\n');

shifter = matA(1,7)
shifter = 0;
for i = 1:NLevels
  C(i,1)= matB(i,3);
  C(i,2)= matB(i,2);
  C(i,3)= matB(i,4) - shifter;
  k     = matB(i,2);
  C(i,4)= 0;
  C(i,5)= matA(k+1,4);
  C(i,6)= matA(k+1,6);
  C(i,7)= matA(k+1,5) - shifter;
  C(i,8)= matA(k+1,7) - shifter;
  C(i,9)= matB(i,7);
  C(i,10)= matB(i,5);
  C(i,11)= matB(i,6);
end


for i = 1:NLevels
    fprintf(fileID,fmt, C(i,:));
end
fclose(fileID);
