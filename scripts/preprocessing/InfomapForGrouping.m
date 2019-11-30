clear all
close all
clc

%% 1. Run "WriteRatesToNetwork_Infomap.m" file


%% 2. Run Infomap with the network just obtained: infomap ./InelNetwork_Directed.net output/ --directed --preferred-number-of-modules 35 -N 100


%% 3. Convert the Output of Infomap into a csv file for Gephi:

filename = './InelNetwork_Directed.tree';
delimiter = ':';
startRow = 3;
formatSpec = '%f%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
Bin = dataArray{:, 1};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

filename = './InelNetwork_Directed.tree';
delimiter = '"';
startRow = 3;
formatSpec = '%*q%*q%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
Idx = dataArray{:, 1};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

File    = strcat('./Groups.csv');
fileID  = fopen(File,'w');
fprintf(fileID,'id,Group\n');
for i=1:length(Bin)
  fprintf(fileID,'%i,%i\n', Idx(i), Bin(i));
end
fclose(fileID);


%% 4. Open the csv file with the inelastic levels together with the Groups.csv file; clean the binning; export the level and the group in a file "InelBinning.csv"


%% 5. Eliminate emtpy bins:

clc
close all

filename = './InelBinning1.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
Group = dataArray{:, 1};
Level = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

Pointer = zeros(max(Group),1);
for iGroup=1:length(NoRep)
  Pointer(NoRep(iGroup)) = iGroup;
end

FileName = strcat('./NewGroup.csv');
fileID1 = fopen(FileName,'w');
fprintf(fileID1,'id,Bin\n');
for iLevel=1:length(Group)
  FinalGroup(iLevel) = Pointer((Group(iLevel)));
  fprintf(fileID1,'%i,%i\n', Level(iLevel), FinalGroup(iLevel));
end 
fclose(fileID1);


%% 6. Open the csv file with the dissociation levels; export the level and the bin in a file "DissBinning.csv"


%% 7. Merge the 2 Binnings:

filename = './DissLevels.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
id_Diss            = dataArray{:, 1};
Level_Diss         = dataArray{:, 2};
v_Diss             = dataArray{:, 3};
Longitude_Diss     = dataArray{:, 4};
Latitude_Diss      = dataArray{:, 5};
rIn_Diss           = dataArray{:, 6};
EeVVib_Diss        = dataArray{:, 7};
EeVRot_Diss        = dataArray{:, 8};
DeltaEBarrier_Diss = dataArray{:, 9};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

filename = './InelLevels.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
id_Inel            = dataArray{:, 1};
Level_Inel         = dataArray{:, 2};
v_Inel             = dataArray{:, 3};
Longitude_Inel     = dataArray{:, 4};
Latitude_Inel      = dataArray{:, 5};
rIn_Inel           = dataArray{:, 6};
EeVVib_Inel        = dataArray{:, 7};
EeVRot_Inel        = dataArray{:, 8};
DeltaEBarrier_Inel = dataArray{:, 9};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;


filename = './DissBinning.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
Bin_Diss       = dataArray{:, 2};
Level_Diss_Bin = dataArray{:, 1};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

filename = './InelBinning.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
Bin_Inel       = dataArray{:, 2};
Level_Inel_Bin = dataArray{:, 1};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;


Levels  = [Level_Diss; Level_Inel];
NLevels = size(Levels,1);

vqn           = zeros(NLevels,1);
Longitude     = zeros(NLevels,1);
Latitude      = zeros(NLevels,1);
rIn           = zeros(NLevels,1);
EeVVib        = zeros(NLevels,1);
EeVRot        = zeros(NLevels,1);
DeltaEBarrier = zeros(NLevels,1);
Bin           = zeros(NLevels,1);

for i=1:size(id_Inel,1)
    iLevels                = Level_Inel(i);
    vqn(iLevels)           = v_Inel(i);
    Longitude(iLevels)     = Longitude_Inel(i);
    Latitude(iLevels)      = Latitude_Inel(i);
    rIn(iLevels)           = rIn_Inel(i);
    EeVVib(iLevels)        = EeVVib_Inel(i);
    EeVRot(iLevels)        = EeVRot_Inel(i);
    DeltaEBarrier(iLevels) = DeltaEBarrier_Inel(i);
    
    iLevels                = Level_Inel_Bin(i);
    Bin(iLevels)           = Bin_Inel(i);
end

for i=1:size(id_Diss,1)
    iLevels                = Level_Diss(i);
    vqn(iLevels)           = v_Diss(i);
    Longitude(iLevels)     = Longitude_Diss(i);
    Latitude(iLevels)      = Latitude_Diss(i);
    rIn(iLevels)           = rIn_Diss(i);
    EeVVib(iLevels)        = EeVVib_Diss(i);
    EeVRot(iLevels)        = EeVRot_Diss(i);
    DeltaEBarrier(iLevels) = DeltaEBarrier_Diss(i);
    
    iLevels                = Level_Diss_Bin(i);
    Bin(iLevels)           = Bin_Diss(i) + max(Bin_Inel);
end


FileName = strcat('./Levels.csv');
fileID1 = fopen(FileName,'w');
fprintf(fileID1,'id,v,Longitude,Latitude,rIn,EeVVib,EeVRot,DeltaEBarrier,Bin\n');
FileName = strcat('./InelDiss.csv');
fileID2 = fopen(FileName,'w');
fprintf(fileID2,'id,Bin\n');
for i=1:NLevels
    fprintf(fileID1,'%i,%i,%e,%e,%e,%e,%e,%e,%i\n', i, vqn(i), Longitude(i), Latitude(i), rIn(i), EeVVib(i), EeVRot(i), DeltaEBarrier(i), Bin(i));
    fprintf(fileID2,'%i,%i\n', i, Bin(i));
end
fclose(fileID1);
fclose(fileID2);