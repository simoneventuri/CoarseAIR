clc

iBinnedMol = 1;
ExpVec(1:NLevels(iBinnedMol),1) = Levelg(1:NLevels(iBinnedMol),1) .* exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0_Vec(1) .* UKb) );


filename = '/home/venturi/WORKSPACE/CoarseAIR/RESULTS/From_NASA/N3_10000_plato/N3_Diss.dat';
delimiter = {',','(',')',':'};
formatSpec = '%*s%*s%*s%f%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
TempRates = dataArray{:, 1};
clearvars filename delimiter formatSpec fileID dataArray ans;

filename = '/home/venturi/WORKSPACE/CoarseAIR/RESULTS/From_NASA/N3_10000_plato/N3_Diss.dat';
delimiter = {',','(',')',':'};
formatSpec = '%*s%f%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
TempLevels = dataArray{:, 1};
clearvars filename delimiter formatSpec fileID dataArray ans;

ProcessesRates(TempLevels,1,1) = TempRates;


filename = '/home/venturi/WORKSPACE/CoarseAIR/RESULTS/From_NASA/RatesMarco_opt2_new/FullSet.ascii.dat';
formatSpec = '%8f%8f%*26*s%*26*s%*26*s%*26f%26f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
Source    = dataArray{:, 1};
Target    = dataArray{:, 2};
TempRates = dataArray{:, 3};
clearvars filename formatSpec fileID dataArray ans;

for i=1:length(Source)
  TempRatesInverse(i)                  = TempRates(i) .* ExpVec(Source(i)) ./ ExpVec(Target(i));
  RatesMatrix(Source(i),Target(i),1,1) = TempRates(i);
  RatesMatrix(Target(i),Source(i),1,1) = TempRatesInverse(i);
end