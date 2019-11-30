close all
clear all
clc

T      = 10000.d0;
gC     = 1.d0;
gO     = 5.d0;
CMass  = 21868.661757d0;
OMass  = 29148.94559d0;
COMass = CMass + OMass;

EEh_to_EeV = 27.2114d0;
UKb        = 1.380658e-23;
Ue         = 1.602191e-19;
DSWtoKg    = 1.d-3/1.8208e+03;
ATMToPa    = 101325.d0;
Plnck      = 6.62607004d-34;
AvN        = 6.0221409e+23;

LevelsInfoFile   = '/Users/sventuri/Desktop/CO2_for_Amal/CO_LevelsInfo.dat'
InelRatesFile    = '/Users/sventuri/Desktop/CO2_for_Amal/CO_InelRates.csv'
NEWInelRatesFile = '/Users/sventuri/Desktop/CO2_for_Amal/CO_InelRates_New.csv'
DissRatesFile    = '/Users/sventuri/Desktop/CO2_for_Amal/CO_DissRates.csv'
NEWDissRatesFile = '/Users/sventuri/Desktop/CO2_for_Amal/CO_DissRates_New.csv'


filename = LevelsInfoFile;
startRow = 2;
formatSpec = '%8f%7f%7f%15f%15f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
Level = dataArray{:, 1};
vqn   = dataArray{:, 2};
jqn   = dataArray{:, 3};
EEh   = dataArray{:, 4};
g     = dataArray{:, 5};
iBin  = dataArray{:, 6};
clearvars filename startRow formatSpec fileID dataArray ans;


QTranC  = (2.d0 .* pi .* CMass  .* DSWtoKg ./ AvN .* UKb .* T ./ Plnck.^2).^(3/2);
QTranO  = (2.d0 .* pi .* OMass  .* DSWtoKg ./ AvN .* UKb .* T ./ Plnck.^2).^(3/2);
QTranCO = (2.d0 .* pi .* COMass .* DSWtoKg ./ AvN .* UKb .* T ./ Plnck.^2).^(3/2);

for iLevel = 1:length(Level)
  KcDiss(iLevel) = gC * gO / g(iLevel) * QTranC * QTranO / QTranCO * exp( EEh(iLevel) * Ue / (UKb * T));
end


for iLevel = 2:length(Level)
  for jLevel = 1:iLevel
    Kc(iLevel,jLevel) = g(jLevel) / g(iLevel) * exp( - (EEh(jLevel) - EEh(iLevel)) * Ue / (UKb * T));
  end
end



filename = DissRatesFile;
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
iLevel = dataArray{:, 1};
jLevel = dataArray{:, 2};
KDiss  = dataArray{:, 3};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

% KRec(:) = KDiss(:) ./ KcDiss(iLevel(:))';
% 

figure
semilogy(iLevel,KDiss,'og')
hold on
semilogy(iLevel,KRec,'or')


fid  = fopen(DissRatesFile);
fwrt = fopen(NEWDissRatesFile,'wt');

tline = fgetl(fid);
tline = strcat(tline,', "K(j,i)"\n');
fprintf(fwrt, tline);

iLine = 1;
while ischar(tline)
    %disp(tline)
    tline = fgetl(fid);
    if ischar(tline)
      newline = strcat(tline,',',num2str(KRec(iLine),'%10.4e\n'),'\n');
      fprintf(fwrt, newline);
      iLine = iLine + 1;
    end
end

fclose(fid);
fclose(fwrt);



filename = InelRatesFile;
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
iLevel = dataArray{:, 1};
jLevel = dataArray{:, 2};
Kf  = dataArray{:, 3};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

for iLine = 1:length(iLevel)
  Kb(iLine) = Kf(iLine) ./ Kc(iLevel(iLine),jLevel(iLine));
end


fid  = fopen(InelRatesFile);
fwrt = fopen(NEWInelRatesFile,'wt');

tline = fgetl(fid);
tline = strcat(tline,', "K(j,i)"\n');
fprintf(fwrt, tline);

iLine = 1;
while ischar(tline)
    %disp(tline)
    tline = fgetl(fid);
    if ischar(tline)
      newline = strcat(tline,',',num2str(Kb(iLine),'%10.4e\n'),'\n');
      fprintf(fwrt, newline);
      iLine = iLine + 1;
    end
end

fclose(fid);
fclose(fwrt);
