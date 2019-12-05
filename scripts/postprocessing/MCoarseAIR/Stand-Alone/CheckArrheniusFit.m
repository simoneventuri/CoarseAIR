close all
clear all
clc


NBins           = 20
TVec            = [1000, 2500, 5000, 6000, 7500, 10000, 15000]
RatesFldr       = '/Users/sventuri/WORKSPACE/CoarseAIR/run_O3_PES9_CG20/Test/O3/O2/Rates'
fileKinetics    = '/Users/sventuri/WORKSPACE/CoarseAIR/run_O3_PES9_CG20/Test/RunHegel/database/kinetics/O3'
SeperateExchFlg = 0


delimiter = {',',';',':'};
startRow = 2;
formatSpec = '%*s%f%f%f%*s%[^\n\r]';
fileID = fopen(fileKinetics,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
C1 = dataArray{:, 1};
C2 = dataArray{:, 2};
C3 = dataArray{:, 3};
clearvars delimiter startRow formatSpec fileID dataArray ans;

delimiter = {'+','_'};
startRow = 2;
formatSpec = '%*s%f%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(fileKinetics,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
IniState = dataArray{:, 1};
clearvars delimiter startRow formatSpec fileID dataArray ans;

delimiter = {'+','_'};
startRow = 2;
formatSpec = '%*s%*s%*s%s%[^\n\r]';
fileID = fopen(fileKinetics,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
  raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
rawData = dataArray{1};
for row=1:size(rawData, 1);
  % Create a regular expression to detect and remove non-numeric prefixes and
  % suffixes.
  regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
  try
    result = regexp(rawData{row}, regexstr, 'names');
    numbers = result.numbers;
    
    % Detected commas in non-thousand locations.
    invalidThousandsSeparator = false;
    if any(numbers==',');
      thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
      if isempty(regexp(numbers, thousandsRegExp, 'once'));
        numbers = NaN;
        invalidThousandsSeparator = true;
      end
    end
    % Convert numeric text to numbers.
    if ~invalidThousandsSeparator;
      numbers = textscan(strrep(numbers, ',', ''), '%f');
      numericData(row, 1) = numbers{1};
      raw{row, 1} = numbers{1};
    end
  catch me
  end
end
R = cellfun(@(x) (~isnumeric(x) && ~islogical(x)) || isnan(x),raw); % Find non-numeric cells
raw(R) = {0.0}; % Replace non-numeric cells
FinState = cell2mat(raw(:, 1));
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;
FinState = FinState + 1;

for l=1:length(C1)
  i        = IniState(l);
  j        = FinState(l);
  ArrFit   = [C1(l), C2(l), C3(l)];
  KKinetics(i,j,:) = exp(log(ArrFit(1)) + ArrFit(2).*log(TVec) - ArrFit(3)./TVec);
end

K     = zeros(NBins,NBins*3+1,length(TVec));
KTemp = zeros(NBins,NBins*3+1,length(TVec));
for iBin = 1:NBins
  for iT = 1:length(TVec)
    TT = TVec(iT);
    filename = strcat(RatesFldr,'/T_',num2str(TT),'_',num2str(TT),'/Bin',num2str(iBin),'.dat');
    delimiter = ' ';
    startRow = 6;
    formatSpec = '%*s%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    FinState = dataArray{:, 1};
    TempRate = dataArray{:, 2};
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    for j=1:length(FinState)
      iFinState = FinState(j);
      KTemp(iBin,iFinState,iT) = TempRate(j);
    end   
  end
  K(iBin,1,:) = KTemp(iBin,1,:);
  kk = 2;
  ii = 1;
  if SeperateExchFlg == 1
    MaxInt = 2;
  else
    MaxInt = 4;
  end
  while ii < MaxInt
    jj=1;
    while jj <= NBins
      K(iBin,jj+1,:) = K(iBin,jj+1,:) + KTemp(iBin,kk,:);
      kk = kk+1;
      jj = jj+1;
    end
    ii = ii+1;
  end
end 

figure; 
semilogy(squeeze(K(20,:,:)),'k')
hold on
semilogy(squeeze(KKinetics(20,:,:)),'r')


figure; 
semilogy(squeeze(K(:,1,:)),'k')
hold on
semilogy(squeeze(KKinetics(:,1,:)),'r')
