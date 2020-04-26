close all
clc


NPESs = 100;

iFigure = 1;

for iA = [120, 140, 160, 180]
  
  figure(iFigure)

  filename = strcat('/Users/sventuri/WORKSPACE/CG-QCT/run_N3_NASA/Test/PlotPES/FixedR1_',num2str(iA),'.csv');
  delimiter = ',';
  startRow = 2;
  formatSpec = '%f%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  RVec = dataArray{:, 1};
  EVec = dataArray{:, 2};
  clearvars filename delimiter startRow formatSpec fileID dataArray ans;

  plot(RVec,EVec)
  hold on

  for iPES=1:NPESs
    
    filename = strcat('/Users/sventuri/WORKSPACE/CG-QCT/run_N3_CG_BNN_100/Test/PlotPES/FixedR1_',num2str(iA),'.csv.', num2str(iPES));
    delimiter = ',';
    startRow = 2;
    formatSpec = '%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    RVec = dataArray{:, 1};
    EVec = dataArray{:, 2};
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    plot(RVec,EVec)
    
  end

  filename = strcat('/Users/sventuri/WORKSPACE/CG-QCT/run_N3_CG_BNN_100/Test/PlotPES/FixedR1_Stats_',num2str(iA),'.csv');
  startRow = 2;
  formatSpec = '%16s%16s%16s%16s%s%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
  for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
  end
  numericData = NaN(size(dataArray{1},1),size(dataArray,2));
  for col=[1,2,3,4,5]
    rawData = dataArray{col};
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
          numericData(row, col) = numbers{1};
          raw{row, col} = numbers{1};
        end
      catch me
      end
    end
  end
  I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),2); % Find rows with non-numeric cells
  raw(I,:)  = [];
  RVec      = cell2mat(raw(:, 1));
  EMeanVec  = cell2mat(raw(:, 2));
  ESDVec    = cell2mat(raw(:, 3));
  EMinusVec = cell2mat(raw(:, 4));
  EPlusVec  = cell2mat(raw(:, 5));
  clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousan
  plot(RVec, EMinusVec)
  plot(RVec, EMeanVec)
  plot(RVec, EPlusVec)
  iFigure = iFigure+1;
  
end


figure(iFigure)

filename = '/Users/sventuri/WORKSPACE/CG-QCT/run_N3_NASA/Test/PlotPES/FixedR1R2_1.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
AngleVec = dataArray{:, 1};
EVec     = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

plot(AngleVec,EVec)
hold on

for iPES=1:NPESs
  
  filename = strcat('/Users/sventuri/WORKSPACE/CG-QCT/run_N3_CG_BNN_100/Test/PlotPES/FixedR1R2_1.csv.', num2str(iPES));
  delimiter = ',';
  startRow = 2;
  formatSpec = '%f%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  AngleVec = dataArray{:, 1};
  EVec     = dataArray{:, 2};
  clearvars filename delimiter startRow formatSpec fileID dataArray ans;
  plot(AngleVec,EVec)
  
end

filename = '/Users/sventuri/WORKSPACE/CG-QCT/run_N3_CG_BNN_100/Test/PlotPES/FixedR1R2_Stats_1.csv';
startRow = 2;
formatSpec = '%16s%16s%16s%16s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
  raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,4,5]
  rawData = dataArray{col};
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
        numericData(row, col) = numbers{1};
        raw{row, col} = numbers{1};
      end
    catch me
    end
  end
end
I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),2); % Find rows with non-numeric cells
raw(I,:)  = [];
AngleVec  = cell2mat(raw(:, 1));
EMeanVec  = cell2mat(raw(:, 2));
ESDVec    = cell2mat(raw(:, 3));
EMinusVec = cell2mat(raw(:, 4));
EPlusVec  = cell2mat(raw(:, 5));
clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousan
plot(AngleVec, EMinusVec)
plot(AngleVec, EMeanVec)
plot(AngleVec, EPlusVec)
figure = figure+1;
