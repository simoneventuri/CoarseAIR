close all
clc
clear all

filename = './N4_AbInitio_UMN.txt';
startRow = 21;
formatSpec = '%12s%12s%12s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
  raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3]
  % Converts text in the input cell array to numbers. Replaced non-numeric
  % text with NaN.
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
raw(I,:) = [];
xTemp = cell2mat(raw(:, 1));
yTemp = cell2mat(raw(:, 2));
zTemp = cell2mat(raw(:, 3));
clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me I J K;


filename = './N4_AbInitio_UMN.txt';
startRow = 21;
formatSpec = '%6s%13s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
  raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2]
  % Converts text in the input cell array to numbers. Replaced non-numeric
  % text with NaN.
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
raw(I,:) = [];
Temp = cell2mat(raw(:, 1));
EEhTemp  = cell2mat(raw(:, 2));
clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me I J K;

NPoints = length(xTemp)/4;

clear x y z EEh
iPoint = 1;
iLine1 = 1;
iLine2 = 1;
R      = [];
while iPoint <= NPoints
  EEh(iPoint) = EEhTemp(iLine1) .* 0.159360144d-2 ;
  iLine1 = iLine1+1;
  for iA=1:4
    x(iA,iPoint) = xTemp(iLine2) ./ 0.52917721092;
    y(iA,iPoint) = yTemp(iLine2) ./ 0.52917721092;
    z(iA,iPoint) = zTemp(iLine2) ./ 0.52917721092;
    iLine2 = iLine2+1;
  end
  R(1,iPoint) = sqrt( (x(1,iPoint) - x(2,iPoint)).^2 + (y(1,iPoint) - y(2,iPoint)).^2 + (z(1,iPoint) - z(2,iPoint)).^2 );
  R(2,iPoint) = sqrt( (x(1,iPoint) - x(3,iPoint)).^2 + (y(1,iPoint) - y(3,iPoint)).^2 + (z(1,iPoint) - z(3,iPoint)).^2 );
  R(3,iPoint) = sqrt( (x(1,iPoint) - x(4,iPoint)).^2 + (y(1,iPoint) - y(4,iPoint)).^2 + (z(1,iPoint) - z(4,iPoint)).^2 );
  R(4,iPoint) = sqrt( (x(2,iPoint) - x(3,iPoint)).^2 + (y(2,iPoint) - y(3,iPoint)).^2 + (z(2,iPoint) - z(3,iPoint)).^2 );
  R(5,iPoint) = sqrt( (x(2,iPoint) - x(4,iPoint)).^2 + (y(2,iPoint) - y(4,iPoint)).^2 + (z(2,iPoint) - z(4,iPoint)).^2 );
  R(6,iPoint) = sqrt( (x(3,iPoint) - x(4,iPoint)).^2 + (y(3,iPoint) - y(4,iPoint)).^2 + (z(3,iPoint) - z(4,iPoint)).^2 );
  iPoint = iPoint + 1;
end

figure(1)
scatter3(x(1,:),y(1,:),z(1,:))
hold on
scatter3(x(2,:),y(2,:),z(2,:))
scatter3(x(3,:),y(3,:),z(3,:))
scatter3(x(4,:),y(4,:),z(4,:))

figure(2)
for iPoint = 1:NPoints
  plot3(x(:,iPoint),y(:,iPoint),z(:,iPoint),'o-')
  pause(0.1)
end

FileName1 = strcat('./N2_N2_Scans.csv');
fileID1   = fopen(FileName1,'w');
fprintf(fileID1,'R1,R2,R3,R4,R5,R6\n');
for iPoints=1:15365
  fprintf(fileID1,'%e,%e,%e,%e,%e,%e\n', R(1,iPoints), R(2,iPoints), R(3,iPoints), R(4,iPoints), R(5,iPoints), R(6,iPoints) );
end
fclose(fileID1);

FileName1 = strcat('./N3_N_Scans.csv');
fileID1   = fopen(FileName1,'w');
fprintf(fileID1,'R1,R2,R3,R4,R5,R6\n');
for iPoints=15366:16382
  fprintf(fileID1,'%e,%e,%e,%e,%e,%e\n', R(1,iPoints), R(2,iPoints), R(3,iPoints), R(4,iPoints), R(5,iPoints), R(6,iPoints) );
end
fclose(fileID1);

FileName1 = strcat('./N3_Scans.csv');
fileID1   = fopen(FileName1,'w');
fprintf(fileID1,'R1,R2,R3,R4,R5,R6\n');
for iPoints=15366:16382
  fprintf(fileID1,'%e,%e,%e,%e,%e,%e\n', R(1,iPoints), R(2,iPoints), R(3,iPoints), R(4,iPoints), R(5,iPoints), R(6,iPoints) );
end
fclose(fileID1);


for iPoints=1:NPoints
  if (max(R(:,iPoints)) > 12.0)
    iPoints
    R(:,iPoints)
    pause
    clc
  end
end