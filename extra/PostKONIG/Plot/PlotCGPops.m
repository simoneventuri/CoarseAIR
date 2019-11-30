clc
%close all


iBinnedMol = 1;
iFigure    = 2001;


filename = '/Users/sventuri/WORKSPACE/neqplasma_QCT/konig/runs/CO2/output/pop_CO.dat';
delimiter = ' ';
startRow = 3;
formatSpec = '%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
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
nTemp = cell2mat(raw(:, 2));
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me I J K;

i = 1;
for iSteps = 1:length(t)
  for iBins = 1:NbBins
    nBinsKONIGOverg(iBins,iSteps) = nTemp(i);
    i = i + 1;
  end
end


TempStpInstants = StpInstants;
for iSteps = TempStpInstants
  
  nKONIG(1:NLevels(iBinnedMol),iSteps)      = Q(1:NLevels(iBinnedMol)) .* nBinsKONIGOverg(LevToBinFinal(1:NLevels(iBinnedMol)),iSteps);
  nKONIGOverG(1:NLevels(iBinnedMol),iSteps) = nKONIG(1:NLevels(iBinnedMol),iSteps) ./ Levelg(1:NLevels(iBinnedMol),iBinnedMol);
  
end
TempStpInstants =[2:10:2000];
for iSteps = TempStpInstants
  
  nKONIG(1:NLevels(iBinnedMol),iSteps)      = Q(1:NLevels(iBinnedMol)) .* nBinsKONIGOverg(LevToBinFinal(1:NLevels(iBinnedMol)),iSteps);
  nKONIGOverG(1:NLevels(iBinnedMol),iSteps) = nKONIG(1:NLevels(iBinnedMol),iSteps) ./ Levelg(1:NLevels(iBinnedMol),iBinnedMol);
  
end

for iSteps = StpInstants
  
  figure(iFigure)
  hold on
  scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),nKONIGOverG(1:NLevels(iBinnedMol),iSteps),20,LevToBin(1:NLevels(iBinnedMol),iBinnedMol,1),'Filled');
  colormap(lines(max(max(max(LevToBin)))))
  cb=colorbar;
  %cb.Ticks = [1, 1.5]; %Create 8 ticks from zero to 1
  %cb.TickLabels = {'1','2'}
  ylab = ylabel(cb, '$Group$');
  ylab.Interpreter = 'latex';
  set(cb,'FontSize',LegendFontSz,'FontName',LegendFontNm,'TickLabelInterpreter','latex');
  hold on
  xt = get(gca, 'XTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  yt = get(gca, 'YTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  str_x = ['Energy [eV]'];
  xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  xlab.Interpreter = 'latex';
  xlim([max(min(LevelEeV(:,iBinnedMol)),MinEvPlot(iBinnedMol)), min(max(LevelEeV(:,iBinnedMol)),MaxEvPlot(iBinnedMol))]); 
  str_y = ['$N_{i} / g_{i} [m^{-3}]$'];
  ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  ylab.Interpreter = 'latex';
  ylim([1.d5, 1.d20]);
  set(gca, 'YScale', 'log')
  iFigure = iFigure + 1;
  
end