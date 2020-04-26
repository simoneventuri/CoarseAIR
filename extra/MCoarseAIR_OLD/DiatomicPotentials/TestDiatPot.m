%close all
clear all
clc

global System T0_Vec
global NBins NTint NAtoms AtomsName MoleculesName NMolecules DegeneracyFactor BinnedMolName ColPartToComp NBinnedMol BinnedMolToComp NComp ...
       CompNames CompColor AtomColor AtomSize AllMoleculesName PairColor AtomMass ComponentMass ColorVec ComponentDeg MoleculeMu RxLxIdx MoleculedDissEn
global BCVec Pair_to_Atoms

System = 'H3'
T0_Vec = []
Pair_to_Atoms = [1,2; 1,3; 2,3];

UpdateAllInput()


R = linspace(0.8,10.0,1000);

figure(1)
for iMol=1:1
  jqn     = 0;
  [V, dV] = DiatPot(100.0, jqn, iMol)
  [V, dV] = DiatPot(R, jqn, iMol);
  plot(R,V)
  hold on
  %plot(R,dV)
end
legend(MoleculesName)


xlim([0.6,10])
ylim([-10,60])


% 
% 
% filename = '/Users/sventuri/WORKSPACE/SPES/spes/Data_DiatPot/CHN/MRCI_NH.dat';
% delimiter = ',';
% startRow = 2;
% % For more information, see the TEXTSCAN documentation.
% formatSpec = '%s%s%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
% for col=1:length(dataArray)-1
%   raw(1:length(dataArray{col}),col) = dataArray{col};
% end
% numericData = NaN(size(dataArray{1},1),size(dataArray,2));
% for col=[1,2]
%   % Converts text in the input cell array to numbers. Replaced non-numeric
%   % text with NaN.
%   rawData = dataArray{col};
%   for row=1:size(rawData, 1);
%     % Create a regular expression to detect and remove non-numeric prefixes and
%     % suffixes.
%     regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
%     try
%       result = regexp(rawData{row}, regexstr, 'names');
%       numbers = result.numbers;
%       
%       % Detected commas in non-thousand locations.
%       invalidThousandsSeparator = false;
%       if any(numbers==',');
%         thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
%         if isempty(regexp(numbers, thousandsRegExp, 'once'));
%           numbers = NaN;
%           invalidThousandsSeparator = true;
%         end
%       end
%       % Convert numeric text to numbers.
%       if ~invalidThousandsSeparator;
%         numbers = textscan(strrep(numbers, ',', ''), '%f');
%         numericData(row, col) = numbers{1};
%         raw{row, col} = numbers{1};
%       end
%     catch me
%     end
%   end
% end
% RPoints = cell2mat(raw(:, 1));
% VPoints = cell2mat(raw(:, 2));
% VPoints = ( VPoints - VPoints(end) ) .* 27.2113839712790;
% clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me;
% 
% hold on
% plot(RPoints,VPoints,'ro')
