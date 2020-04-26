close all
clc

iBinnedMol = 1;


filename = '//Users/sventuri/Desktop/SSSS/Inelastic/output1/InelRates1.clu';
delimiter = ' ';
startRow = 3;
formatSpec = '%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
LevelNb     = dataArray{:, 1};
GroupNbTemp = dataArray{:, 2};
Rank        = dataArray{:, 3};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
GroupNbOrig(LevelNb(:)) = GroupNbTemp(:);
GroupNbBegin = GroupNbOrig(1);

% for iLevels = 1:NLevels
%   if GroupNbOrig(iLevels) > GroupNbBegin
%     GroupNbOrig(iLevels) = GroupNbOrig(iLevels) - 1;
%   elseif GroupNbOrig(iLevels) == GroupNbBegin
%     GroupNbOrig(iLevels) = 0;
%   end
% end

MaxGroupNb = max(GroupNbOrig);

GroupNb = GroupNbOrig;



% filename = '/Users/sventuri/Dropbox/TempRes/output/InelRates.tree';
% delimiter = {' ',':'};
% startRow = 3;
% formatSpec = '%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% A1 = dataArray{:, 1};
% A2 = dataArray{:, 2};
% A3 = dataArray{:, 3};
% A4 = dataArray{:, 4};
% clearvars filename delimiter startRow formatSpec fileID dataArray ans;
% 
% filename = '/Users/sventuri/Dropbox/TempRes/output/InelRates.tree';
% delimiter = ' ';
% startRow = 3;
% formatSpec = '%*q%*q%*q%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% LevelNb = dataArray{:, 1};
% clearvars filename delimiter startRow formatSpec fileID dataArray ans;
% GroupNbTemp(LevelNb(:)) = (A1(:)-1).*max(A2) + A2(:);
% %GroupNb(LevelNb(:)) = (A1(:)-1).*max(A2) + (A2(:)-1).*max(A3) + A3(:);
% 
% MaxGroupNb = max(A1);
% for iLevels = 1:NLevels(iBinnedMol)
%   if LevelNb(iLevels) == 1
%     GroupNbBegin = A1(iLevels)
%   end
% end
% 
% for iLevels = 1:NLevels(iBinnedMol)
%   if A1(iLevels) > GroupNbBegin
%     A1(iLevels) = A1(iLevels) - 1;
%   elseif A1(iLevels) == GroupNbBegin
%     A1(iLevels) = 0;
%   end
% end
% 
% for iLevels = 1:NLevels(iBinnedMol)
%   if A1(iLevels) == 0
%     GroupNb(LevelNb(iLevels)) = MaxGroupNb - 1 + A2(iLevels);
%   else
%     GroupNb(LevelNb(iLevels)) = A1(iLevels);
%   end
% end


FileName = strcat('/Users/sventuri/Desktop/SSSS/InelLevelsNew.csv');
fileID = fopen(FileName,'w');
fprintf(fileID,'id,v,Longitude,Latitude,rIn,EeVVib,EeVRot,Group\n');
for i = 1:NLevels(1)
  fprintf(fileID,'%i,%i,%e,%e,%e,%e,%e,%i\n', i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), LevelEeVVib0(i), LevelEeVRot(i), GroupNb(i));
end
fclose(fileID);
