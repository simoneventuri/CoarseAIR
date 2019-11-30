function [LevToBin, BinEnergy, BinEnergyMax, qnToBin] = ReadLevToBin(NLevels, Levelvqn, Leveljqn, LevelEeV, DeltaEintDiss, StpInstants, GroupNb)

  global SystemPath BinnedMolName NBinnedMol Ue T0_Vec UKb
  
  LevToBin         = zeros(size(Levelvqn,1),2,length(StpInstants));
  BinEnergy        = 0.d0;
  BinEnergyMax     = 0.d0;
  
  for iMol = 1:1%NBinnedMol
    
%  for iSteps = 1:length(StpInstants)-1
  ii = 0;
  for iSteps = StpInstants

    ii = ii + 1;
%     %filename = strcat('/Users/sventuri/Dropbox/TempRes/Results/Inel/WOExchange/EqDeviation/',num2str(iSteps),'.csv')
%     filename = '/Users/sventuri/Dropbox/Gephi_Results/CO_Inelastic_NoExchange_10000K/Modularity0.01/0.01.csv'
%     delimiter = ',';
%     startRow = 2;
%     formatSpec = '%f%s%s%f%f%f%f%f%f%f%[^\n\r]';
%     %formatSpec = '%f%s%s%f%f%f%f%f%f%[^\n\r]';
%     fileID = fopen(filename,'r');
%     dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%     fclose(fileID);
%     Id = dataArray{:, 1};
%     modularity_class = dataArray{:, 10};
%     %modularity_class = dataArray{:, 9};
%     clearvars filename delimiter startRow formatSpec fileID dataArray ans;
%     length(Id)
%     if length(Id) == NLevels(iMol)
%       LevToBin(Id(:),iMol,ii) = modularity_class(:) + 1;
%     else
%       LevToBin(1:NLevels(iMol),iMol,ii) = 0.d0 .* Levelvqn(1:NLevels(iMol),iMol) + 2;
%       LevToBin(Id(:),iMol,ii) = modularity_class(:) + 1;
%     end
    
      LevToBin(1:NLevels(iMol),iMol,ii) = GroupNb(1:NLevels(iMol));

  end
  
    
%     NLevToBin = LevToBin.*0.d0;
%     for i = 1:size(LevToBin,1)
%       NLevToBin(LevToBin(i)) = NLevToBin(LevToBin(i)) + 1;
%     end
%     for i = 1:size(LevToBin,1)
%       if NLevToBin(LevToBin(i)) <= 1
%         LevToBin(i) = 0;
%       end
%     end
%   
    
%     for i = 1:size(LevToBin,1)-2
%       if LevToBin(i,1) > 50
%         LevToBin(i,1) = 50;
%       end
%       qnToBin(Levelvqn(i,iMol)+1,Leveljqn(i,iMol)+1) = LevToBin(i,iMol);
%     end
%     
  end
  
end