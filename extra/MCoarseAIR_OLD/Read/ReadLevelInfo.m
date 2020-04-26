% -- MATLAB --
%%==============================================================================================================
% 
% Coarse-Grained method for Quasi-Classical Trajectories (CG-QCT) 
% 
% Copyright (C) 2018 Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
%
% Based on "VVTC" (Vectorized Variable stepsize Trajectory Code) by David Schwenke (NASA Ames Research Center). 
% 
% This program is free software; you can redistribute it and/or modify it under the terms of the 
% Version 2.1 GNU Lesser General Public License as published by the Free Software Foundation. 
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU Lesser General Public License for more details. 
% 
% You should have received a copy of the GNU Lesser General Public License along with this library; 
% if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
% 
%---------------------------------------------------------------------------------------------------------------
%%==============================================================================================================

function [NLevels, Levelvqn, Leveljqn, LevelEh, LevelEeV, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelg, LevelToBin, LevelQ, vToLevel, DeltaEintDiss, DeltaEintDepth, rIn, rOut, rMin, rMax, VMin, VMax, Tau, Egam, dVIn, ddVIn, dVOut, ddVOut, dVdJIn, dVdJOut, ErrorVec, DeltaEintDissJ] = ReadLevelInfo(iT, EeV, ErrorVec);
      
    
  %% (METHOD DEPENDENT)

  global T0_Vec System SystemPath 
  global NBins BinnedMolName NBinnedMol KeV KinMthd EhToeV MoleculesName AtomMass Pair_to_Atoms ParaViewFlg

  
  for iBinnedMol=1:NBinnedMol       
    ErrorVec = [];
    
 
    filename = strcat(SystemPath,'/',BinnedMolName(iBinnedMol,:),'/Bins_',num2str(NBins(iBinnedMol)),'/QNsEnBin.csv')
    opts = delimitedTextImportOptions("NumVariables", 6);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["Level", "vqn", "jqn", "EEh", "Degeneracy", "Bin"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    tbl = readtable(filename, opts);
    NLevels(iBinnedMol) = size(tbl.Level,1);
    Levelvqn(1:NLevels(iBinnedMol),iBinnedMol)   = tbl.vqn;
    Leveljqn(1:NLevels(iBinnedMol),iBinnedMol)   = tbl.jqn;
    LevelEeV(1:NLevels(iBinnedMol),iBinnedMol)   = tbl.EEh;
    Levelg(1:NLevels(iBinnedMol),iBinnedMol)     = tbl.Degeneracy;
    LevelToBin(1:NLevels(iBinnedMol),iBinnedMol) = tbl.Bin;

    %if sum(KinMthd(iBinnedMol,:) =='STS')==3
    
%         filename = strcat(SystemPath,'/',BinnedMolName(iBinnedMol,:),'/levels_cut.inp')
%         startRow = 16;
%         formatSpec = '%6f%5f%15f%[^\n\r]';
%         fileID = fopen(filename,'r');
%         dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', 0.0, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%         fclose(fileID);
%         vqnVec = dataArray{:, 1};
%         jqnVec = dataArray{:, 2};
%         EintVec = dataArray{:, 3};
%         clearvars filename startRow formatSpec fileID dataArray ans;
% 
%         Egam(1:NLevels(iBinnedMol),iBinnedMol) = 0.0;
% 
%         filename = strcat(SystemPath,'/',BinnedMolName(iBinnedMol,:),'/levels_cut.inp')
%         startRow = 16;
%         formatSpec = '%*41s%15f%15f%15f%15f%15f%15f%16f%[^\n\r]';
%         fileID = fopen(filename,'r');
%         dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', 0.0, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%         fclose(fileID);
%         rMin = dataArray{:, 1};
%         rMax = dataArray{:, 2};
%         VMin = dataArray{:, 3}.*EhToeV;
%         VMax = dataArray{:, 4}.*EhToeV;
%         Tau = dataArray{:, 5};
%         rIn = dataArray{:, 6};
%         rOut = dataArray{:, 7};
%         clearvars filename startRow formatSpec fileID dataArray ans;
      
      filename = strcat(SystemPath,'/',BinnedMolName(iBinnedMol,:),'/levels_cut.inp')
      startRow = 16;
      formatSpec = '%6f%5f%15f%15f%15f%15f%15f%15f%15f%15f%f%[^\n\r]';
      fileID = fopen(filename,'r');
      dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
      fclose(fileID);
      vqnVec  = dataArray{:, 1};
      jqnVec  = dataArray{:, 2};
      EintVec = dataArray{:, 3};
      Egam(1:NLevels(iBinnedMol),iBinnedMol) = max(dataArray{:, 4},1.d-300);
      rMin(1:NLevels(iBinnedMol),iBinnedMol) = dataArray{:, 5};
      rMax(1:NLevels(iBinnedMol),iBinnedMol) = dataArray{:, 6};
      VMin(1:NLevels(iBinnedMol),iBinnedMol) = dataArray{:, 7}.*EhToeV;
      VMax(1:NLevels(iBinnedMol),iBinnedMol) = dataArray{:, 8}.*EhToeV;
      Tau(1:NLevels(iBinnedMol),iBinnedMol)  = dataArray{:, 9};
      rIn(1:NLevels(iBinnedMol),iBinnedMol)  = dataArray{:, 10};
      rOut(1:NLevels(iBinnedMol),iBinnedMol) = dataArray{:, 11};
      clearvars filename startRow formatSpec fileID dataArray ans;
      
      
      LevelEeV(1:NLevels(iBinnedMol),iBinnedMol)   = LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) - VMax(1,iBinnedMol);
      LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol)  = LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) - min(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol));
      LevelEh(1:NLevels(iBinnedMol),iBinnedMol)    = LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) ./ EhToeV;
      
      VMax(1:NLevels(iBinnedMol),iBinnedMol)       = VMax(1:NLevels(iBinnedMol),iBinnedMol) - VMax(1,iBinnedMol);
      VMin(1:NLevels(iBinnedMol),iBinnedMol)       = VMin(1:NLevels(iBinnedMol),iBinnedMol) - VMax(1,iBinnedMol);

% %       if iBinnedMol == 1
% %         figure
%         VMin(1:NLevels(iBinnedMol),iBinnedMol) = 0.d0;
%         VMax(1:NLevels(iBinnedMol),iBinnedMol) = 0.d0;
%         [rMin0, VMin0, rMax0, VMax0, VInfTemp] = FindMinMaxDiatPot(0, iBinnedMol, rMin(1), rMax(1));
%         rMinOld = rMin0;
%         rMaxOld = rMax0;
%         for jqn = 1:max(Leveljqn(:,iBinnedMol))
%           jqn
%           [rMinTemp(jqn), VMinTemp(jqn), rMaxTemp(jqn), VMaxTemp(jqn), VInfTemp] = FindMinMaxDiatPot(jqn, iBinnedMol, rMinOld, rMaxOld);
%           rMinOld = rMinTemp(jqn);
%           rMaxOld = rMaxTemp(jqn);
%         end
        dVIn(:,iBinnedMol)    = rMin .* 0.d0;
        ddVIn(:,iBinnedMol)   = rMin .* 0.d0;
        dVOut(:,iBinnedMol)   = rMin .* 0.d0;
        ddVOut(:,iBinnedMol)  = rMin .* 0.d0;
        dVdJIn(:,iBinnedMol)  = rMin .* 0.d0;
        dVdJOut(:,iBinnedMol) = rMin .* 0.d0;
%         iii = 0;
% 
%         for iLevels = 1:NLevels(iBinnedMol)
%           %[dVIn(iLevels,iBinnedMol), ddVIn(iLevels,iBinnedMol), dVdJIn(iLevels,iBinnedMol)]    = FindDerivativeDiatPot(rIn(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), iBinnedMol, Leveljqn);
%           %[dVOut(iLevels,iBinnedMol), ddVOut(iLevels,iBinnedMol), dVdJOut(iLevels,iBinnedMol)] = FindDerivativeDiatPot(rOut(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), iBinnedMol, Leveljqn);
%           if Leveljqn(iLevels,iBinnedMol) ~= 0
%             if VMaxTemp(Leveljqn(iLevels,iBinnedMol)) > LevelEeV(iLevels,iBinnedMol) 
%               rMin(iLevels,iBinnedMol) = rMinTemp(Leveljqn(iLevels,iBinnedMol)); 
%               rMax(iLevels,iBinnedMol) = rMaxTemp(Leveljqn(iLevels,iBinnedMol)); 
%               VMin(iLevels,iBinnedMol) = VMinTemp(Leveljqn(iLevels,iBinnedMol));%.*EhToeV; 
%               VMax(iLevels,iBinnedMol) = VMaxTemp(Leveljqn(iLevels,iBinnedMol));%.*EhToeV;
%               [rIn(iLevels,iBinnedMol), rOut(iLevels,iBinnedMol), ErrorFlg] = TurningPoints( [1.20d0, rMin(iLevels,iBinnedMol)], [rMin(iLevels,iBinnedMol), rMax(iLevels,iBinnedMol)], LevelEeV(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), iBinnedMol);
%             else
%               iii=iii+1;
%               ErrorVec(iii) = iLevels;
%             end
%           else
%             rMin(iLevels,iBinnedMol) = rMin0; 
%             rMax(iLevels,iBinnedMol) = rMax0; 
%             VMin(iLevels,iBinnedMol) = VMin0;%.*EhToeV; 
%             VMax(iLevels,iBinnedMol) = 0.d0;%VMax0.*EhToeV;
%             [rIn(iLevels,iBinnedMol), rOut(iLevels,iBinnedMol), ErrorFlg] = TurningPoints( [1.20d0, rMin(iLevels,iBinnedMol)], [rMin(iLevels,iBinnedMol), rMax(iLevels,iBinnedMol)], LevelEeV(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), iBinnedMol);
%           end
%         end    
%         
%         
% %         filename = '/home/venturi/Downloads/O2_Temp_Sorted.inp';
% %         startRow = 16;
% %         formatSpec = '%6s%5s%15s%15s%15s%15s%15s%15s%15s%15s%s%[^\n\r]';
% %         fileID = fopen(filename,'r');
% %         dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% %         fclose(fileID);
% %         raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
% %         for col=1:length(dataArray)-1
% %           raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
% %         end
% %         numericData = NaN(size(dataArray{1},1),size(dataArray,2));
% %         for col=[1,2,3,4,5,6,7,8,9,10,11]
% %           % Converts text in the input cell array to numbers. Replaced non-numeric
% %           % text with NaN.
% %           rawData = dataArray{col};
% %           for row=1:size(rawData, 1)
% %             % Create a regular expression to detect and remove non-numeric prefixes and
% %             % suffixes.
% %             regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
% %             try
% %               result = regexp(rawData(row), regexstr, 'names');
% %               numbers = result.numbers;
% % 
% %               % Detected commas in non-thousand locations.
% %               invalidThousandsSeparator = false;
% %               if numbers.contains(',')
% %                 thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
% %                 if isempty(regexp(numbers, thousandsRegExp, 'once'))
% %                   numbers = NaN;
% %                   invalidThousandsSeparator = true;
% %                 end
% %               end
% %               % Convert numeric text to numbers.
% %               if ~invalidThousandsSeparator
% %                 numbers = textscan(char(strrep(numbers, ',', '')), '%f');
% %                 numericData(row, col) = numbers{1};
% %                 raw{row, col} = numbers{1};
% %               end
% %             catch
% %               raw{row, col} = rawData{row};
% %             end
% %           end
% %         end
% %         R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
% %         raw(R) = {NaN}; % Replace non-numeric cells
% %         vqnTemp  = cell2mat(raw(:, 1));
% %         jqnTemp  = cell2mat(raw(:, 2));
% %         for iLevels=1:length(Levelvqn)
% %           for iTemp=1:length(vqnTemp)
% %             if Levelvqn(iLevels) == vqnTemp(iTemp) && Leveljqn(iLevels) == jqnTemp(iTemp)
% %               LevelEh(iLevels)   = cell2mat(raw(iTemp, 3)) - cell2mat(raw(1, 8));
% %               Egam(iLevels)      = cell2mat(raw(iTemp, 4));
% %               rMin(iLevels)      = cell2mat(raw(iTemp, 5));
% %               rMax(iLevels)      = cell2mat(raw(iTemp, 6));
% %               VMin(iLevels)      = (cell2mat(raw(iTemp, 7)) - cell2mat(raw(1, 8)) .* EhToeV);
% %               VMax(iLevels)      = (cell2mat(raw(iTemp, 8)) - cell2mat(raw(1, 8)) .* EhToeV);
% %               Tau(iLevels)       = cell2mat(raw(iTemp, 9));
% %               rIn(iLevels)       = cell2mat(raw(iTemp, 10));
% %               rOut(iLevels)      = cell2mat(raw(iTemp, 11));
% %               LevelEeV(iLevels)  = LevelEh(iLevels) .* EhToeV;
% %               LevelEeV0(iLevels) = LevelEeV(iLevels) - LevelEeV(1);
% %             end
% %           end
% %         end
% %         clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;
% 
%         
%         %%% Compute Vibrational Energy-Distance From Centrifugal Barrier
        DeltaEintDiss(1:NLevels(iBinnedMol),iBinnedMol)  = LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) - VMax(1:NLevels(iBinnedMol),iBinnedMol); 
        DeltaEintDepth(1:NLevels(iBinnedMol),iBinnedMol) = LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) - VMin(1:NLevels(iBinnedMol),iBinnedMol); 
      
%     %end
    
    if sum(KinMthd(iBinnedMol,:) == 'CGM')==3 || sum(KinMthd(iBinnedMol,:) == 'VIB')==3
      LevelQ(1:NLevels(iBinnedMol),iBinnedMol) = Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* exp( - (LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol)) ./ (KeV*T0_Vec(iT)) );
    elseif sum(KinMthd(iBinnedMol,:) =='STS')==3
      %LevelQ(1:NLevels(iBinnedMol),iBinnedMol) = Levelg(1:NLevels(iBinnedMol),iBinnedMol);
      LevelQ(1:NLevels(iBinnedMol),iBinnedMol) = Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* exp( - (LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol)) ./ (KeV*T0_Vec(iT)) );
    end 
    LevelQ(1:NLevels(iBinnedMol),iBinnedMol) = LevelQ(1:NLevels(iBinnedMol),iBinnedMol) ./ sum(LevelQ(1:NLevels(iBinnedMol),iBinnedMol),1);  
    
    for iLevels = 1:NLevels(iBinnedMol)
      if (Leveljqn(iLevels,iBinnedMol) == 0) 
        vToLevel(Levelvqn(iLevels,iBinnedMol)+1,iBinnedMol) = iLevels; 
        vEeVVib(Levelvqn(iLevels,iBinnedMol)+1,iBinnedMol)  = LevelEeV(iLevels,iBinnedMol); 
        vEeVVib0(Levelvqn(iLevels,iBinnedMol)+1,iBinnedMol) = LevelEeV0(iLevels,iBinnedMol);         
      end
      LevelEeVVib0(iLevels,iBinnedMol) = LevelEeV0(vToLevel(Levelvqn(iLevels,iBinnedMol)+1,iBinnedMol),iBinnedMol);
      LevelEeVRot(iLevels,iBinnedMol)  = LevelEeV0(iLevels,iBinnedMol) - LevelEeVVib0(iLevels,iBinnedMol);      
    end
    
%     for iLevels = 1:NLevels(iBinnedMol)
%       iLevels
%       LevelvTry(iLevels,iBinnedMol) = 0;
%       while LevelvTry(iLevels,iBinnedMol) < max(max(Levelvqn)) &&  LevelEeV(iLevels,iBinnedMol) > vEeVVib(LevelvTry(iLevels,iBinnedMol)+1,iBinnedMol)
%         LevelvTry(iLevels,iBinnedMol) = LevelvTry(iLevels,iBinnedMol) + 1;
%       end
%     end

    %%% Compute Rotational Energy-Distance From Centrifugal Barrier
    DeltaEintDissVTemp = 0.0 .* LevelEeV;
%     for ivqn = 0:max(Levelvqn(:,iBinnedMol))
%       ivqn
%       
%       clear TEMPP TempEj TempDelta
%       TempEv = [];
%       for iLevel=1:NLevels(iBinnedMol)
%         if Levelvqn(iLevel,iBinnedMol) == ivqn
%           TempEj(Leveljqn(iLevel,iBinnedMol)+1,1)    = LevelEeV(iLevel,iBinnedMol);
%           TempDelta(Leveljqn(iLevel,iBinnedMol)+1,1) = DeltaEintDiss(iLevel,iBinnedMol);
%         end
%       end
% 
%       fitDelta = fit(TempEj,TempDelta,'smoothingspline');
%       fun      = fitDelta;
%       x0       = [TempEj(end), TempEj(end) + (max(Levelvqn(:,iBinnedMol))+2-ivqn).*(TempEj(end) - TempEj(end-1))];
%       x        = fzero(fun,x0);
%       fun(x)
% 
% %       figure
% %       plot(TempEj, TempDelta,'ko')
% %       hold on
% %       TEMPP = linspace(TempEj(1), TempEj(end) + (max(Levelvqn)+2-ivqn).*(TempEj(end) - TempEj(end-1)), 1000);
% %       plot(TEMPP, fitDelta(TEMPP),'-')
% %       plot(x,fun(x),'ro')
% 
%       DeltaEintDissVTemp(ivqn+1) = x;
%     end
    
    for iLevels = 1:NLevels(iBinnedMol)
      DeltaEintDissJ(iLevels,iBinnedMol) = -DeltaEintDissVTemp(Levelvqn(iLevels,iBinnedMol)+1) + LevelEeV(iLevels,iBinnedMol);
    end
    clear DeltaEintDissVTemp
    
    
    if ParaViewFlg == 1
      fileName = strcat('./',MoleculesName(iBinnedMol,:),'_LevelsProperties.csv');
      fileID   = fopen(fileName,'w');
      fprintf(fileID,'Variables = "rIn", "rOut", "J", "v", "E", "DeltaEDiss", "DeltaEintDissJ", "ERot", "EVib", "g", "Q", "dVIn", "dVOut", "rMin", "VMin", "rMax", "VMax", "Tau", "Egam"\n');
      for i = 1:NLevels(iBinnedMol)
        fprintf(fileID,'%15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e\n', ...
        rIn(i,iBinnedMol), rOut(i,iBinnedMol), Leveljqn(i,iBinnedMol), Levelvqn(i,iBinnedMol), EeV(i,iBinnedMol)+LevelEeV(1,iBinnedMol), DeltaEintDiss(i,iBinnedMol), DeltaEintDissJ(i,iBinnedMol), LevelEeVRot(i,iBinnedMol), LevelEeVVib0(i,iBinnedMol), Levelg(i,iBinnedMol), LevelQ(i,iBinnedMol), ...
        dVIn(i,iBinnedMol), dVOut(i,iBinnedMol), rMin(i,iBinnedMol), VMin(i,iBinnedMol), rMax(i,iBinnedMol), VMax(i,iBinnedMol), Tau(i,iBinnedMol), Egam(i,iBinnedMol));
      end
      fclose(fileID);
    end
    
    
    
  end
    
end
