% -- MATLAB --
%%==============================================================================================================
% 
% Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
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

function [NLevels, Levelvqn, Leveljqn, LevelEh, LevelEeV, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelg, LevelToBin, LevelQ, vToLevel, DeltaEintDiss, DeltaEintDepth, rIn, rOut, rMin, rMax, VMin, VMax, Tau, Egam, dVIn, ddVIn, dVOut, ddVOut, dVdJIn, dVdJOut] = ReadLevelInfo(iT, EeV);
      
    
  %% (METHOD DEPENDENT)

  global T0_Vec System SystemPath 
  global NBins BinnedMolName NBinnedMol KeV KinMthd EhToeV MoleculesName AtomMass Pair_to_Atoms ParaViewFlg

  
  for iBinnedMol=1:NBinnedMol

    if sum(KinMthd(iBinnedMol,:) == 'CGM')==3 || sum(KinMthd(iBinnedMol,:) == 'VIB')==3
      filename = strcat(SystemPath,'/',BinnedMolName(iBinnedMol,:),'/',BinnedMolName(iBinnedMol,:),'_',num2str(NBins(iBinnedMol)),'/qnsEnBin.dat')
    elseif sum(KinMthd(iBinnedMol,:) =='STS')==3
      filename = strcat(SystemPath,'/',BinnedMolName(iBinnedMol,:),'/',BinnedMolName(iBinnedMol,:),'_',num2str(NBins(iBinnedMol)),'/qnsEnBin.dat')
    end
    startRow = 2;
%     formatSpec = '%8f%7f%7f%15f%15f%f%[^\n\r]';
%     fileID = fopen(filename,'r');
%     dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%     formatSpec = '%8f%7f%7f%15f%15f%f%[^\n\r]';
%     fileID = fopen(filename,'r');
%     dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');        
    formatSpec = '%8f%7f%7f%15f%15f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    NLevels(iBinnedMol)=size(dataArray{:,1},1)
    Levelvqn(1:NLevels(iBinnedMol),iBinnedMol)   = dataArray{:, 2};
    Leveljqn(1:NLevels(iBinnedMol),iBinnedMol)   = dataArray{:, 3};
    LevelEeV(1:NLevels(iBinnedMol),iBinnedMol)   = dataArray{:, 4};
    Levelg(1:NLevels(iBinnedMol),iBinnedMol)     = dataArray{:, 5};
    LevelToBin(1:NLevels(iBinnedMol),iBinnedMol) = dataArray{:, 6};
    clearvars filename startRow formatSpec fileID dataArray ans;
    
    %if sum(KinMthd(iBinnedMol,:) =='STS')==3
      
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

%       if iBinnedMol == 1
%         figure
%         VMin(1:NLevels(iBinnedMol),iBinnedMol) = 0.d0;
%         VMax(1:NLevels(iBinnedMol),iBinnedMol) = 0.d0;
%         [rMin0, VMin0, rMax0, VMax0, VInfTemp] = FindMinMaxDiatPot(0, iBinnedMol);
%         for jqn = 1:max(Leveljqn(:,iBinnedMol))
%            jqn
%           [rMinTemp(jqn), VMinTemp(jqn), rMaxTemp(jqn), VMaxTemp(jqn), VInfTemp] = FindMinMaxDiatPot(jqn, iBinnedMol);
%         end
         dVIn    = rMin .* 0.d0;
         ddVIn   = rMin .* 0.d0;
         dVOut   = rMin .* 0.d0;
         ddVOut  = rMin .* 0.d0;
         dVdJIn  = rMin .* 0.d0;
         dVdJOut = rMin .* 0.d0;
%         for iLevels = 1:NLevels(iBinnedMol)
%           [dVIn(iLevels,iBinnedMol), ddVIn(iLevels,iBinnedMol), dVdJIn(iLevels,iBinnedMol)]    = FindDerivativeDiatPot(rIn(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), iBinnedMol, Leveljqn);
%           [dVOut(iLevels,iBinnedMol), ddVOut(iLevels,iBinnedMol), dVdJOut(iLevels,iBinnedMol)] = FindDerivativeDiatPot(rOut(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), iBinnedMol, Leveljqn);
%           if Leveljqn(iLevels,iBinnedMol) ~= 0
%             rMin(iLevels,iBinnedMol) = rMinTemp(Leveljqn(iLevels,iBinnedMol)); 
%             rMax(iLevels,iBinnedMol) = rMaxTemp(Leveljqn(iLevels,iBinnedMol)); 
%             VMin(iLevels,iBinnedMol) = VMinTemp(Leveljqn(iLevels,iBinnedMol)).*EhToeV; 
%             VMax(iLevels,iBinnedMol) = VMaxTemp(Leveljqn(iLevels,iBinnedMol)).*EhToeV;
%           else
%             rMin(iLevels,iBinnedMol) = rMin0; 
%             rMax(iLevels,iBinnedMol) = rMax0; 
%             VMin(iLevels,iBinnedMol) = VMin0.*EhToeV; 
%             VMax(iLevels,iBinnedMol) = 0.d0;%VMax0.*EhToeV;
%           end
%           [rIn(iLevels,iBinnedMol), rOut(iLevels,iBinnedMol), ErrorFlg] = TurningPoints( [1.50d0, rMin(iLevels,iBinnedMol)], [rMin(iLevels,iBinnedMol), rMax(iLevels,iBinnedMol)], LevelEh(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), iBinnedMol);
%         end
        DeltaEintDiss(1:NLevels(iBinnedMol),iBinnedMol)  = LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) - VMax(1:NLevels(iBinnedMol),iBinnedMol); 
        DeltaEintDepth(1:NLevels(iBinnedMol),iBinnedMol) = LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) - VMin(1:NLevels(iBinnedMol),iBinnedMol); 
%       end
      
    %end
    
    if sum(KinMthd(iBinnedMol,:) == 'CGM')==3 || sum(KinMthd(iBinnedMol,:) == 'VIB')==3
      LevelQ(1:NLevels(iBinnedMol),iBinnedMol) = Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* exp( - (LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol)) ./ (KeV*T0_Vec(iT)) );
    elseif sum(KinMthd(iBinnedMol,:) =='STS')==3
      %LevelQ(1:NLevels(iBinnedMol),iBinnedMol) = Levelg(1:NLevels(iBinnedMol),iBinnedMol);
      LevelQ(1:NLevels(iBinnedMol),iBinnedMol) = Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* exp( - (LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol)) ./ (KeV*T0_Vec(iT)) );
    end 
    LevelQ = LevelQ ./ sum(LevelQ,1);  
    
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
    
    
    if ParaViewFlg == 1
      fileName = strcat('./',MoleculesName(iBinnedMol,:),'_LevelsProperties.csv');
      fileID   = fopen(fileName,'w');
      fprintf(fileID,'Variables = "rIn", "rOut", "J", "v", "E", "DeltaEDiss", "ERot", "EVib", "g", "Q", "dVIn", "dVOut", "rMin", "VMin", "rMax", "VMax", "Tau", "Egam"\n');
      for i = 1:NLevels(iBinnedMol)
        fprintf(fileID,'%15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e\n', ...
        rIn(i,iBinnedMol), rOut(i,iBinnedMol), Leveljqn(i,iBinnedMol), Levelvqn(i,iBinnedMol), EeV(i,iBinnedMol)+LevelEeV(1,iBinnedMol), DeltaEintDiss(i,iBinnedMol), LevelEeVRot(i,iBinnedMol), LevelEeVVib0(i,iBinnedMol), Levelg(i,iBinnedMol), LevelQ(i,iBinnedMol), ...
        dVIn(i,iBinnedMol), dVOut(i,iBinnedMol), rMin(i,iBinnedMol), VMin(i,iBinnedMol), rMax(i,iBinnedMol), VMax(i,iBinnedMol), Tau(i,iBinnedMol), Egam(i,iBinnedMol));
      end
      fclose(fileID);
    end
    
    
  end
    
end
