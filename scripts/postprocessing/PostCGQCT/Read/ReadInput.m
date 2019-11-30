%% The Function reads the input variables that have been used in the CoarseAIR simulations 
%
%  Input Global Var: - PathToOutput: The path to the output folder (e.g.: ../Test/ )
%

function ReadInput()

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
  
  global PathToOutput MoleculesName KinMthd

  clear MoleculesName
  clear KinMthd
  
  global System NTtra T0_Vec NTint TInt_Vec NMolecules StartBin FinalBin NBins NBinnedMol BinnedMolName

  filename = strcat(PathToOutput,'/InputForBash.inp');
  delimiter = ' ';
  formatSpec = '%s%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
  fclose(fileID);
  InputForBash = [dataArray{1:end-1}];
  clearvars filename delimiter formatSpec fileID dataArray ans;

  iTot = 1;
  System = char(InputForBash(iTot,:))

  iTot = iTot + 1;
  TranFlg = str2num(char(InputForBash(iTot,:)))
  
  iTot = iTot + 1;
  NTtra  = str2num(char(InputForBash(iTot,:)))
  for iTtra = 1: NTtra
    iTot = iTot + 1;
    T0_Vec(iTtra) = str2num(char(InputForBash(iTot,:)))
  end

  iTot = iTot + 1;
  NTint = str2num(char(InputForBash(iTot,:)))
  if NTint > 1 
    for iTtra = 1:NTint
      iTot = iTot + 1;
      TInt_Vec(iTint) = str2num(char(InputForBash(iTot,:)))
    end
  end

  iTot = iTot + 1;
  NMolecules = str2num(char(InputForBash(iTot,:)))
  for iMol = 1:NMolecules
    iTot = iTot + 1;
    MoleculesName(iMol,:) = char(InputForBash(iTot,:))
    iTot = iTot + 1;
    SortMthdTemp          = char(InputForBash(iTot,:))
    if SortMthdTemp(1:3) == 'Sta'
       KinMthd(iMol,:) = 'STS'
       iTot = iTot + 1;
       StartBin = str2num(char(InputForBash(iTot,:)))
       iTot = iTot + 1;
       FinalBin = str2num(char(InputForBash(iTot,:)))
       iTot = iTot + 1;
    elseif SortMthdTemp(1:3) == 'Vib'
       KinMthd(iMol,:) = 'VIB'
       iTot = iTot + 1;
       StartBin = str2num(char(InputForBash(iTot,:)))
       iTot = iTot + 1;
       FinalBin = str2num(char(InputForBash(iTot,:)))
       iTot = iTot + 1;
    elseif sum(SortMthdTemp(1:3) == 'RoV') == 3 || sum(SortMthdTemp(1:3) == 'Fro') == 3
       KinMthd(iMol,:) = 'CGM'
       iTot = iTot + 1;
       StartBin = str2num(char(InputForBash(iTot,:)))
       iTot = iTot + 1;
       FinalBin = str2num(char(InputForBash(iTot,:)))
       iTot = iTot + 1;
       NBins    = str2num(char(InputForBash(iTot,:)))
    end
  end
  
  BinnedMolName = MoleculesName;
  NBinnedMol    = NMolecules;
  
end
