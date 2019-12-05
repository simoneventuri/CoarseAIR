%% The Function writes the rates in CVS Files
%
%  Input Arguments:  - NLevels:        Number of Levels/Bins
%                    - RatesMatrix:    Matrix of StS/BintoBin Rates K(iLevel,jLevel)
%                    - ProcessesRates: Matrix of Processes Rates K(iLevel,iProcess)
%                    - LevelEeV:       Levels/Bins' energy in eV, with the 0 corresponding to the dissociation energy
%
%  Input Global Var: - EeVSpace:       Flag 0/1. If =0 writes the levels indexes; if =1 it writes the levels energies in eV

function WriteCSVFiles(NLevels, RatesMatrix, ProcessesRates, LevelEeV)

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

  global System Pair_To_Molecule MoleculesName EeVSpace
  
  WritingDir         = strcat('./CGQCT-',System,'-RatesCSV/');
  [status,msg,msgID] = mkdir(WritingDir);
  
  
  RatesMatrixTemp = zeros(size(RatesMatrix,1),size(RatesMatrix,2),size(MoleculesName,1));
  for iP = 1:3
    
    clear xCol yCol zCol n m
    n                                       = NLevels(1)
    m                                       = NLevels(Pair_To_Molecule(iP))
    NLevelsTemp(Pair_To_Molecule(iP)) = m;
    RatesMatrixTemp(1:n,1:m,Pair_To_Molecule(iP)) = RatesMatrixTemp(1:n,1:m,Pair_To_Molecule(iP)) + RatesMatrix(1:n,1:m,iP,1);

    % xCol = zeros(n*(n-1)/2,1);
    % yCol = zeros(n*(n-1)/2,1);
    % zCol = zeros(n*(n-1)/2,1);
    iTot = 1;
    for i=1:n
      for j=1:m
        if LevelEeV(j,Pair_To_Molecule(iP)) < LevelEeV(i,1)
          zTemp      = RatesMatrix(i,j,iP,1);
          if zTemp > 1.d-15
            if EeVSpace == 1
              xCol(iTot) = LevelEeV(i,1);
              yCol(iTot) = LevelEeV(j,Pair_To_Molecule(iP));
            else
              xCol(iTot) = i;
              yCol(iTot) = j;
            end
            zCol(iTot) = zTemp;
            iTot       = iTot+1;
          end 
        end
      end
    end
    NameStr = strcat(WritingDir, MoleculesName(1,1:2), '_', MoleculesName(Pair_To_Molecule(iP),1:2), '_Pair', num2str(iP), '.csv');
    csvwrite(NameStr,[xCol',yCol',zCol'])
    
  end
  
  
  for iMol = 1:size(MoleculesName,1)
    
    clear xCol yCol zCol
    n    = NLevels(1)
    m    = NLevelsTemp(iMol)
    % xCol = zeros(n*(n-1)/2,1);
    % yCol = zeros(n*(n-1)/2,1);
    % zCol = zeros(n*(n-1)/2,1);
    iTot = 1;
    for i=1:n
      for j=1:m
        if LevelEeV(j,iMol) < LevelEeV(i,1)
          zTemp      = RatesMatrixTemp(i,j,iMol,1);
          if zTemp > 1.d-15
            if EeVSpace == 1
              xCol(iTot) = LevelEeV(i,1);
              yCol(iTot) = LevelEeV(j,iMol);
            else
              xCol(iTot) = i;
              yCol(iTot) = j;
            end
            zCol(iTot) = zTemp;
            iTot       = iTot+1;
          end 
        end
      end
    end
    NameStr = strcat(WritingDir, MoleculesName(1,1:2), '_', MoleculesName(iMol,1:2), '.csv');
    csvwrite(NameStr,[xCol',yCol',zCol'])
    
  end
  
  
  clear xCol yCol zCol
  n    = NLevels(1);
  % xCol = zeros(n*(n-1)/2,1);
  % yCol = zeros(n*(n-1)/2,1);
  % zCol = zeros(n*(n-1)/2,1);
  iTot = 1;
  for i=1:n
    zTemp      = ProcessesRates(i,1,1);
    if zTemp > 1.d-15
      if EeVSpace == 1
        xCol(iTot) = LevelEeV(i,1);
        yCol(iTot) = 0;
      else
        xCol(iTot) = i;
        yCol(iTot) = j;
      end
      zCol(iTot) = zTemp;
      iTot       = iTot+1;
    end 
  end
  NameStr = strcat(WritingDir, MoleculesName(1,1:2), '_DissRates.csv');
  csvwrite(NameStr,[xCol',yCol',zCol'])
  
  
end