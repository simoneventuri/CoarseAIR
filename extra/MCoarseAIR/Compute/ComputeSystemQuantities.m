%% The Function computes the parameters of the chemical system
%
%  Input Global Var: - NBins          Nb of Levels/Bins
%                    - T0_Vec:        Vector of Translational Temperatures (e.g.: [10000])
%                    - NAtoms         Number of Atoms in the system
%                    - AtomsName      Vector of strings containing the names of the atoms
%                    - MoleculesName  Vector of strings containing the names of the molecules
%                       

function ComputeSystemQuantities()    

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

  global NAtoms AtomsName NBins MoleculesName T0_Vec 
  global xDim yDim TotDim PlotPair ProcToLevIP TempChar TempChar2 Pair_Name Pair_To_Molecule Pair_to_Atoms iPInternal iPExternal

  if NAtoms == 3
      Pair_to_Atoms(1,:)=[1,2];
      Pair_to_Atoms(2,:)=[1,3];
      Pair_to_Atoms(3,:)=[2,3];
  end


  xDim
  TotDim(1)=1;
  for iP=1:3
      if (AtomsName(Pair_to_Atoms(iP,1)) == AtomsName(Pair_to_Atoms(iP,2)))
          Pair_Name(iP,:)=strcat(AtomsName(Pair_to_Atoms(iP,1)),'2');
      else
          Pair_Name(iP,:)=strcat(AtomsName(Pair_to_Atoms(iP,1)),AtomsName(Pair_to_Atoms(iP,2)));
      end
      Pair_Name
      for iMolecules=1:size(MoleculesName,1)
          if sum((Pair_Name(iP,:) == MoleculesName(iMolecules,:))) == 2 || sum((Pair_Name(iP,2:-1:1) == MoleculesName(iMolecules,:))) == 2
              Pair_To_Molecule(iP)=iMolecules;
              TotDim(iP+1) = TotDim(iP) + NBins(iMolecules);
              yDim(iP)     = NBins(iMolecules);
              break
          else
              TotDim(iP+1) = TotDim(iP) + 1;
              yDim(iP)     = 1;
          end
      end
  end
  TotDim
  yDim


  PlotPair(1:3)=1;
  for iP=2:3
      for jP=iP-1,1
          if (Pair_Name(iP,:) == Pair_Name(jP,:))
              PlotPair(iP) = 0;
              break
          end
      end
  end
  PlotPair


  iPInternal = [1];
  iPExternal = [];
  for iP=2:3
    if PlotPair(iP) == 0
       iPInternal = [iPInternal, iP];
    else
       iPExternal = [iPExternal, iP];
    end
  end 
  iPInternal
  iPExternal


  iP               = 1;
  iLevel           = 1;
  ProcToLevIP      = zeros(TotDim(end),2);
  for iProcess = 2:TotDim(end)
    if iProcess > TotDim(iP+1)
      iP     = iP+1;
      iLevel = 1;
    end
    ProcToLevIP(iProcess,:) = [ iLevel, iP];
    iLevel                  = iLevel+1;
  end


  TempChar  = ['          ';'          ';'          ';'          ';'          ';];
  TempChar2 = ['          ';'          ';'          ';'          ';'          ';]; 
  TempChar(1,:)  = ['T = ', num2str(T0_Vec(1),'%10.5i'), 'K']
  TempChar2(1,:) = strcat(TempChar(1,:))
  for iT = 2:length(T0_Vec)
      TempChar(iT,:)  = ['T = ', num2str(T0_Vec(iT),'%10.5i'), 'K']
      TempChar2(iT,:) = strcat(TempChar(iT,:))
  end

end
