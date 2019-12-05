%% The Function checks if the rates files exist and if they are not empty
%
%  Input Arguments:  - iT:            Index for the current Translational Temperature
%                    - iBinsStart:    (Optional) Initial Level/Bin to read. If not given, it is assumed equal to be equal to 1
%                    - iBinsEnd:      (Optional) Final Level/Bin to read. If not given, it is assumed equal to be equal to NBins
%
%  Input Global Var: - T0_Vec:        Vector of Translational Temperatures (e.g.: [10000])
%                    - RatesPath:     The path to the output folder (e.g.: ../Test/N3/N2/Rates)
%                    - MoleculesName: A vector of strings containing the name of the molecules present in the system
%                    - NBins          Nb of Levels/Bins
%

function [] = CheckEmptyRates(iT, iBinsStart, iBinsEnd)    

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

  global System NPESs NBins T0_Vec RatesPath ProcToLevIP MoleculesName NMolecules NLevels

  for iMol = 1:NMolecules

    str = strcat('Molecule= ', MoleculesName(iMol,:),'; Temperature= ', num2str(T0_Vec(iT)));
    disp( str )

    for iLevel = iBinsStart:iBinsEnd

      for iPES = 1:NPESs

        filename = strcat(RatesPath, '/T_' ,num2str(T0_Vec(iT)), '_', num2str(T0_Vec(iT)), '/Bin',num2str(iLevel),'.dat');

        if exist(filename, 'file') == 2

          startRow = 6;
          formatSpec = '%*20s%20f%20f%20f%[^\n\r]';
          fileID = fopen(filename,'r');
          dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
          fclose(fileID);
          ProcVec = dataArray{:, 1};
          RateVec = dataArray{:, 2};
          SDVec = dataArray{:, 3};
          clearvars filename startRow formatSpec fileID dataArray ans;

          if length(ProcVec) < 1
            str = strcat('Molecule= ', MoleculesName(iMol,:),'; Temperature= ', num2str(T0_Vec(iT)), '; PES= ', num2str(iPES), '; iLevel= ', num2str(iLevel),':   NO RATES FOUND IN THE FILE!!!');
            disp( str )
          end

        else
          
          str = strcat('Molecule= ', MoleculesName(iMol,:),'; Temperature= ', num2str(T0_Vec(iT)), '; PES= ', num2str(iPES), '; iLevel= ', num2str(iLevel),':   NO RATES FILE FOUND !!!');
          disp( str )
          
        end

      end

    end

  end
  
end