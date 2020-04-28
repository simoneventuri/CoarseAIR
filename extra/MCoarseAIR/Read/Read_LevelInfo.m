%% The Function Reads the Molecules' Level Info from the list used for QCT
%
%  Required Variables: - Input.SystNameLong
%        
function ReadLevelInfo()
      
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

    
    global Input Syst Param
  
    
    for iMol = 1:Syst.NMolecules       
    
        
        %% Reading levels_cut.inp
        filename = strcat(Input.Paths.ToQCTFldr, '/', Syst.Name, '/', Syst.Molecule(iMol).Name, '/levels_cut.inp')
        startRow = 16;
        formatSpec = '%6f%5f%15f%15f%15f%15f%15f%15f%15f%15f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileID);
        Syst.Molecule(iMol).EEhDiss   = min(dataArray{:, 8});
        
        Syst.Molecule(iMol).Levelvqn  = dataArray{:, 1};
        Syst.Molecule(iMol).Leveljqn  = dataArray{:, 2};
        Syst.Molecule(iMol).LevelEEh  = dataArray{:, 3} - Syst.Molecule(iMol).EEhDiss;
        Syst.Molecule(iMol).LevelEgam = dataArray{:, 4};
        Syst.Molecule(iMol).LevelrMin = dataArray{:, 5};
        Syst.Molecule(iMol).LevelrMax = dataArray{:, 6};
        Syst.Molecule(iMol).LevelVMin = dataArray{:, 7} - Syst.Molecule(iMol).EEhDiss;
        Syst.Molecule(iMol).LevelVMax = dataArray{:, 8} - Syst.Molecule(iMol).EEhDiss;
        Syst.Molecule(iMol).LevelTau  = dataArray{:, 9};
        Syst.Molecule(iMol).LevelrIn  = dataArray{:, 10};
        Syst.Molecule(iMol).LevelrOut = dataArray{:, 11};
        clearvars filename startRow formatSpec fileID dataArray ans;

        Syst.Molecule(iMol).NLevels   = length(Syst.Molecule(iMol).LevelEEh);
        Syst.Molecule(iMol).LevelEeV  = Syst.Molecule(iMol).LevelEEh * Param.EhToeV;
        Syst.Molecule(iMol).EEhRef    = min(Syst.Molecule(iMol).LevelEEh(1))
        Syst.Molecule(iMol).EeVRef    = min(Syst.Molecule(iMol).LevelEeV(1))
        Syst.Molecule(iMol).LevelEEh0 = Syst.Molecule(iMol).LevelEEh - min(Syst.Molecule(iMol).LevelEEh(1));
        Syst.Molecule(iMol).LevelEeV0 = Syst.Molecule(iMol).LevelEeV - min(Syst.Molecule(iMol).LevelEeV(1));
        Syst.Molecule(iMol).Nvqn      = max(Syst.Molecule(iMol).Levelvqn) + 1;
        Syst.Molecule(iMol).Njqn      = max(Syst.Molecule(iMol).Leveljqn) + 1;

        
        %% Computing Equivalent Number of States = Levels/Groups
        if strcmp(Syst.Molecule(iMol).KinMthdIn, 'StS')
            Syst.Molecule(iMol).EqNStatesIn = Syst.Molecule(iMol).NLevels;
        elseif strcmp(Syst.Molecule(iMol).KinMthdIn, 'VSM')
            Syst.Molecule(iMol).EqNStatesIn = Syst.Molecule(iMol).Nvqn;
        elseif strcmp(Syst.Molecule(iMol).KinMthdIn, 'CGM')
            Syst.Molecule(iMol).EqNStatesIn = Syst.Molecule(iMol).NGroupsIn;
        end
        
        Syst.NTotProc = 1;
        for iP = 1:3
           jP            = Param.iPOpp(iP);
           iMol          = Syst.Pair(iP).ToMol;
           jMol          = Syst.Pair(jP).ToMol;
           Syst.NTotProc = Syst.NTotProc + (Syst.Molecule(iMol).EqNStatesIn + 1) * (Syst.Molecule(jMol).EqNStatesIn + 1);
        end
        
        
        %% Reading QNsEnBin
        opts = delimitedTextImportOptions("NumVariables", 6);
        opts.DataLines = [2, Inf];
        opts.Delimiter = ",";
        opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Levelg", "LevelToBin"];
        opts.SelectedVariableNames = ["Levelg", "LevelToBin"];
        opts.VariableTypes = ["string", "string", "string", "string", "double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4"], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4"], "EmptyFieldRule", "auto");
        tbl = readtable(strcat(Input.Paths.ToQCTFldr, '/', Syst.Name, '/', Syst.Molecule(iMol).Name, '/Bins_', num2str(Syst.Molecule(iMol).EqNStatesIn), '/QNsEnBin.csv'), opts);
        Syst.Molecule(iMol).Levelg     = tbl.Levelg;
        Syst.Molecule(iMol).LevelToBin = tbl.LevelToBin;
        clear opts tbl
        
    end
    
end
