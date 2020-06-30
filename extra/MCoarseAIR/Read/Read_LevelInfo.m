%% The Function Reads the Molecules' Level Info from the list used for QCT
%        
function [Syst] = Read_LevelInfo(Syst)
      
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

    
    global Input Param

    fprintf('= Read_LevelInfo ===================================\n')
    fprintf('====================================================\n')
    
    
    for iMol = 1:Syst.NMolecules       
        fprintf('Reading Level Quantities for Molecule Nb %i \n',  iMol )
        
        LevelsFile = strcat(Syst.HDF5_File);
        fprintf(['  Checking if HDF5 File is Present: ' LevelsFile '\n'] )
                
        if (1==2)%isfile(strcat(LevelsFile))
            
            %% Reading HDF5 File
            fprintf(['  Reading From File: ' LevelsFile '\n'] )
            
            VarChar                        = strcat('/', Syst.Molecule(iMol).Name, '/LevelVMax');
            TempVar                        = h5read(Syst.HDF5_File, VarChar);
            Syst.Molecule(iMol).EEhDiss    = min(TempVar);
            
            VarChar                        = strcat('/', Syst.Molecule(iMol).Name, '/Levelvqn');
            Syst.Molecule(iMol).Levelvqn   = h5read(Syst.HDF5_File, VarChar);            
            VarChar                        = strcat('/', Syst.Molecule(iMol).Name, '/Leveljqn');
            Syst.Molecule(iMol).Leveljqn   = h5read(Syst.HDF5_File, VarChar);
            VarChar                        = strcat('/', Syst.Molecule(iMol).Name, '/LevelEEh');
            Syst.Molecule(iMol).LevelEEh   = h5read(Syst.HDF5_File, VarChar) - Syst.Molecule(iMol).EEhDiss;            
            VarChar                        = strcat('/', Syst.Molecule(iMol).Name, '/LevelEgam');
            Syst.Molecule(iMol).LevelEgam  = h5read(Syst.HDF5_File, VarChar);            
            VarChar                        = strcat('/', Syst.Molecule(iMol).Name, '/LevelrMin');
            Syst.Molecule(iMol).LevelrMin  = h5read(Syst.HDF5_File, VarChar);            
            VarChar                        = strcat('/', Syst.Molecule(iMol).Name, '/LevelrMax');
            Syst.Molecule(iMol).LevelrMax  = h5read(Syst.HDF5_File, VarChar);            
            VarChar                        = strcat('/', Syst.Molecule(iMol).Name, '/LevelVMin');
            Syst.Molecule(iMol).LevelVMin  = h5read(Syst.HDF5_File, VarChar) - Syst.Molecule(iMol).EEhDiss;
            VarChar                        = strcat('/', Syst.Molecule(iMol).Name, '/LevelVMax');
            Syst.Molecule(iMol).LevelVMax  = h5read(Syst.HDF5_File, VarChar) - Syst.Molecule(iMol).EEhDiss;
            VarChar                        = strcat('/', Syst.Molecule(iMol).Name, '/LevelTau');
            Syst.Molecule(iMol).LevelTau   = h5read(Syst.HDF5_File, VarChar);
            VarChar                        = strcat('/', Syst.Molecule(iMol).Name, '/LevelrIn');
            Syst.Molecule(iMol).LevelrIn   = h5read(Syst.HDF5_File, VarChar);
            VarChar                        = strcat('/', Syst.Molecule(iMol).Name, '/LevelrOut');
            Syst.Molecule(iMol).LevelrOut  = h5read(Syst.HDF5_File, VarChar);
            
            VarChar                            = strcat('/', Syst.Molecule(iMol).Name, '/Levelg');
            Syst.Molecule(iMol).Levelg         = h5read(Syst.HDF5_File, VarChar);
            VarChar                            = strcat('/', Syst.Molecule(iMol).Name, '/LevelToGroupIn');
            Syst.Molecule(iMol).LevelToGroupIn = h5read(Syst.HDF5_File, VarChar);
            
        else
        
            %% Reading levels_cut.inp
            filename = strcat(Input.Paths.ToQCTFldr, '/', Syst.Name, '/', Syst.Molecule(iMol).Name, '/levels_cut.inp');
            fprintf(strcat('Reading From File: ', filename, ' \n'))  

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
            
            
            %% Reading QNsEnBin.csv
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
            Syst.Molecule(iMol).Levelg         = tbl.Levelg;
            Syst.Molecule(iMol).LevelToGroupIn = tbl.LevelToBin;
            clear opts tbl
            
        end
        
        Syst.Molecule(iMol).NLevels   = length(Syst.Molecule(iMol).LevelEEh);
        Syst.Molecule(iMol).LevelEeV  = Syst.Molecule(iMol).LevelEEh  * Param.EhToeV;
        Syst.Molecule(iMol).LevelVMin = Syst.Molecule(iMol).LevelVMin * Param.EhToeV;
        Syst.Molecule(iMol).LevelVMax = Syst.Molecule(iMol).LevelVMax * Param.EhToeV;
        Syst.Molecule(iMol).EEhRef    = min(Syst.Molecule(iMol).LevelEEh(1));
        Syst.Molecule(iMol).EeVRef    = min(Syst.Molecule(iMol).LevelEeV(1));
        Syst.Molecule(iMol).LevelEEh0 = Syst.Molecule(iMol).LevelEEh - min(Syst.Molecule(iMol).LevelEEh(1));
        Syst.Molecule(iMol).LevelEeV0 = Syst.Molecule(iMol).LevelEeV - min(Syst.Molecule(iMol).LevelEeV(1));
        Syst.Molecule(iMol).Nvqn      = max(Syst.Molecule(iMol).Levelvqn) + 1;
        Syst.Molecule(iMol).Njqn      = max(Syst.Molecule(iMol).Leveljqn) + 1;
        
        Syst.Molecule(iMol).LevelECB  = Syst.Molecule(iMol).LevelVMax - Syst.Molecule(iMol).LevelEeV; 

        vToLevel = zeros(Syst.Molecule(iMol).Nvqn,1);
        for iLevels = 1:Syst.Molecule(iMol).NLevels
            if (Syst.Molecule(iMol).Leveljqn(iLevels) == 0) 
                vToLevel(Syst.Molecule(iMol).Levelvqn(iLevels)+1,1)                   = iLevels; 
                Syst.Molecule(iMol).vEeVVib0(Syst.Molecule(iMol).Levelvqn(iLevels)+1) = Syst.Molecule(iMol).LevelEeV0(iLevels);         
            end
        end
                
        for iLevels = 1:Syst.Molecule(iMol).NLevels
            Syst.Molecule(iMol).LevelEeVVib0(iLevels) = Syst.Molecule(iMol).LevelEeV0(vToLevel(Syst.Molecule(iMol).Levelvqn(iLevels)+1));
            Syst.Molecule(iMol).LevelEeVRot(iLevels)  = Syst.Molecule(iMol).LevelEeV0(iLevels) - Syst.Molecule(iMol).LevelEeVVib0(iLevels);
        end
        
    end
    
    
    fprintf('Grouping Molecules in Input \n')
    Syst = Group_In(Syst)
    
    
    
    for iMol = 1:Syst.NMolecules       
        fprintf('Reading Level Quantities for Molecule Nb %i \n',  iMol )
        
        
        %% Computing Equivalent Number of States = Levels/Groups
        if strcmp(Syst.Molecule(iMol).KinMthdIn, 'StS')
            Syst.Molecule(iMol).EqNStatesIn = Syst.Molecule(iMol).NLevels;
        elseif strcmp(Syst.Molecule(iMol).KinMthdIn, 'VSM')
            Syst.Molecule(iMol).EqNStatesIn = Syst.Molecule(iMol).Nvqn;
        elseif strcmp(Syst.Molecule(iMol).KinMthdIn, 'CGM')
            Syst.Molecule(iMol).EqNStatesIn = Syst.Molecule(iMol).NGroupsIn;
        end
        
        
        %% Computing the Overall number of possible final processes
        Syst.NTotProc = 1;
        if Syst.NAtoms == 3
            for kP = 1:3
               kMol          = Syst.Pair(kP).ToMol;
               Syst.NTotProc = Syst.NTotProc + (Syst.Molecule(kMol).EqNStatesIn + 1);
            end
        else
            for kP = 1:3
               lP            = Param.iPOpp(kP);
               kMol          = Syst.Pair(kP).ToMol;
               lMol          = Syst.Pair(lP).ToMol;
               Syst.NTotProc = Syst.NTotProc + (Syst.Molecule(kMol).EqNStatesIn + 1) * (Syst.Molecule(lMol).EqNStatesIn + 1);
            end
        end
      
        
    end
    
    
    fprintf('====================================================\n\n')  
    
end
