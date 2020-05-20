%% The Function Computes /Reads the Mapping Level -> Group
%        
function Compute_GroupsOutMapping()
      
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

    
    global Input Syst Rates Param Temp 

    fprintf('= Compute_GroupsOutMapping =========================\n')
    fprintf('====================================================\n')

    
    for iMol = 1:Syst.NMolecules       
        fprintf('Mapping Levels to Output Group for Molecule Nb %i \n',  iMol )

        if strcmp(Input.Kin.MolResolutionOut(iMol,:), 'VSM')
            fprintf('The Molecule will be Grouped Vibrationally Specific\n' )
            
            Syst.Molecule(iMol).LevelToGroupOut = Syst.Molecule(1).Levelvqn+1;
            
        else
            fprintf(strcat('The Molecule will be Grouped Based on the Mapping From File:', '{}', Input.Kin.PathToMappingOut(iMol,:), '\n') );

            opts = delimitedTextImportOptions("NumVariables", 2);
            opts.DataLines = [2, Inf];
            opts.Delimiter = ",";
            opts.VariableNames = ["Var1", "Bin"];
            opts.SelectedVariableNames = "Bin";
            opts.VariableTypes = ["string", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
            opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");
            tbl = readtable(Input.Kin.PathToMappingOut(iMol,:), opts);
            Syst.Molecule(iMol).LevelToGroupOut = tbl.Bin;
            clear opts tbl
    
        end
        
        Syst.Molecule(iMol).NGroupsOut = max(Syst.Molecule(iMol).LevelToGroupOut);
        
    end

    
    fprintf('====================================================\n\n')
    
end
