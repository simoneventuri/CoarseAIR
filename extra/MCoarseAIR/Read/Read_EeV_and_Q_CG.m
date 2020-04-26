%% The Function Reads the Grouped Molecules' Partion Functions and Energies
%
%  Initializing Input Global Var: - Syst.Paths.ToReadFldr: The path to the output folder (e.g.: .../CoarseAIR/run_O3/Test/ )
%
function Read_EeV_and_Q_CG()     

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

    global Syst Input Temp

 
    for iMol = 1:Syst.NMolecules        
        
        
        %% Read T5000.csv
        filename = strcat(Input.Paths.ToQCTFldr, '/', Syst.Name, '/', Syst.Molecule(iMol).Name, '/Bins_', num2str(Syst.Molecule(iMol).EqNStatesIn), '/T', Temp.TNowChar, '.csv')
        opts = delimitedTextImportOptions("NumVariables", 3);
        opts.DataLines = [2, Inf];
        opts.Delimiter = ",";
        opts.VariableNames = ["LevelBinRatio", "PartFunc", "EnergyeV"];
        opts.VariableTypes = ["double", "double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        tbl = readtable(filename, opts);
        Syst.Molecule(iMol).GroupsIn.QRatio = tbl.LevelBinRatio;
        Syst.Molecule(iMol).GroupsIn.Q      = tbl.PartFunc;
        Syst.Molecule(iMol).GroupsIn.EeV    = tbl.EnergyeV + Syst.Molecule(iMol).EeVRef;
        clear opts tbl 
        
        
    end
    
end
