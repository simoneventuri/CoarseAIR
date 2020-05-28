%% The Function Reads the Grouped Molecules' Partion Functions and Energies
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

    global Temp Param Syst
    %global Input
 
    fprintf('= Read_EeV_and_Q_CG ==================== T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    
    
    for iMol = 1:Syst.NMolecules      
        
        fprintf('Computing Grouped Quantities for Molecule Nb %i \n',  iMol )
        
%         %% Read T5000.csv
%         filename = strcat(Input.Paths.ToQCTFldr, '/', Syst.Name, '/', Syst.Molecule(iMol).Name, '/Bins_', num2str(Syst.Molecule(iMol).EqNStatesIn), '/T', Temp.TNowChar, '.csv')
%         opts = delimitedTextImportOptions("NumVariables", 3);
%         opts.DataLines = [2, Inf];
%         opts.Delimiter = ",";
%         opts.VariableNames = ["LevelBinRatio", "PartFunc", "EnergyeV"];
%         opts.VariableTypes = ["double", "double", "double"];
%         opts.ExtraColumnsRule = "ignore";
%         opts.EmptyLineRule = "read";
%         tbl = readtable(filename, opts);
%         Syst.Molecule(iMol).T(Temp.iT).GroupsIn.QRatio = tbl.LevelBinRatio;
%         Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q      = tbl.PartFunc;
%         Syst.Molecule(iMol).T(Temp.iT).GroupsIn.EeV    = tbl.EnergyeV + Syst.Molecule(iMol).EeVRef;
%         clear opts tbl 
   

        %% Compute Group Energy and Part Function
        Syst.Molecule(iMol).T(Temp.iT).Levelq = Syst.Molecule(iMol).Levelg .* exp( -  Syst.Molecule(iMol).LevelEeV .* Param.Ue ./ (Temp.TNow .* Param.UKb) );
        
        Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q   = zeros(Syst.Molecule(iMol).EqNStatesIn,1);
        Syst.Molecule(iMol).T(Temp.iT).GroupsIn.EeV = zeros(Syst.Molecule(iMol).EqNStatesIn,1);
        Syst.Molecule(iMol).T(Temp.iT).QRot         = 0.0;
        for iLevel = 1:Syst.Molecule(iMol).NLevels
            iBin = Syst.Molecule(iMol).LevelToGroupIn(iLevel);
            Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q(iBin)     = Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q(iBin)   + Syst.Molecule(iMol).T(Temp.iT).Levelq(iLevel);
            Syst.Molecule(iMol).T(Temp.iT).GroupsIn.EeV(iBin)   = Syst.Molecule(iMol).T(Temp.iT).GroupsIn.EeV(iBin) + Syst.Molecule(iMol).T(Temp.iT).Levelq(iLevel) * Syst.Molecule(iMol).LevelEeV(iLevel);          
            Syst.Molecule(iMol).T(Temp.iT).QRot                 = Syst.Molecule(iMol).T(Temp.iT).QRot               + Syst.Molecule(iMol).T(Temp.iT).Levelq(iLevel);
        end
        Syst.Molecule(iMol).T(Temp.iT).GroupsIn.EeV   = Syst.Molecule(iMol).T(Temp.iT).GroupsIn.EeV ./ Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q;
        
        if strcmp(Syst.Molecule(iMol).KinMthdIn, 'StS')
            Syst.Molecule(iMol).T(Temp.iT).GroupsIn.g = Syst.Molecule(iMol).Levelg;
        else
            Syst.Molecule(iMol).T(Temp.iT).GroupsIn.g = Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q(iBin);
        end
        
        Syst.Molecule(iMol).T(Temp.iT).GroupsIn.QTot   = sum(Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q);
        Syst.Molecule(iMol).T(Temp.iT).GroupsIn.QRatio = Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q   ./ Syst.Molecule(iMol).T(Temp.iT).GroupsIn.QTot;
        
        
    end
    
    
    fprintf('====================================================\n\n')

end
