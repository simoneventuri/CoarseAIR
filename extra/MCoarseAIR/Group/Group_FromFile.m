%% Group the Levels Based on Centrifugal Barrier
%
function [LevelToGroup] = Group_FromFile(Syst, Controls, iMol)    
    
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

    global Input Kin Param Temp Rates

    fprintf('  = Group_FromFile ===================================\n')
    fprintf('  ====================================================\n')
    
    fprintf(['  Reading from File: ', char(Controls.FilePath(iMol)), '\n'])

    opts = delimitedTextImportOptions("NumVariables", 2);
    opts.DataLines = [2, Inf];
    opts.Delimiter = [" ", ","];
    opts.VariableNames = ["Id", "Group"];
    opts.VariableTypes = ["double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    tbl   = readtable(char(Controls.FilePath(iMol)), opts);
    Idx   = tbl.Id;
    Group = tbl.Group;
    clear opts tbl
    
    LevelToGroup = zeros(Syst.Molecule(iMol).NLevels,1);
    for iLevels = 1:Syst.Molecule(iMol).NLevels
      LevelToGroup(Idx(iLevels)) = Group(iLevels);
    end    
    
    fprintf('  Nb Bins = %i \n', max(LevelToGroup) );
    
    fprintf('  ====================================================\n\n')
    
end