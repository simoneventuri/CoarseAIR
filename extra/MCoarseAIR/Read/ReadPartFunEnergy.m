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

function [QBins, EeV] = ReadPartFunEnergy(iT)     

    %% (METHOD DEPENDENT)

    global T0_Vec BinnedMolName NBinnedMol NBins DatabasePath OutputPath SystemPath KinMthd

    for iBinnedMol=1:NBinnedMol

        % Reading Binned Molecules' Bins Partion Functions and Bins Energies
        
        
        filename = strcat(SystemPath,'/',BinnedMolName(iBinnedMol,:),'/Bins_',num2str(NBins(iBinnedMol)),'/T',num2str(T0_Vec(iT)),'.csv')
        opts = delimitedTextImportOptions("NumVariables", 3);
        opts.DataLines = [2, Inf];
        opts.Delimiter = ",";
        opts.VariableNames = ["LevelBinRatio", "PartFunc", "EnergyeV"];
        opts.VariableTypes = ["double", "double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        tbl = readtable(filename, opts);
        Temp                          = length(tbl.LevelBinRatio);
        QBinsRatio(1:Temp,iBinnedMol) = tbl.LevelBinRatio;
        QBins(1:Temp,iBinnedMol)      = tbl.PartFunc;
        EeV(1:Temp,iBinnedMol)        = tbl.EnergyeV;
        clear opts tbl 
        
    end
    
end
