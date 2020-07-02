%% The Function reads the Rates from the File for PLATO
%
%  Input Global Var: - Temp.TNowChar
%                    - Syst.HDF5_File
%
function [Rates] = Read_Rates_FromPLATO(Rates, Syst, RatesFldr)    

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
    
    global Input Temp Param  
    
    
    fprintf('  = Read_Rates_FromPLATO ================= T = %i K\n', Temp.TNow)
    fprintf('  ====================================================\n')
    fprintf(['  Reading Rates in PLATO Format from Folder: ' RatesFldr '\n'])
    
    iMol=1;
    
    if (Syst.NAtoms == 3)
        
        if (Input.Kin.Proc.DissFlg > 0)
            if (Input.Kin.Proc.DissFlg == 1)
               FileName = strcat(RatesFldr, '/', 'Diss.dat');
            elseif (Input.Kin.Proc.DissFlg == 2)
               FileName = strcat(RatesFldr, '/', 'Diss_Corrected.dat');
            end
            fprintf(['  Reading Dissociation Rates from: ' FileName '\n'] )
            opts = delimitedTextImportOptions("NumVariables", 14);
            opts.DataLines = [1, Inf];
            opts.Delimiter = ["(", ")", "+", ",", ":"];
            opts.VariableNames = ["Var1", "VarName2", "Var3", "Var4", "Var5", "Var6", "e12", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14"];
            opts.SelectedVariableNames = ["VarName2", "e12"];
            opts.VariableTypes = ["string", "double", "string", "string", "string", "string", "double", "string", "string", "string", "string", "string", "string", "string"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            opts = setvaropts(opts, ["Var1", "Var3", "Var4", "Var5", "Var6", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14"], "WhitespaceRule", "preserve");
            opts = setvaropts(opts, ["Var1", "Var3", "Var4", "Var5", "Var6", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14"], "EmptyFieldRule", "auto");
            tbl = readtable(FileName, opts);
            iIdx      = tbl.VarName2;
            RatesTemp = tbl.e12;
            clear opts tbl
            for iProc = 1:length(RatesTemp)
                Rates.T(Temp.iT).Diss(iIdx(iProc)) = RatesTemp(iProc);
            end
            clear iIdx RatesTemp
        end
        
        FileName = strcat(RatesFldr, '/Inel.dat');
        fprintf(['  Reading Inelastic Rates from: ' FileName '\n'] )
        opts = delimitedTextImportOptions("NumVariables", 15);
        opts.DataLines = [1, Inf];
        opts.Delimiter = ["(", ")", "+", ",", ":"];
        opts.VariableNames = ["Var1", "VarName2", "Var3", "Var4", "VarName5", "Var6", "Var7", "e11", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"];
        opts.SelectedVariableNames = ["VarName2", "VarName5", "e11"];
        opts.VariableTypes = ["string", "double", "string", "string", "double", "string", "string", "double", "string", "string", "string", "string", "string", "string", "string"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts = setvaropts(opts, ["Var1", "Var3", "Var4", "Var6", "Var7", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, ["Var1", "Var3", "Var4", "Var6", "Var7", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"], "EmptyFieldRule", "auto");
        tbl = readtable(FileName, opts);
        iIdx      = tbl.VarName2;
        jIdx      = tbl.VarName5;
        RatesTemp = tbl.e11;
        clear opts tbl
        for iProc = 1:length(RatesTemp)
            Rates.T(Temp.iT).Inel(iIdx(iProc),jIdx(iProc)) = RatesTemp(iProc);
        end
        
        NLevels=size(Rates.T(Temp.iT).Inel,1);
        for jLevel=1:NLevels
            for iLevel=1:NLevels
                if (Syst.Molecule(iMol).T(Temp.iT).GroupsIn.EeV(iLevel) < Syst.Molecule(iMol).T(Temp.iT).GroupsIn.EeV(jLevel))
                    Rates.T(Temp.iT).Inel(iLevel,jLevel) = Rates.T(Temp.iT).Inel(jLevel,iLevel) * Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q(jLevel) / Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q(iLevel);
                end
            end
        end
        clear iIdx jIdx RatesTemp

        
        
        for iExch=1:size(Syst.ExchToMol,1)
            FileName = strcat(RatesFldr, '/Exch_Type', num2str(iExch), '.dat');
            fprintf(['  Reading Exchange Rates Nb ' num2str(iExch) ' from: ' FileName '\n'] )
            jMol = Syst.ExchToMol(iExch);
            
            opts = delimitedTextImportOptions("NumVariables", 15);
            opts.DataLines = [1, Inf];
            opts.Delimiter = ["(", ")", "+", ",", ":"];
            opts.VariableNames = ["Var1", "VarName2", "Var3", "Var4", "VarName5", "Var6", "Var7", "e11", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"];
            opts.SelectedVariableNames = ["VarName2", "VarName5", "e11"];
            opts.VariableTypes = ["string", "double", "string", "string", "double", "string", "string", "double", "string", "string", "string", "string", "string", "string", "string"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            opts = setvaropts(opts, ["Var1", "Var3", "Var4", "Var6", "Var7", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"], "WhitespaceRule", "preserve");
            opts = setvaropts(opts, ["Var1", "Var3", "Var4", "Var6", "Var7", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"], "EmptyFieldRule", "auto");
            tbl = readtable(FileName, opts);
            iIdx      = tbl.VarName2;
            jIdx      = tbl.VarName5;
            RatesTemp = tbl.e11;
            clear opts tbl
            for iProc = 1:length(RatesTemp)
                Rates.T(Temp.iT).ExchType(iExch).Exch(iIdx(iProc),jIdx(iProc)) = RatesTemp(iProc);
            end
            
            iNLevels=size(Rates.T(Temp.iT).ExchType(iExch).Exch,1);
            jNLevels=size(Rates.T(Temp.iT).ExchType(iExch).Exch,2);
            fprintf(['  Matrix Size: (' num2str(iNLevels) ',' num2str(jNLevels) ')\n'] )
            for jLevel=1:jNLevels
                for iLevel=1:iNLevels
                    if (Syst.Molecule(iMol).T(Temp.iT).GroupsIn.EeV(iLevel) < Syst.Molecule(jMol).T(Temp.iT).GroupsIn.EeV(jLevel)) && (max(iLevel,jLevel)<=min(iNLevels,jNLevels))
                        Rates.T(Temp.iT).ExchType(iExch).Exch(iLevel,jLevel) = Rates.T(Temp.iT).ExchType(iExch).Exch(jLevel,iLevel) * Syst.Molecule(jMol).T(Temp.iT).GroupsIn.Q(jLevel) / Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q(iLevel);
                    end
                end
            end
            clear iIdx jIdx RatesTemp
        end
        
        
    else

        
    end
    
    
    fprintf('  ====================================================\n\n')

    
end