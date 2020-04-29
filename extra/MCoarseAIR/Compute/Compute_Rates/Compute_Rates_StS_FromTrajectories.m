%% The Function reads the Rates from the HD5 File
%
%  Input Global Var: - Temp.TNowChar
%                    - Syst.HDF5_File
%
function Compute_Rates_StS_FromTrajectories()    

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
    
    global Rates Syst Temp Kin Input
    
    
    opts = delimitedTextImportOptions("NumVariables", 4);
    opts.DataLines = [2, Inf];
    opts.Delimiter = " ";
    opts.VariableNames = ["vv", "Var2", "Var3", "Var4"];
    opts.SelectedVariableNames = "vv";
    opts.VariableTypes = ["double", "string", "string", "string"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";
    opts = setvaropts(opts, ["Var2", "Var3", "Var4"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var2", "Var3", "Var4"], "EmptyFieldRule", "auto");
    FileName = strcat(Input.Paths.ToQCTFldr, '/Velocity_', Temp.TNowChar, '.0.dat');
    tbl = readtable("/home/venturi/WORKSPACE/CoarseAIR/CO2_ALL/Test/Velocity_20000.0.dat", opts);
    vv  = tbl.vv;
    clear opts tbl
    
    
    if (Syst.NAtoms == 3)
    
        
        
    else 
        
        NMol    = Syst.NMolecules;
        iMol    = Syst.Pair(1).ToMol;
        jMol    = Syst.Pair(6).ToMol;
        iNBins  = Syst.Molecule(iMol).EqNStatesIn;
        jNBins  = Syst.Molecule(iMol).EqNStatesIn;    

        iProc = 1
        for iBin = 1:iNBins
            for jBin = 1:jNBins
                if (jBin >= iBin)
                    fprintf('i = %i; j = %i\n', iBin, jBin)
     

                    opts = delimitedTextImportOptions("NumVariables", 14);
                    opts.DataLines = [2, Inf];
                    opts.Delimiter = ",";
                    opts.VariableNames = ["Var1", "Var2", "Var3", "b_i", "j1_i", "v1_i", "j2_i", "v2_i", "arr_i", "j1_f", "v1_f", "j2_f", "v2_f", "arr_f"]
                    opts.SelectedVariableNames = ["b_i", "j1_i", "v1_i", "j2_i", "v2_i", "arr_i", "j1_f", "v1_f", "j2_f", "v2_f", "arr_f"];
                    opts.VariableTypes = ["string", "string", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
                    opts.ExtraColumnsRule = "ignore";
                    opts.EmptyLineRule = "read";
                    opts = setvaropts(opts, ["Var1", "Var2", "Var3"], "WhitespaceRule", "preserve");
                    opts = setvaropts(opts, ["Var1", "Var2", "Var3"], "EmptyFieldRule", "auto");
                    FileName = strcat(Input.Paths.ToQCTFldr, '/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Bins_', num2str(iBin), '_', num2str(jBin), '/trajectories.csv');
                    tbl      = readtable(FileName, opts);
                    b_i   = tbl.b_i;
                    j1_i  = tbl.j1_i;
                    v1_i  = tbl.v1_i;
                    j2_i  = tbl.j2_i;
                    v2_i  = tbl.v2_i;
                    arr_i = tbl.arr_i;
                    j1_f  = round(tbl.j1_f  - 0.5);
                    v1_f  = round(tbl.v1_f  - 0.5);
                    j2_f  = round(tbl.j2_f  - 0.5);
                    v2_f  = round(tbl.v2_f  - 0.5);
                    arr_f = round(tbl.arr_f - 0.5);
                    clear opts tbl

                    TotTraj       = zeros(Nb,1); 
                    DissVecTemp_1 = zeros(Nb,1); 
                    DissMatTemp_1 = zeros(NBins, Nb); 
                    InelMatTemp_1 = zeros(NBins, NBins, Nb); 
                    ExchMatTemp_1 = zeros(NBins, NBins, Nb); 
                    for iTraj=1:length(b_i)
                        ib          = sum(b_i(iTraj) >= bVec);
                        TotTraj(ib) = TotTraj(ib) + 1;

                        iP = floor(arr_f(iTraj) / 16.0);
                        kk =   mod(arr_f(iTraj) , 16.0);
                        ii =   mod(          kk ,  4.0);
                        jj = floor(          kk /  4.0);
                        i  = v1_f(iTraj) + 1;
                        j  = v2_f(iTraj) + 1;

                        if (ii>1) && (jj>1)
                            DissVecTemp_1(ib,1)   = DissVecTemp_1(ib,1)   + 1;
                        elseif (ii>1)
                            DissMatTemp_1(j,ib)   = DissMatTemp_1(j,ib)   + 1;
                        elseif (jj>1)
                            DissMatTemp_1(i,ib)   = DissMatTemp_1(i,ib)   + 1;
                        elseif (iP==1)
                            InelMatTemp_1(i,j,ib) = InelMatTemp_1(i,j,ib) + 1;
                        else
                            ExchMatTemp_1(i,j,ib) = ExchMatTemp_1(i,j,ib) + 1;
                        end
                    end

                    Rates.T(Temp.iT).Diss(iBin,jBin,1)       = sum(DissVecTemp_1 ./ TotTraj .* bVecRing) .* vv;
                    for kBin = 1:iNBins
                        Rates.T(Temp.iT).DissInel(iBin,jBin,kBin,lBin)    = sum(squeeze(DissMatTemp_1(kBin,:))'     ./ TotTraj .* bVecRing) .* vv;
                        for lBin = 1:jNBins
                            Rates.T(Temp.iT).InelMat(iBin,jBin,kBin,lBin) = sum(squeeze(InelMatTemp_1(kBin,lBin,:)) ./ TotTraj .* bVecRing) .* vv;
                            Rates.T(Temp.iT).ExchMat(iBin,jBin,kBin,lBin) = sum(squeeze(ExchMatTemp_1(kBin,lBin,:)) ./ TotTraj .* bVecRing) .* vv;
                        end
                    end


                end
                iProc=iProc+1;
            end
        end

    end
    
    
end