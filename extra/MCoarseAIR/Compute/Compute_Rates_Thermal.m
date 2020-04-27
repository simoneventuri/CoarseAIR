%% The Function reads the Rates from the HD5 File
%
%  Input Global Var: - Temp.TNowChar
%                    - Syst.HDF5_File
%
function Compute_Rates_Thermal()    

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
    
    global Rates Syst Temp Kin
    
    
    if (Syst.NAtoms == 3)
    
        
        
    else 

        iMol   = Syst.Pair(1).ToMol;
        jMol   = Syst.Pair(6).ToMol;

        iQTemp = Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q ./ Syst.Molecule(iMol).T(Temp.iT).GroupsIn.QTot;
        jQTemp = Syst.Molecule(jMol).T(Temp.iT).GroupsIn.Q ./ Syst.Molecule(jMol).T(Temp.iT).GroupsIn.QTot;

        Rates.T(Temp.iT).DissTh          = 0.0;
        Rates.T(Temp.iT).DissTh_Diss     = 0.0;
        Rates.T(Temp.iT).DissTh_DissInel = 0.0;

        for iBin = 1:Syst.Molecule(iMol).EqNStatesIn
            for jBin = iBin:Syst.Molecule(jMol).EqNStatesIn
                Rates.T(Temp.iT).DissTh_Diss     = Rates.T(Temp.iT).DissTh_Diss     + 2.0 * Rates.T(Temp.iT).Diss(iBin,jBin,1)           * iQTemp(iBin) * jQTemp(jBin); 
            end
        end

        for iBin = 1:Syst.Molecule(iMol).EqNStatesIn
            for jBin = iBin:Syst.Molecule(jMol).EqNStatesIn
                Rates.T(Temp.iT).DissTh_DissInel = Rates.T(Temp.iT).DissTh_DissInel + sum(sum(Rates.T(Temp.iT).DissInel(iBin,jBin,:,1))) * iQTemp(iBin) * jQTemp(jBin); 
            end
        end
        
        Rates.T(Temp.iT).DissTh = Rates.T(Temp.iT).DissTh_Diss + Rates.T(Temp.iT).DissTh_DissInel;
        
        Rates.T(Temp.iT).DissTh
    end

end