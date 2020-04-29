%% The Function reads the Rates from the HD5 File
%
%  Input Global Var: - Temp.TNowChar
%                    - Syst.HDF5_File
%
function Compute_Rates_Global()    

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
        
        
        iMol   = Syst.Pair(1).ToMol;
        iNBins = Syst.Molecule(iMol).EqNStatesIn;
        
        ExchTot = zeros(iNBins, Syst.NProc-2);
        for iExch = 1:Syst.NProc-2
            ExchTot(:,iExch) = sum(Rates.T(Temp.iT).ExchType(iExch).Exch, 2);
        end
        
        Rates.T(Temp.iT).DissGlobal = zeros(Kin.T(Temp.iT).NSteps,1);
        Rates.T(Temp.iT).ExchGlobal = zeros(Kin.T(Temp.iT).NSteps,Syst.NProc-2);
        
        for iStep = 1:Kin.T(Temp.iT).NSteps
            Rates.T(Temp.iT).DissGlobal(iStep)           = sum( Rates.T(Temp.iT).Diss(:,1) .* Kin.T(Temp.iT).Molecule(iMol).DF(iStep,:)' ); 
            for iExch = 1:Syst.NProc-2
                Rates.T(Temp.iT).ExchGlobal(iStep,iExch) = sum( ExchTot(:,iExch)           .* Kin.T(Temp.iT).Molecule(iMol).DF(iStep,:)' ); 
            end
        end
        
        fprintf('= Compute_Rates_Global ================= T = %i K\n', Temp.TNow)
        fprintf('====================================================\n')
        fprintf('Eq. Dissociation    Rate = %e cm^3/s\n', Rates.T(Temp.iT).DissGlobal(end) ) 
        for iExch = 1:Syst.NProc-2
            fprintf('Eq. Exchange (Nb %i) Rate = %e cm^3/s\n', iExch, Rates.T(Temp.iT).ExchGlobal(end,iExch) ) 
        end
        fprintf('====================================================\n\n')

        
    else
        
        
        iMol   = Syst.Pair(1).ToMol;
        jMol   = Syst.Pair(6).ToMol;

        Rates.T(Temp.iT).DissGlobal          = zeros(Kin.T(Temp.iT).NSteps,1);
        Rates.T(Temp.iT).DissGlobal_Diss     = zeros(Kin.T(Temp.iT).NSteps,1);
        Rates.T(Temp.iT).DissGlobal_DissInel = zeros(Kin.T(Temp.iT).NSteps,1);
        
        for iStep = 1:Kin.T(Temp.iT).NSteps

            for iBin = 1:Syst.Molecule(iMol).EqNStatesIn
                for jBin = iBin:Syst.Molecule(jMol).EqNStatesIn
                    Rates.T(Temp.iT).DissGlobal_Diss(iStep)     = Rates.T(Temp.iT).DissGlobal_Diss(iStep)     + 2.0 * Rates.T(Temp.iT).Diss(iBin,jBin,1)            * Kin.T(Temp.iT).Molecule(iMol).DF(iStep,iBin) * Kin.T(Temp.iT).Molecule(jMol).DF(iStep,jBin); 
                end
            end

            for iBin = 1:Syst.Molecule(iMol).EqNStatesIn
                for jBin = iBin:Syst.Molecule(jMol).EqNStatesIn
                    Rates.T(Temp.iT).DissGlobal_DissInel(iStep) = Rates.T(Temp.iT).DissGlobal_DissInel(iStep) + sum(sum(Rates.T(Temp.iT).DissInel(iBin,jBin,:,1)))  * Kin.T(Temp.iT).Molecule(iMol).DF(iStep,iBin) * Kin.T(Temp.iT).Molecule(jMol).DF(iStep,jBin); 
                end
            end

            Rates.T(Temp.iT).DissGlobal(iStep)                  = Rates.T(Temp.iT).DissGlobal_Diss(iStep)     + Rates.T(Temp.iT).DissGlobal_DissInel(iStep);

        end
        
        fprintf('= Compute_Rates_Global ================= T = %i K\n', Temp.TNow)
        fprintf('====================================================\n')
        fprintf('Eq. Dissociation    Rate = %e cm^3/s\n', Rates.T(Temp.iT).DissGlobal(end) ) 
%         for iExch = 1:Syst.NProc-2
%             fprintf('Eq. Exchange (Nb %i) Rate = %e cm^3/s\n', iExch, Rates.T(Temp.iT).ExchGlobal(end,iExch) ) 
%         end
        fprintf('====================================================\n\n')
        
        
    end
    
    
end