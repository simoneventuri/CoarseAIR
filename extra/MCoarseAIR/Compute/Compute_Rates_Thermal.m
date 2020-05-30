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
    
    global Rates Syst Temp Input
    

    fprintf('= Compute_Rates_Thermal ================ T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')    

    
    if (Syst.NAtoms == 3)
    
        
        iMol   = Syst.Pair(1).ToMol;
        iNBins = Syst.Molecule(iMol).EqNStatesIn;
        iQTemp = Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q ./ Syst.Molecule(iMol).T(Temp.iT).GroupsIn.QTot;

        ExchTot = zeros(iNBins, Syst.NProc-2);
        for iExch = 1:Syst.NProc-2
            ExchTot(:,iExch) = sum(Rates.T(Temp.iT).ExchType(iExch).Exch, 2);
        end
        
        Rates.T(Temp.iT).DissTh = 0.0;
        Rates.T(Temp.iT).ExchTh = zeros(1,Syst.NProc-2);
        for iBin = 1:Syst.Molecule(iMol).EqNStatesIn
            Rates.T(Temp.iT).DissTh              = sum( Rates.T(Temp.iT).Diss(:,1) .* iQTemp(:) ); 
            for iExch = 1:Syst.NProc-2
                Rates.T(Temp.iT).ExchTh(1,iExch) = sum( ExchTot(:,iExch)           .* iQTemp(:) ); 
            end
        end
        
        
        fprintf('Eq. Dissociation    Rate = %e cm^3/s\n',  Rates.T(Temp.iT).DissTh )
        for iExch = 1:Syst.NProc-2
            fprintf('Eq. Exchange (Nb %i) Rate = %e cm^3/s\n', iExch, Rates.T(Temp.iT).ExchTh(iExch) ) 
        end

        
        %% Writing Dissociation and Exchange Values at Equilibrium and QSS 
        %
        [status,msg,msgID] = mkdir(Input.Paths.SaveDataFldr);
        FileName           = strcat(Input.Paths.SaveDataFldr, '/KEq.csv');
        if exist(FileName, 'file')
            fileID1  = fopen(FileName,'a');
        else
            fileID1  = fopen(FileName,'w');
            if Syst.NProc == 3
                HeaderStr = strcat('# T [K], K^D Eq, K_{', Syst.Molecule(Syst.ExchToMol(1)).Name, '}^E Eq \n');
            else
                HeaderStr = strcat('# T [K], K^D Eq, K_{', Syst.Molecule(Syst.ExchToMol(1)).Name, '}^E Eq, K_{', Syst.Molecule(Syst.ExchToMol(2)).Name, '}^E Eq \n');
            end
            fprintf(fileID1,HeaderStr);
        end
        if Syst.NProc == 3
            fprintf(fileID1,'%e,%e,%e\n',    Temp.TNow, Rates.T(Temp.iT).DissTh, Rates.T(Temp.iT).ExchTh(1,1) );
        else
            fprintf(fileID1,'%e,%e,%e,%e\n', Temp.TNow, Rates.T(Temp.iT).DissTh, Rates.T(Temp.iT).ExchTh(1,1), Rates.T(Temp.iT).ExchTh(1,2) );
        end
        fclose(fileID1);
        
        
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
        
        
        fprintf('Eq. Dissociation    Rate = %e cm^3/s\n',  Rates.T(Temp.iT).DissTh )
%         for iExch = 1:Syst.NProc-2
%             fprintf('Eq. Exchange (Nb %i) Rate = %e cm^3/s\n', iExch, Rates.T(Temp.iT).ExchTh(iExch) ) 
%         end
        
        
    end
    
    
    fprintf('====================================================\n\n')        

end