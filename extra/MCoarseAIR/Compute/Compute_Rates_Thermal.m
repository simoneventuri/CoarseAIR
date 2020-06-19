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
    
    global Rates Syst Temp Input OtherSyst
    

    fprintf('= Compute_Rates_Thermal ================ T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')    

    
    if (Syst.NAtoms == 3)
    
        
        for iMol = 1:length(Input.Kin.ReadOtherSyst)+1
            if iMol == 1
                ComputeFlg = true;
                NProc      = Syst.NProc;
            elseif (Input.Kin.ReadOtherSyst(iMol-1))
                ComputeFlg = true;
                NProc      = OtherSyst(iMol-1).Syst.NProc;
            end

            if (ComputeFlg)
                iNBins = Syst.Molecule(iMol).EqNStatesIn;
                iQTemp = Syst.Molecule(iMol).T(Temp.iT).GroupsIn.Q ./ Syst.Molecule(iMol).T(Temp.iT).GroupsIn.QTot;

                Rates.T(Temp.iT).Molecule(iMol).DissTh = 0.0;
                Rates.T(Temp.iT).Molecule(iMol).ExchTh = zeros(1,NProc-2);
                for iBin = 1:Syst.Molecule(iMol).EqNStatesIn
                    Rates.T(Temp.iT).Molecule(iMol).DissTh              = sum( Rates.T(Temp.iT).Molecule(iMol).Overall(:,1)       .* iQTemp(:) ); 
                    for iExch = 1:NProc-2
                        Rates.T(Temp.iT).Molecule(iMol).ExchTh(1,iExch) = sum( Rates.T(Temp.iT).Molecule(iMol).Overall(:,2+iExch) .* iQTemp(:) ); 
                    end
                end


                fprintf('Eq. Dissociation    Rate = %e cm^3/s\n',  Rates.T(Temp.iT).Molecule(iMol).DissTh )
                for iExch = 1:NProc-2
                    fprintf('Eq. Exchange (Nb %i) Rate = %e cm^3/s\n', iExch, Rates.T(Temp.iT).Molecule(iMol).ExchTh(iExch) ) 
                end


                %% Writing Dissociation and Exchange Values at Equilibrium and QSS 
                %
                [status,msg,msgID] = mkdir(Input.Paths.SaveDataFldr);
                FileName           = strcat(Input.Paths.SaveDataFldr, '/KEq_', Input.Kin.Proc.OverallFlg, '_', Syst.Molecule(iMol).Name , '.csv');
                if exist(FileName, 'file')
                    fileID1  = fopen(FileName,'a');
                else
                    fileID1  = fopen(FileName,'w');
                    if NProc == 2
                        HeaderStr = strcat('# T [K], K^D Eq\n');
                    elseif NProc == 3
                        HeaderStr = strcat('# T [K], K^D Eq, K_1^E Eq \n');
                    else
                        HeaderStr = strcat('# T [K], K^D Eq, K_1^E Eq, K_2^E Eq \n');
                    end
                    fprintf(fileID1,HeaderStr);
                end
                if NProc == 2
                    fprintf(fileID1,'%e,%e\n',       Temp.TNow, Rates.T(Temp.iT).Molecule(iMol).DissTh );
                elseif NProc == 3
                    fprintf(fileID1,'%e,%e,%e\n',    Temp.TNow, Rates.T(Temp.iT).Molecule(iMol).DissTh, Rates.T(Temp.iT).Molecule(iMol).ExchTh(1,1) );
                else
                    fprintf(fileID1,'%e,%e,%e,%e\n', Temp.TNow, Rates.T(Temp.iT).Molecule(iMol).DissTh, Rates.T(Temp.iT).Molecule(iMol).ExchTh(1,1), Rates.T(Temp.iT).Molecule(iMol).ExchTh(1,2) );
                end
                fclose(fileID1);
            
            end

        end
        
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