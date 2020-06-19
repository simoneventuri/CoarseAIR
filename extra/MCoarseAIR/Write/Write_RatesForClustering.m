%% The Function writes the Rates in the Format for being read by Amal's Clustering Algorithm
%
function Write_RatesForClustering(Controls)

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

    global Input Kin Param Syst Temp Rates

    
    fprintf('= Write_RatesForClustering ============= T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')

    
    WriteFldr = strcat(Controls.WriteFldr);
    [status,msg,msgID] = mkdir(WriteFldr);
    WriteFldr = strcat(WriteFldr, '/', Syst.NameLong, '/');
    [status,msg,msgID] = mkdir(WriteFldr);
    WriteFldr = strcat(WriteFldr, '/T', Temp.TNowChar, 'K/');
    [status,msg,msgID] = mkdir(WriteFldr);

    cm3_To_m3 = 1.d-6;
    MinRate   = Controls.MinRate .* cm3_To_m3;         %%% For Cutting out Rates < MinRate
    fprintf('Cutting Value for Rates = %e m3/s\n', MinRate)


    for iMol = 1:1
        fprintf(['Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );
        MinState = max(Controls.MinState, 1);
        MaxState = min(Controls.MaxState, Syst.Molecule(iMol).NLevels);

        
        
        fprintf('Writing Level Properties\n')
        FileName1 = strcat(WriteFldr, '/LevelsInfo.dat');
        fileID1   = fopen(FileName1,'w');
        jLevel    = 0;
        Mapping   = -1. * ones(Syst.Molecule(iMol).NLevels,1);
        for iLevel = MinState:MaxState
            WritingFlg = ( (sum(Rates.T(Temp.iT).Inel(iLevel,:)) > Controls.MinRate*Syst.Molecule(iMol).NLevels) && ...
                           (sum(Rates.T(Temp.iT).Inel(:,iLevel)) > Controls.MinRate*Syst.Molecule(iMol).NLevels) );
            
            if (WritingFlg)
                Mapping(iLevel) = jLevel;
                fprintf(fileID1,'%i %e %e %e %i %i %e %i %i\n',    jLevel,                                      ...
                                                                   Syst.Molecule(iMol).LevelEeV(iLevel),        ...
                                                                   Syst.Molecule(iMol).Levelg(iLevel),          ...
                                                                   Syst.Molecule(iMol).LevelrIn(iLevel),        ...
                                                                   Syst.Molecule(iMol).Levelvqn(iLevel),        ...
                                                                   Syst.Molecule(iMol).Leveljqn(iLevel),        ...
                                                                   Syst.Molecule(iMol).LevelECB(iLevel),        ...
                                                                   Syst.Molecule(iMol).LevelToGroupOut(iLevel), ...
                                                                   jLevel);
                jLevel = jLevel+1;
            else
                fprintf('Excluded Level Nb %i\n', iLevel)
            end
            
        end
        fclose(fileID1);
        
        
        KDiss = Rates.T(Temp.iT).Diss(:,1)   .* cm3_To_m3;
        KRec  = Rates.T(Temp.iT).Recomb(:,1) .* cm3_To_m3;
        
        fprintf('Writing Dissociation Rate Coefficients\n')
        FileName2 = strcat(WriteFldr, '/DissRates.dat');
        fileID2    = fopen(FileName2,'w');
        for iLevel = MinState:MaxState
            if (KDiss(iLevel) > MinRate) && (Mapping(iLevel) > -1)
              fprintf(fileID2,'%i %e %e\n', Mapping(iLevel), KDiss(iLevel), KRec(iLevel) );
            end
        end
        fclose(fileID2);

        
        
        Kij   = Rates.T(Temp.iT).Inel;
        Kij   = Kij .* cm3_To_m3;
        Kji   = Kij'; 

        tic;
        fprintf('Writing Inelastic Rate Coefficients\n')
        FileName3 = strcat(WriteFldr, '/InelRates_WOExch.dat');
        fileID3    = fopen(FileName3,'w');
        for jLevel = MinState:MaxState
            for iLevel = MinState:jLevel-1
                if (Kji(iLevel,jLevel) > MinRate) && (Mapping(iLevel) > -1) && (Mapping(jLevel) > -1)
                    fprintf(fileID3,'%i %i %e %e\n', Mapping(iLevel), Mapping(jLevel), Kji(iLevel,jLevel), Kij(iLevel,jLevel)); 
                                                   %        i,        j,                Kij,                Kji 
                                                   %        where dni/dt  = - sum_j Kij ni nC + sum_j kji nj nC
                end
            end
        end
        fclose(fileID3);
        timee = toc;
        fprintf('Inelastic Rates Matrix written in %e s\n', timee)
        
        
        
        AddedInelFlg = false;
        Kij          = Rates.T(Temp.iT).Inel;
        for iExch = 1:Syst.NProc-2
            jMol = Syst.ExchToMol(iExch);
            if (jMol==1) && (Controls.IncludeExch)
                fprintf('Adding homogeneous Exchange to the Inelastic Processes\n')
                Kij          = Kij + Rates.T(Temp.iT).ExchType(iExch).Exch;
                AddedInelFlg = true;
            end
        end
        Kij   = Kij .* cm3_To_m3;
        Kji   = Kij'; 
    
        if (AddedInelFlg)
            tic;
            fprintf('Writing Inelastic + Exchange Rate Coefficients\n')
            FileName3 = strcat(WriteFldr, '/InelRates_WExch.dat');
            fileID3    = fopen(FileName3,'w');
            for jLevel = MinState:MaxState
                for iLevel = MinState:jLevel-1
                    if (Kji(iLevel,jLevel) > MinRate) && (Mapping(iLevel) > -1)  && (Mapping(jLevel) > -1)
                        fprintf(fileID3,'%i %i %e %e\n', Mapping(iLevel), Mapping(jLevel), Kji(iLevel,jLevel), Kij(iLevel,jLevel)); 
                                                       %        i,        j,                Kij,                Kji 
                                                       %        where dni/dt  = - sum_j Kij ni nC + sum_j kji nj nC
                    end
                end
            end
            fclose(fileID3);
            timee = toc;
            fprintf('Inelastic + Exchange Rates Matrix written in %e s\n', timee)
        end

        
        
    end

    
    fprintf('====================================================\n\n')

end