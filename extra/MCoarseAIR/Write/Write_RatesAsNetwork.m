%% Group the Levels Based on Centrifugal Barrier
%
function Write_RatesAsNetwork(Controls)    
    
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

    
    fprintf('= Write_RatesAsNetwork ================= T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')

    
    WriteFldr = strcat(Controls.WriteFldr);
    [status,msg,msgID] = mkdir(WriteFldr);
    WriteFldr = strcat(WriteFldr, '/', Syst.NameLong, '/');
    [status,msg,msgID] = mkdir(WriteFldr);
    WriteFldr = strcat(WriteFldr, '/T', Temp.TNowChar, 'K/');
    [status,msg,msgID] = mkdir(WriteFldr);
    
    
    MinState = max(1,Controls.MinState);
    MaxState = min(Syst.Molecule(1).NLevels,Controls.MaxState);
    fprintf('Writing Rates for Levels %i-to-%i in a Network File\n', MinState, MaxState)

    cm3_To_m3 = 1.d-6;
    MinRate   = Controls.MinRate .* cm3_To_m3;         %%% For Cutting out Rates < MinRate
    fprintf('Cutting Value for Rates = %e m3/s\n', MinRate)
    
    
    MinState = max(1,Controls.MinState);
    MaxState = min(Syst.Molecule(1).NLevels,Controls.MaxState);
    fprintf('Writing Rates for Levels %i-to-%i in a Network File\n', MinState, MaxState)

    cm3_To_m3 = 1.d-6;
    MinRate   = Controls.MinRate .* cm3_To_m3;         %%% For Cutting out Rates < MinRate
    fprintf('Cutting Value for Rates = %e m3/s\n', MinRate)
    

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
    Kij   = Kij - diag(diag(Kij));
    WOut  = sum(Kij,2);
    Kij   = Kij - diag(WOut);
    Kji   = Kij'; 
    
    

    P0    = 9940.737d0;
    X0    = [0.05d0, 0.95d0];
    T0    = 300;

    V0      = 1.d0;
    R       = 8.3144598;
    n0      = (P0 * V0 / (Param.Ru * T0)) * Param.AvN;
    nO0     = n0 * X0(1);
    nCOTot0 = n0 * X0(2);

    TempMat = Kji .* nO0;      
    NEdges  = sum(sum(abs(Kji) > MinRate));
    fprintf('Found %e Edges\n', NEdges)


    NLevelTemp = 0;
    LevelToNew = [];
    for i = MinState:MaxState
        NLevelTemp    = NLevelTemp + 1;
        LevelToNew(i) = NLevelTemp;
    end


    tic
    FileName = strcat(WriteFldr, './InelRates.net');
    fileID   = fopen(FileName,'w');
    fprintf(fileID,'# A network in Pajeks .net format\n');

    fprintf(fileID,'*Vertices %i\n', NLevelTemp);
    for i = MinState:MaxState
        StrNb = ['"v=', num2str(Syst.Molecule(1).Levelvqn(i)), ',J=', num2str(Syst.Molecule(1).Leveljqn(i)), '"'];
        fprintf(fileID,'%i %s\n', LevelToNew(i), StrNb);
    end


    fprintf(fileID,'*Arcs %i\n',NEdges);
    for i = MinState:MaxState
        for j = MinState:MaxState
          %if i == j
            if abs(Kji(j,i)) >= MinRate
              %fprintf(fileID,'%i,%i,%e,"direct"\n', j, i, Temp(i,j) );
              fprintf(fileID,'%i %i %e\n', LevelToNew(i), LevelToNew(j), abs(TempMat(j,i)) );
            end
          %end
        end
    end
    fclose(fileID);
    clear TempMat;
    timee = toc;
    fprintf('  Inelastic Rates Matrix written in %e s\n', timee)


    FileName1 = strcat(WriteFldr, './LevelsForInelastic.csv');
    fileID1   = fopen(FileName1,'w');
    fprintf(fileID1,'Id,vqn,jqn,EeV,rIn,iLevel\n');
    for i = MinState:MaxState
        fprintf(fileID1,'%i,%i,%i,%e,%e,%i\n', LevelToNew(i), Syst.Molecule(1).Levelvqn(i), Syst.Molecule(1).Leveljqn(i), Syst.Molecule(1).LevelEeV(i), Syst.Molecule(1).LevelrIn(i), i);
    end
    fclose(fileID1);
    
    
    fprintf('====================================================\n\n')

end