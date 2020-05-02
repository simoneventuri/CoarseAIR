%% Computing the Energy Depleated from each Internal Mode
%
function Compute_EnergyDepletions(Controls)
    
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


    fprintf('= Compute_EnergyDepletions ============= T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    
    iProj = Controls.ProjTarg(1);
    iTarg = Controls.ProjTarg(2);
    

    for iMol = Controls.MoleculesOI
        fprintf(['Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );

        LevelToBin   = Syst.Molecule(iMol).LevelToBin;
        Levelvqn     = Syst.Molecule(iMol).Levelvqn;
        LevelEeV     = Syst.Molecule(iMol).LevelEeV;
        LevelEeV0    = Syst.Molecule(iMol).LevelEeV0;
        LevelEeVRot  = Syst.Molecule(iMol).LevelEeVRot;
        LevelEeVVib0 = Syst.Molecule(iMol).LevelEeVVib0;
        Nvqn         = Syst.Molecule(iMol).Nvqn;
        NLevels      = Syst.Molecule(iMol).NLevels;
    
        
        KRemoval = ones(NLevels,1) .* 1.d-30;
        for iProc = Controls.RemovalProc
            KRemoval(:) = KRemoval(:) + Rates.T(Temp.iT).Overall(:,iProc);
        end
        
        
        TempMat  = zeros(NLevels, Kin.T(Temp.iT).NSteps);
        jStep    = 1;
        for iStep = 1:Kin.T(Temp.iT).NSteps
            
            if strcmp(Syst.Molecule(iMol).KinMthdIn, 'StS')
                LevelPop(:,1) = Kin.T(Temp.iT).Molecule(iMol).Pop(iStep,:);            
            else
                LevelPop(:,1) = Kin.T(Temp.iT).Molecule(iMol).PopOverg(iStep,LevelToBin(:)) .* Kin.T(Temp.iT).Molecule(iMol).Levelq(:);
            end
            PopTot = sum(LevelPop);
            

            rhoA    = Kin.T(Temp.iT).nd(iStep) * Kin.T(Temp.iT).MolFracs(iStep,iProj);
            rhoM    = Kin.T(Temp.iT).nd(iStep) * Kin.T(Temp.iT).MolFracs(iStep,iTarg);
            PopA    = PopTot / Kin.T(Temp.iT).MolFracs(iStep,iTarg) * Kin.T(Temp.iT).MolFracs(iStep,iProj);
            rhoI(:) = rhoM .* LevelPop(:) ./ PopTot;
            
            TempVec(:,1)  = KRemoval(:) .* PopA .* LevelPop(:) .* ( rhoA.^2 - rhoI(:) );
            Determ        = sum( TempVec );
            
            CDInt(jStep)  = sum( TempVec .* LevelEeV0 )     ./ Determ;
            CDVib(jStep)  = sum( TempVec .* LevelEeVVib0' ) ./ Determ;
            CDRot(jStep)  = sum( TempVec .* LevelEeVRot' )  ./ Determ;

            
            jStep = jStep + 1;
        end

        Kin.T(Temp.iT).Molecule(iMol).CDInt = CDInt ./ abs(LevelEeV(1));
        Kin.T(Temp.iT).Molecule(iMol).CDVib = CDVib ./ abs(LevelEeV(1));
        Kin.T(Temp.iT).Molecule(iMol).CDRot = CDRot ./ abs(LevelEeV(1));
        
        CDIntEq  = Kin.T(Temp.iT).Molecule(iMol).CDInt(end);
        CDRotEq  = Kin.T(Temp.iT).Molecule(iMol).CDRot(end);
        CDVibEq  = Kin.T(Temp.iT).Molecule(iMol).CDVib(end);

        fprintf('At Eq., Int. Energy Depletion Coefficient = %e \n',          CDIntEq );
        fprintf('At Eq., Rot. Energy Depletion Coefficient = %e (%e%%) \n',   CDRotEq, CDRotEq/CDIntEq*100 );
        fprintf('At Eq., Vib. Energy Depletion Coefficient = %e (%e%%) \n\n', CDVibEq, CDVibEq/CDIntEq*100 );
        
        
        CDIntQSS = Kin.T(Temp.iT).Molecule(iMol).CDInt(Kin.T(Temp.iT).QSS.i);
        CDRotQSS = Kin.T(Temp.iT).Molecule(iMol).CDRot(Kin.T(Temp.iT).QSS.i);
        CDVibQSS = Kin.T(Temp.iT).Molecule(iMol).CDVib(Kin.T(Temp.iT).QSS.i);

        fprintf('At QSS, Int. Energy Depletion Coefficient = %e \n',          CDIntQSS );
        fprintf('At QSS, Rot. Energy Depletion Coefficient = %e (%e%%) \n',   CDRotQSS, CDRotQSS/CDIntEq*100 );
        fprintf('At QSS, Vib. Energy Depletion Coefficient = %e (%e%%) \n\n', CDVibQSS, CDVibQSS/CDIntEq*100 );

        
        [status,msg,msgID] = mkdir(Input.Paths.SaveDataFldr);
        FileName          = strcat(Input.Paths.SaveDataFldr, '/EDCoeffs_', Syst.Molecule(iMol).Name, '_', Input.Kin.Proc.OverallFlg, '.csv');
        if exist(FileName, 'file')
            fileID1  = fopen(FileName,'a');
        else
            fileID1  = fopen(FileName,'w');
            HeaderStr = strcat('# T [K], C_Int Eq, C_Rot Eq, C_Vib Eq, C_Int QSS, C_Rot QSS, C_Vib QSS\n');
            fprintf(fileID1,HeaderStr);
        end
        fprintf(fileID1,'%e,%e,%e,%e,%e,%e,%e\n', Temp.TNow, CDIntEq, CDRotEq, CDVibEq, CDIntQSS, CDRotQSS, CDVibQSS );
        fclose(fileID1);
       
        
    end

    fprintf('====================================================\n\n')        

    

end