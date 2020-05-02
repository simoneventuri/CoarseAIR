%% Computing the Energy Stored in each of the Internal Modes
%
function Compute_Energies(Controls)    

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

    
    fprintf('= Compute_Energies ===================== T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    

    Pressure = Kin.T(Temp.iT).MolFracs(end, Syst.ColPartToComp) * Kin.T(Temp.iT).P(end) / 101325.0;
    
    
    for iMol = Controls.MoleculesOI
        fprintf(['Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );
        
        LevelToBin  = Syst.Molecule(iMol).LevelToBin;
        Levelvqn    = Syst.Molecule(iMol).Levelvqn;
        LevelEeV    = Syst.Molecule(iMol).LevelEeV;
        LevelEeV0   = Syst.Molecule(iMol).LevelEeV0;
        LevelEeVRot = Syst.Molecule(iMol).LevelEeVRot;
        vEeVVib0    = Syst.Molecule(iMol).vEeVVib0;
        Nvqn        = Syst.Molecule(iMol).Nvqn;
        NLevels     = Syst.Molecule(iMol).NLevels;
        
        
        eInt = zeros(Kin.T(Temp.iT).NSteps,1);
        eRot = zeros(Kin.T(Temp.iT).NSteps,1);
        eVib = zeros(Kin.T(Temp.iT).NSteps,1);        
        vPop = zeros(Nvqn,1);
        for iStep = 1:Kin.T(Temp.iT).NSteps

            
            if strcmp(Syst.Molecule(iMol).KinMthdIn, 'StS')
                LevelPop(:,1) = Kin.T(Temp.iT).Molecule(iMol).Pop(iStep,:);            
            else
                LevelPop(:,1) = Kin.T(Temp.iT).Molecule(iMol).PopOverg(iStep,LevelToBin(:)) .* Kin.T(Temp.iT).Molecule(iMol).Levelq(:);
            end
            PopTot = sum(LevelPop);
            
            vPop = vPop.*0.0;
            for iLevels = 1:NLevels
                vPop(Levelvqn(iLevels)+1) = vPop(Levelvqn(iLevels)+1) + LevelPop(iLevels);
            end      

            
            eInt(iStep) = sum( LevelEeV0    .* LevelPop ) ./ PopTot;
            eRot(iStep) = sum( LevelEeVRot' .* LevelPop ) ./ PopTot;
            eVib(iStep) = sum( vEeVVib0'    .* vPop )     ./ PopTot;
            
            
        end
        
        
        Kin.T(Temp.iT).Molecule(iMol).eInt = eInt;
        Kin.T(Temp.iT).Molecule(iMol).eRot = eRot;
        Kin.T(Temp.iT).Molecule(iMol).eVib = eVib;
        
        
        eIntLim = (eInt(end) - eInt(1)) * 0.632 + eInt(1);
        iInt=1;
        while eInt(iInt) < eIntLim
          iInt = iInt+1;
        end
        tauInt = (Kin.T(Temp.iT).t(iInt) + Kin.T(Temp.iT).t(iInt-1)) / 2.d0;
        [xData, yData] = prepareCurveData( Kin.T(Temp.iT).t, eInt - eIntLim );
        ft = 'splineinterp';
        [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
        Kin.T(Temp.iT).Molecule(iMol).tauInt = fzero(fitresult, tauInt);
        
        eRotLim = (eRot(end) - eRot(1)) * 0.632 + eRot(1);
        iRot=1;
        while eRot(iRot) < eRotLim
          iRot = iRot+1;
        end
        tauRot = (Kin.T(Temp.iT).t(iRot) + Kin.T(Temp.iT).t(iRot-1)) / 2.d0;
        [xData, yData] = prepareCurveData( Kin.T(Temp.iT).t, eRot - eRotLim );
        ft = 'splineinterp';
        [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
        Kin.T(Temp.iT).Molecule(iMol).tauRot = fzero(fitresult, tauRot);
        
        eVibLim = (eVib(end) - eVib(1)) * 0.632 + eVib(1);
        iVib=1;
        while eVib(iVib) < eVibLim
          iVib = iVib+1;
        end
        tauVib = (Kin.T(Temp.iT).t(iVib) + Kin.T(Temp.iT).t(iVib-1)) / 2.d0;
        [xData, yData] = prepareCurveData( Kin.T(Temp.iT).t, eVib - eVibLim );
        ft = 'splineinterp';
        [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
        Kin.T(Temp.iT).Molecule(iMol).tauVib = fzero(fitresult, tauVib);
        
        
        Kin.T(Temp.iT).Molecule(iMol).tauIntP = Kin.T(Temp.iT).Molecule(iMol).tauInt * Pressure;
        Kin.T(Temp.iT).Molecule(iMol).tauRotP = Kin.T(Temp.iT).Molecule(iMol).tauRot * Pressure;
        Kin.T(Temp.iT).Molecule(iMol).tauVibP = Kin.T(Temp.iT).Molecule(iMol).tauVib * Pressure;
        
        fprintf('P*tau_Int = %e [Pa*s]\n',   Kin.T(Temp.iT).Molecule(iMol).tauIntP );
        fprintf('P*tau_Rot = %e [Pa*s]\n',   Kin.T(Temp.iT).Molecule(iMol).tauRotP );
        fprintf('P*tau_Vib = %e [Pa*s]\n\n', Kin.T(Temp.iT).Molecule(iMol).tauVibP );


       [status,msg,msgID] = mkdir(Input.Paths.SaveDataFldr);
        FileName          = strcat(Input.Paths.SaveDataFldr, '/Taus_', Syst.Molecule(iMol).Name, '_', Input.Kin.Proc.OverallFlg, '.csv');
        if exist(FileName, 'file')
            fileID1  = fopen(FileName,'a');
        else
            fileID1  = fopen(FileName,'w');
            HeaderStr = strcat('# T [K], P [Pa], tau_Int, tau_Rot, tau_Vib\n');
            fprintf(fileID1,HeaderStr);
        end
        fprintf(fileID1,'%e,%e,%e,%e,%e\n', Temp.TNow, Pressure, Kin.T(Temp.iT).Molecule(iMol).tauIntP, Kin.T(Temp.iT).Molecule(iMol).tauRotP, Kin.T(Temp.iT).Molecule(iMol).tauVibP  );
        fclose(fileID1);
        
        
    end

    fprintf('====================================================\n\n')        

    
end