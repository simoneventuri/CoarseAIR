%% Compute Initial and Final Instants of Quasi-Steady State
%
%  Input Global Var: - Temp.TNowChar
%                    - Syst.HDF5_File
%
function Compute_QSS()

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
    
    global Input Rates Syst Temp Param Kin

    fprintf('= Compute_QSS ========================== T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    
    
    KQSS_Eps = 1.e-9;
    EpsT     = 1.e-4;%5.e-8;
    EpsTT    = 1.e-4;
    InpPerc  = 0.1;

    
    yy       = Rates.T(Temp.iT).DissGlobal+1e-100;
    NSteps   = size(yy,1);
    ExitFlg  = false;
    
    figure(1234)
    loglog(Kin.T(Temp.iT).t, yy, 'k-', 'LineWidth', 2);
    hold on
    
    it = NSteps;
    while ( abs(log10(yy(it)) - log10(yy(end)))     < 1.e-2 )
       it = it - 1;
    end
    itFinal = it;
    tFinal  = Kin.T(Temp.iT).t(itFinal);
    loglog(Kin.T(Temp.iT).t(itFinal), yy(itFinal),'r+', 'MarkerSize', 8);
    if (itFinal < 2)
       ExitFlg = true; 
    end
    
    while ( abs(log10(yy(it)) - log10(yy(itFinal))) < 1.e-1 )
        it = it - 1;
    end
    itEnd   = it;
    tEnd    = Kin.T(Temp.iT).t(itEnd);
    loglog(Kin.T(Temp.iT).t(itEnd), yy(itEnd),'r+', 'MarkerSize', 8);

    while (it >= 1) && ( abs(log10(yy(it)) - log10(yy(it+1)))    > 1.e-4 )
        it = it - 1;
    end
    if (it < 1)
       ExitFlg = true; 
       itQSS   = NSteps;
    else 
       itQSS = it;
    end
    tQSS  = Kin.T(Temp.iT).t(itQSS);
    loglog(Kin.T(Temp.iT).t(itQSS), yy(itQSS),'r+', 'MarkerSize', 8);
    
    if (~ ExitFlg)
        it = itQSS;
        while ( abs(log10(yy(it)) - log10(yy(itQSS)))    < 1.e-1 )
            it = it + 1;
        end
        itEnd = it;
        tEnd  = Kin.T(Temp.iT).t(itEnd);
        loglog(Kin.T(Temp.iT).t(itEnd), yy(itEnd),'r+', 'MarkerSize', 8);
    end
    
    if (~ ExitFlg)
        it = itQSS;
        while ( abs(log10(yy(it)) - log10(yy(itQSS)))   < 1.e-1 ) && (it > 1)
            it = it - 1;
        end
        if it <= 1
            ExitFlg = true;
            itStart = NSteps;
        else
            itStart = it;
        end
        tStart  = Kin.T(Temp.iT).t(itStart);
        loglog(Kin.T(Temp.iT).t(itStart), yy(itStart),'r+', 'MarkerSize', 8);
    end
    
    if (~ ExitFlg)
        while ( abs(log10(yy(it)) - log10(yy(itStart))) < 1.e-1 )
            it = it - 1;
        end
        itIni   = it;
        tIni    = Kin.T(Temp.iT).t(itIni);
        loglog(Kin.T(Temp.iT).t(itIni), yy(itIni),'r+', 'MarkerSize', 8);
    end
    
    if (~ ExitFlg)

        clear fitresult yyy
        yyy(itIni:itFinal) = yy(itIni:itFinal) - yy(itQSS)*(1.0-InpPerc);
        [xData, yData]   = prepareCurveData( Kin.T(Temp.iT).t(itIni:itFinal), yyy(itIni:itFinal) );
        ft = 'splineinterp';
        [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
        x0 = Kin.T(Temp.iT).t(itStart);
        Kin.T(Temp.iT).QSS.tStart = fzero(fitresult, x0);
        it = 1;
        while Kin.T(Temp.iT).t(it) < Kin.T(Temp.iT).QSS.tStart
            it = it + 1;
        end
        Kin.T(Temp.iT).QSS.iStart = it - 1; 
        if (Kin.T(Temp.iT).QSS.iStart < 1) || (Kin.T(Temp.iT).QSS.iStart > NSteps)
            ExitFlg = true;
        end

        clear fitresult yyy xData yData
        if (yy(itQSS) > yy(end))
            yyy(itStart:itEnd) = yy(itStart:itEnd) - yy(itQSS)*(1.0-InpPerc);
        else
            yyy(itStart:itEnd) = yy(itStart:itEnd) - yy(itQSS)*(1.0+InpPerc);
        end
        [xData, yData]   = prepareCurveData( Kin.T(Temp.iT).t(itStart:itEnd), yyy(itStart:itEnd) );
        ft = 'splineinterp';
        [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
        x0 = Kin.T(Temp.iT).t(itEnd);
        Kin.T(Temp.iT).QSS.tEnd = fzero(fitresult, x0);
        it = 1;
        while (Kin.T(Temp.iT).t(it) < Kin.T(Temp.iT).QSS.tEnd) && (it < NSteps)
            it = it + 1;
        end       
        Kin.T(Temp.iT).QSS.iEnd   = it;

        tQSS = 10.0^((log10(Kin.T(Temp.iT).QSS.tStart) + log10(Kin.T(Temp.iT).QSS.tEnd))/2.0);
        itQSS = 1;
        while Kin.T(Temp.iT).t(itQSS) < tQSS
            itQSS = itQSS + 1;
        end
        Kin.T(Temp.iT).QSS.i = itQSS-1;
        Kin.T(Temp.iT).QSS.t = Kin.T(Temp.iT).t(Kin.T(Temp.iT).QSS.i);
    end 
    
    
    if (ExitFlg)
        Kin.T(Temp.iT).QSS.iStart = NSteps;
        Kin.T(Temp.iT).QSS.i      = NSteps;
        Kin.T(Temp.iT).QSS.iEnd   = NSteps;

        Kin.T(Temp.iT).QSS.tStart = Kin.T(Temp.iT).t(NSteps);
        Kin.T(Temp.iT).QSS.t      = Kin.T(Temp.iT).t(NSteps);
        Kin.T(Temp.iT).QSS.tEnd   = Kin.T(Temp.iT).t(NSteps);
    end

    loglog(Kin.T(Temp.iT).t(Kin.T(Temp.iT).QSS.iStart), yy(Kin.T(Temp.iT).QSS.iStart), 'ko', 'MarkerSize', 8);
    loglog(Kin.T(Temp.iT).t(Kin.T(Temp.iT).QSS.iEnd), yy(Kin.T(Temp.iT).QSS.iEnd), 'ko', 'MarkerSize', 8);
    loglog(Kin.T(Temp.iT).t(Kin.T(Temp.iT).QSS.i), yy(Kin.T(Temp.iT).QSS.i), 'ko', 'MarkerSize', 8);       fprintf('QSS Start Time       = %e s\n',  Kin.T(Temp.iT).QSS.tStart );
    
    fprintf('QSS Start Time Step  = %i \n',   Kin.T(Temp.iT).QSS.iStart );
    fprintf('QSS       Time       = %e s\n',  Kin.T(Temp.iT).QSS.t );
    fprintf('QSS       Time Step  = %i \n',   Kin.T(Temp.iT).QSS.i );
    fprintf('QSS End   Time       = %e s\n',  Kin.T(Temp.iT).QSS.tEnd );
    fprintf('QSS End   Time Step  = %i \n\n', Kin.T(Temp.iT).QSS.iEnd );
    
 
    
    KDissEq   = Rates.T(Temp.iT).DissGlobal(end);
    KDissQSS  = Rates.T(Temp.iT).DissGlobal(itQSS);
    Kin.T(Temp.iT).QSS.Diss      = KDissQSS;

    KExch1Eq  = Rates.T(Temp.iT).ExchGlobal(end,1);
    KExch1QSS = Rates.T(Temp.iT).ExchGlobal(itQSS,1);
    Kin.T(Temp.iT).QSS.Exch1     = KExch1QSS;
    
    if Syst.NProc == 4
        KExch2Eq  = Rates.T(Temp.iT).ExchGlobal(end,2);
        KExch2QSS = Rates.T(Temp.iT).ExchGlobal(itQSS,2);
        Kin.T(Temp.iT).QSS.Exch2 = KExch2QSS;
    end

    
    fprintf('QSS KDiss            = %e [cm^3/s]\n', Kin.T(Temp.iT).QSS.Diss );
    for iExch = 1:Syst.NProc-2
        fprintf('QSS KExch, Exch Nb %i = %e [cm^3/s] \n', iExch,  Rates.T(Temp.iT).ExchGlobal(itQSS,iExch) );
    end
    fprintf('\n')
    
    
    %% Writing Dissociation and Exchange Values at Equilibrium and QSS 
    %
    [status,msg,msgID] = mkdir(Input.Paths.SaveDataFldr);
    FileName           = strcat(Input.Paths.SaveDataFldr, '/KGlobal_', Input.Kin.Proc.OverallFlg, '.csv');
    if exist(FileName, 'file')
        fileID1  = fopen(FileName,'a');
    else
        fileID1  = fopen(FileName,'w');
        if Syst.NProc == 3
            HeaderStr = strcat('# T [K], K^D Eq, K_{', Syst.Molecule(Syst.ExchToMol(1)).Name, '}^E Eq, K^D QSS, K_{', Syst.Molecule(Syst.ExchToMol(1)).Name,'}^E QSS \n');
        else
            HeaderStr = strcat('# T [K], K^D Eq, K_{', Syst.Molecule(Syst.ExchToMol(1)).Name, '}^E Eq, K_{', Syst.Molecule(Syst.ExchToMol(2)).Name, '}^E Eq, K^D QSS, K_{', Syst.Molecule(Syst.ExchToMol(1)).Name,'}^E QSS, K_{', Syst.Molecule(Syst.ExchToMol(2)).Name, '}^E QSS \n');
        end
        fprintf(fileID1,HeaderStr);
    end
    if Syst.NProc == 3
        fprintf(fileID1,'%e,%e,%e,%e,%e\n',       Temp.TNow, KDissEq, KExch1Eq, KDissQSS, KExch1QSS );
    else
        fprintf(fileID1,'%e,%e,%e,%e,%e,%e,%e\n', Temp.TNow, KDissEq, KExch1Eq, KExch2Eq, KDissQSS, KExch1QSS, KExch2QSS );
    end
    fclose(fileID1);

    
    %% Writing Percentage of Molecule Depletion Happening During QSS
    %
    for iComp = 1:Syst.NComp
        if Syst.CFDComp(iComp).ToMol > 0
            iMol = Syst.CFDComp(iComp).ToMol;
            fprintf(['Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );

            MolFrac    = Kin.T(Temp.iT).MolFracs(:,iComp);
            MolFracIn  = MolFrac(1);
            MolFracS   = MolFrac(Kin.T(Temp.iT).QSS.iStart);
            MolFracE   = MolFrac(Kin.T(Temp.iT).QSS.iEnd);
            MolFracFin = MolFrac(end);
            
            
            Kin.T(Temp.iT).QSS.Molecule(iMol).PercPre  = abs(MolFracIn - MolFracS)   / abs(MolFracIn - MolFracFin) .* 100.0;
            Kin.T(Temp.iT).QSS.Molecule(iMol).PercAt   = abs(MolFracE  - MolFracS)   / abs(MolFracIn - MolFracFin) .* 100.0;
            Kin.T(Temp.iT).QSS.Molecule(iMol).PercPost = abs(MolFracE  - MolFracFin) / abs(MolFracIn - MolFracFin) .* 100.0;
            fprintf('%e%% of Dissociation happens Pre-QSS\n',    Kin.T(Temp.iT).QSS.Molecule(iMol).PercPre );
            fprintf('%e%% of Dissociation happens During QSS\n', Kin.T(Temp.iT).QSS.Molecule(iMol).PercAt );
            fprintf('%e%% of Dissociation happens Post-QSS\n',   Kin.T(Temp.iT).QSS.Molecule(iMol).PercPost );


            [status,msg,msgID] = mkdir(Input.Paths.SaveDataFldr);
            FileName           = strcat(Input.Paths.SaveDataFldr, '/DuringQSS_', Syst.Molecule(iMol).Name ,'_', Input.Kin.Proc.OverallFlg, '.csv');
            if exist(FileName, 'file')
                fileID1  = fopen(FileName,'a');
            else
                fileID1  = fopen(FileName,'w');
                HeaderStr = strcat('# T [K], % Pre-QSS, % During QSS, % Post QSS \n');
                fprintf(fileID1,HeaderStr);
            end
            fprintf(fileID1,'%e,%e,%e,%e\n', Temp.TNow, Kin.T(Temp.iT).QSS.Molecule(iMol).PercPre, Kin.T(Temp.iT).QSS.Molecule(iMol).PercAt, Kin.T(Temp.iT).QSS.Molecule(iMol).PercPost );
            fclose(fileID1);
        
        end
    end


    fprintf('====================================================\n\n')        
    
end