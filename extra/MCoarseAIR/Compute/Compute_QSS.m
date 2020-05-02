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
    EpsT     = 5.e-4;%5.e-8;
    EpsTT    = 1.e-3;
    InpPerc  = 0.1;

    
    yy       = Rates.T(Temp.iT).DissGlobal;


    iQSS = 1;
    while (yy(iQSS) <= yy(1) + 0.1*(yy(end) - yy(1)))
        iQSS = iQSS + 1;
    end
    while ( abs(log10(yy(iQSS)) - log10(yy(iQSS+1))) / abs(log10(yy(iQSS))) > EpsT ) && ( iQSS < length(yy)-3)
        iQSS = iQSS + 1;
    end
    iQSS_Start = iQSS - 1;
    while ( abs(log10(yy(iQSS)) - log10(yy(iQSS+1))) / abs(log10(yy(iQSS))) <= EpsTT ) && ( iQSS < length(yy)-3)
        iQSS = iQSS + 1;
    end
    iQSS_End = iQSS - 1;
    
    iQSS                 = floor( (iQSS_Start + iQSS_End) / 2.0 );
    Kin.T(Temp.iT).QSS.i = iQSS;
    Kin.T(Temp.iT).QSS.t = Kin.T(Temp.iT).t(Kin.T(Temp.iT).QSS.i);
    
    iStart = 1;
    while (yy(iStart) < yy(iQSS)/6)
        iStart = iStart + 1;
    end
    iEnd = size(yy,1);%iQSS;
    %   while (abs(yy(iEnd) - yy(iQSS))/yy(iQSS) < (1.0+InpPerc+0.05))
    %     iEnd = iEnd + 1;
    %   end

    
    clear fitresult yyy
    yyy(iStart:iEnd)  = yy(iStart:iEnd) - yy(iQSS)*(1.0-InpPerc);
    [xData, yData] = prepareCurveData( Kin.T(Temp.iT).t(iStart:iEnd), yyy(iStart:iEnd) );
    ft = 'splineinterp';
    [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
    x0 = Kin.T(Temp.iT).t(iQSS_Start);
    Kin.T(Temp.iT).QSS.tStart = fzero(fitresult, x0);

    it = 1;
    while Kin.T(Temp.iT).t(it) < Kin.T(Temp.iT).QSS.tStart
        it = it + 1;
    end
    Kin.T(Temp.iT).QSS.iStart = it - 1;    
    
%     figure(1234)
%     h1 = plot( fitresult, xData, yData );
%     hold on
%     h2 = semilogx(x0,                        fitresult(x0),                        'ro', 'MarkerSize', 10 );
%     h3 = semilogx(Kin.T(Temp.iT).QSS.tStart, fitresult(Kin.T(Temp.iT).QSS.tStart), 'ko', 'MarkerSize', 10 );
%     set(gca, 'XScale', 'log')

    
    clear fitresult yyy xData yData
    yyy(iStart:iEnd)  = yy(iStart:iEnd) - yy(iQSS)*(1.0+InpPerc);
    [xData, yData] = prepareCurveData( Kin.T(Temp.iT).t(iStart:iEnd), yyy(iStart:iEnd) );
    ft = 'splineinterp';
    [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
    x0 = Kin.T(Temp.iT).t(iQSS)*5;
    Kin.T(Temp.iT).QSS.tEnd   = fzero(fitresult, x0);

    it = 1;
    while Kin.T(Temp.iT).t(it) < Kin.T(Temp.iT).QSS.tEnd
        it = it + 1;
    end       
    Kin.T(Temp.iT).QSS.iEnd   = it;

%     figure(1234)
%     h1 = plot( fitresult, xData, yData );
%     hold on
%     h2 = semilogx(x0,                        fitresult(x0),                    'ro', 'MarkerSize', 10 );
%     h3 = semilogx(Kin.T(Temp.iT).QSS.tEnd, fitresult(Kin.T(Temp.iT).QSS.tEnd), 'ko', 'MarkerSize', 10 );
%     set(gca, 'XScale', 'log')
    
    fprintf('QSS Start Time       = %e s\n', Kin.T(Temp.iT).QSS.tStart );
    fprintf('QSS Start Time Step  = %i \n',  Kin.T(Temp.iT).QSS.iStart );
    fprintf('QSS       Time       = %e s\n', Kin.T(Temp.iT).QSS.t );
    fprintf('QSS       Time Step  = %i \n',  Kin.T(Temp.iT).QSS.i );
    fprintf('QSS End   Time       = %e s\n', Kin.T(Temp.iT).QSS.tEnd );
    fprintf('QSS End   Time Step  = %i \n\n',  Kin.T(Temp.iT).QSS.iEnd );
    
    
    
    KDissEq   = Rates.T(Temp.iT).DissGlobal(end);
    KDissQSS  = Rates.T(Temp.iT).DissGlobal(iQSS);
    Kin.T(Temp.iT).QSS.Diss      = KDissQSS;

    KExch1Eq  = Rates.T(Temp.iT).ExchGlobal(end,1);
    KExch1QSS = Rates.T(Temp.iT).ExchGlobal(iQSS,1);
    Kin.T(Temp.iT).QSS.Exch1     = KExch1QSS;
    
    if Syst.NProc == 4
        KExch2Eq  = Rates.T(Temp.iT).ExchGlobal(end,2);
        KExch2QSS = Rates.T(Temp.iT).ExchGlobal(iQSS,2);
        Kin.T(Temp.iT).QSS.Exch2 = KExch2QSS;
    end

    
    fprintf('QSS KDiss            = %e [cm^3/s]\n', Kin.T(Temp.iT).QSS.Diss );
    for iExch = 1:Syst.NProc-2
        fprintf('QSS KExch, Exch Nb %i = %e [cm^3/s] \n', iExch,  Rates.T(Temp.iT).ExchGlobal(iQSS,iExch) );
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