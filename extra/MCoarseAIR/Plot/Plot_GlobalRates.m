%% The Function plots the Global Rates (Dissociation and Exchanges)
%
function Plot_GlobalRates(Controls)    
    
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

    global Input Kin Param Syst OtherSyst Temp Rates

    fprintf('= Plot_GlobalRates ===================== T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    
    
    figure(Input.iFig)
    fig = gcf;
    screensize   = get( groot, 'Screensize' );
    %fig.Position = screensize;
    %fig.Color='None';
    
    iC        = 1;
    ProcNames = [];
    for iMol = Controls.MoleculesOI
        if (Syst.NAtoms == 4)
            NExch = 0;
        else
            if (iMol == 1)
                NExch = Syst.NProc-2;
            else
                NExch = OtherSyst(iMol-1).Syst.NProc-2;
            end
        end    
        loglog(Kin.T(Temp.iT).t, Rates.T(Temp.iT).Molecule(iMol).DissGlobal(:,1), 'Color', Param.CMat(iC,:), 'linestyle', char(Param.linS(1)), 'LineWidth', Param.LineWidth)
        iC        = iC + 1;
        ProcNames = [ProcNames, {[Syst.Molecule(iMol).Name, ', $k_{Global}^D$']}];
        hold on
        for iExch = 1:NExch
            if (iMol == 1)
                TempStr = [Syst.Molecule(iMol).Name, ', $k_{Global}^{E_{', Syst.Molecule(Syst.ExchToMol(iExch)).Name,'}}$'];
            else
                TempStr = [Syst.Molecule(iMol).Name, ', $k_{Global}^{E_{', OtherSyst(iMol-1).Syst.Molecule(OtherSyst(iMol-1).Syst.ExchToMol(iExch)).Name,'}}$'];
            end
            loglog(Kin.T(Temp.iT).t, Rates.T(Temp.iT).Molecule(iMol).ExchGlobal(:,iExch), 'Color', Param.CMat(iC,:), 'linestyle', char(Param.linS(1)), 'LineWidth', Param.LineWidth)
            iC        = iC + 1;
            ProcNames = [ProcNames, TempStr];
        end
    end
    hold on

    xt = get(gca, 'XTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

    clab             = legend(ProcNames, 'Location', 'Best');
    clab.Interpreter = 'latex';
    set(clab,'FontSize', Param.LegendFontSz, 'FontName', Param.LegendFontNm, 'Interpreter', 'latex');

    str_x = ['t [s]'];
    xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    xlab.Interpreter = 'latex';
    %xlim(XLimPlot);

    str_y = ['$k_{Global}$~$[cm^3/s]$'];
    ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    ylab.Interpreter = 'latex';
    %ylim(YLimPlot);


  
    pbaspect([1 1 1])
    
    if Input.SaveFigsFlgInt > 0
        [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr);
        FolderPath = strcat(Input.Paths.SaveFigsFldr, '/T_', Temp.TNowChar, 'K_', Input.Kin.Proc.OverallFlg, '/');
        [status,msg,msgID] = mkdir(FolderPath);
        if Input.SaveFigsFlgInt == 1
            FileName   = strcat(FolderPath, 'GlobalRates');
            export_fig(FileName, '-pdf');
        elseif Input.SaveFigsFlgInt == 2
            FileName   = strcat(FolderPath, 'GlobalRates.fig');
            savefig(FileName);
        end
        close
    end
    Input.iFig = Input.iFig + 1;

    
    fprintf('====================================================\n\n')

end