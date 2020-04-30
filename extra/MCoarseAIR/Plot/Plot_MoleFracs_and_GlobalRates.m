%% The Function plots the Mole Fractions of the Chemical System's Components 
%%    on Top of The Global Rates (Dissociation and Exchanges)
%
function Plot_MoleFracs_and_GlobalRates(Controls)    
    
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

    
    
    figure(Input.iFig)
    fig = gcf;
    screensize   = get( groot, 'Screensize' );
    %fig.Position = screensize;
    %fig.Color='None';
    left_color  = Syst.CFDComp(Controls.CompStart).Color;
    right_color = Param.KCVec;
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
        
    
    
    yyaxis right

    ProcNames = [];
    loglog(Kin.T(Temp.iT).t, Rates.T(Temp.iT).DissGlobal, 'Color', Param.CMat(1,:), 'linestyle', char(Param.linS(1)), 'LineWidth', Param.LineWidth)
    ProcNames = {'$\bar{K}^D$'};
    hold on
    for iExch = 1:Syst.NProc-2
        loglog(Kin.T(Temp.iT).t, Rates.T(Temp.iT).ExchGlobal(:,iExch), 'Color', Param.CMat(iExch+1,:), 'linestyle', char(Param.linS(iExch+1)), 'LineWidth', Param.LineWidth)
        ProcNames = [ProcNames, strcat('$\bar{K}_{', Syst.Molecule(Syst.ExchToMol(iExch)).Name,'}^E$')];
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

    str_y = ['$\bar{K}$~$[cm^3/s]$'];
    ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    ylab.Interpreter = 'latex';
    %ylim(YLimPlot);
    
    set(0,'DefaultLegendAutoUpdate','off')
    
    
    
    yyaxis left
    
    CompNames = [];
    for iComp = Controls.CompStart:Controls.CompEnd
      semilogx(Kin.T(Temp.iT).t(:), Kin.T(Temp.iT).MolFracs(:,iComp), 'Color', Syst.CFDComp(iComp).Color, 'linestyle', Syst.CFDComp(iComp).LineStyle, 'LineWidth', Param.LineWidth)
      hold on
      CompNames = [CompNames, Syst.CFDComp(iComp).Name];
    end
    hold on
    
    xt = get(gca, 'XTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

%     clab             = legend(CompNames, 'Location', 'Best');
%     clab.Interpreter = 'latex';
%     set(clab,'FontSize', Param.LegendFontSz, 'FontName', Param.LegendFontNm, 'Interpreter', 'latex');

    str_x = ['t [s]'];
    xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    xlab.Interpreter = 'latex';
    %xlim(XLimPlot);

    str_y = ['Mole Fraction'];
    ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    ylab.Interpreter = 'latex';
    %ylim(YLimPlot);
    
    
    
    semilogx([Kin.T(Temp.iT).QSS.tStart, Kin.T(Temp.iT).QSS.tStart], [0, 1], ':k', 'LineWidth',2)
    semilogx([Kin.T(Temp.iT).QSS.tEnd,   Kin.T(Temp.iT).QSS.tEnd],   [0, 1], ':k', 'LineWidth',2)
    
    
    
    pbaspect([1 1 1])

    if Input.SaveFigsFlgInt > 0
        [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
        FolderPath = strcat(Input.Paths.SaveFigsFldr, '/T_', Temp.TNowChar, 'K_', Input.Kin.Proc.OverallFlg, '/');
        [status,msg,msgID] = mkdir(FolderPath);
        if Input.SaveFigsFlgInt == 1
            FileName   = strcat(FolderPath, 'MoleFractionsAndGlobalRates');
            export_fig(FileName, '-pdf')
        elseif Input.SaveFigsFlgInt == 2
            FileName   = strcat(FolderPath, 'MoleFractionsAndGlobalRates.fig');
            savefig(FileName)
        end
        close
    end
    Input.iFig = Input.iFig + 1;

    
end