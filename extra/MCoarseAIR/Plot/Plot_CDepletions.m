%% The Function reads and plots the Vectors of Dissociation and Exchange Overall Rate Coefficients (at Equilibrium and QSS)
%
function Plot_CDepletions()

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

    global Input Param

    fprintf('= Plot_CDepletions =================================\n')
    fprintf('====================================================\n')

    
    opts = delimitedTextImportOptions("NumVariables", 7);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["TVec", "CDIEq", "CDREq", "CDVEq", "CDIQSS", "CDRQSS", "CDVQSS"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    tbl = readtable(Input.Paths.CDepletions, opts);
    TVec   = tbl.TVec;
    CDIEq  = tbl.CDIEq;
    CDREq  = tbl.CDREq;
    CDVEq  = tbl.CDVEq;
    CDIQSS = tbl.CDIQSS;
    CDRQSS = tbl.CDRQSS;
    CDVQSS = tbl.CDVQSS;
    clear opts tbl


    
    figure(Input.iFig)
    fig = gcf;
    screensize   = get( groot, 'Screensize' );
    %fig.Position = screensize;
    %fig.Color='None';

    h1 = plot(TVec, CDIEq, 'Color', Param.KCVec, 'linestyle', char(Param.linS(1)), 'LineWidth', Param.LineWidth);
    hold on
    h2 = scatter(TVec,  CDIEq, 200, '.', 'MarkerEdgeColor', Param.KCVec, 'MarkerFaceColor', Param.KCVec, 'LineWidth', 1.5);
    h3 = plot(TVec, CDVEq, 'Color', Param.RCVec, 'linestyle', char(Param.linS(1)), 'LineWidth', Param.LineWidth);
    h4 = scatter(TVec,  CDVEq, 60, '+', 'MarkerEdgeColor', Param.RCVec, 'MarkerFaceColor', Param.RCVec, 'LineWidth', 1.5);
    h5 = plot(TVec, CDREq, 'Color', Param.BCVec, 'linestyle', char(Param.linS(1)), 'LineWidth', Param.LineWidth);
    h6 = scatter(TVec,  CDREq, 60, 'x', 'MarkerEdgeColor', Param.BCVec, 'MarkerFaceColor', Param.BCVec, 'LineWidth', 1.5);
    PlotNames = [{'$C_{Eq}^{DI}$'}, {strcat('$C_{Eq}^{DV}$')}, {'$C_{Eq}^{DR}$'}];

%     xt = get(gca, 'XTick');
%     set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
%     yt = get(gca, 'YTick');
%     set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
% 
%     clab             = legend(PlotNames, 'Location', 'Best');
%     clab.Interpreter = 'latex';
%     set(clab,'FontSize', Param.LegendFontSz, 'FontName', Param.LegendFontNm, 'Interpreter', 'latex');
% 
%     str_x = ['T [K]'];
%     xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
%     xlab.Interpreter = 'latex';
%     %xlim(XLimPlot);
% 
%     str_y = ['CE Coefficient'];
%     ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
%     ylab.Interpreter = 'latex';
%     %ylim(YLimPlot);

%     if Input.SaveFigsFlgInt > 0
%         [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
%         FolderPath = strcat(Input.Paths.SaveFigsFldr, '/');
%         [status,msg,msgID] = mkdir(FolderPath);
%         if Input.SaveFigsFlgInt == 1
%             FileName   = strcat(FolderPath, 'CDepletion_Eq');
%             export_fig(FileName, '-pdf')
%         elseif Input.SaveFigsFlgInt == 2
%             FileName   = strcat(FolderPath, 'CDepletion_Eq.fig');
%             savefig(FileName)
%         end
%         %close
%     end
%     Input.iFig = Input.iFig + 1;
% 
% 
%     
%     figure(Input.iFig)
%     fig = gcf;
%     screensize   = get( groot, 'Screensize' );
%     %fig.Position = screensize;
%     %fig.Color='None';

    h11 = plot(TVec, CDIQSS, 'Color', Param.KCVec, 'linestyle', char(Param.linS(2)), 'LineWidth', Param.LineWidth);
    hold on
    h22 = scatter(TVec,  CDIQSS, 200, '.', 'MarkerEdgeColor', Param.KCVec, 'MarkerFaceColor', Param.KCVec, 'LineWidth', 1.5);
    h33 = plot(TVec, CDVQSS, 'Color', Param.RCVec, 'linestyle', char(Param.linS(2)), 'LineWidth', Param.LineWidth);
    h44 = scatter(TVec,  CDVQSS, 60, '+', 'MarkerEdgeColor', Param.RCVec, 'MarkerFaceColor', Param.RCVec, 'LineWidth', 1.5);
    h55 = plot(TVec, CDRQSS, 'Color', Param.BCVec, 'linestyle', char(Param.linS(2)), 'LineWidth', Param.LineWidth);
    h66 = scatter(TVec,  CDRQSS, 60, 'x', 'MarkerEdgeColor', Param.BCVec, 'MarkerFaceColor', Param.BCVec, 'LineWidth', 1.5);
    PlotNames = [PlotNames, {'$C_{QSS}^{DI}$'}, {strcat('$C_{QSS}^{DV}$')}, {'$C_{QSS}^{DR}$'}];

    xt = get(gca, 'XTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

    clab             = legend([h1,h3,h5,h11,h33,h55], PlotNames, 'Location', 'Best');
    clab.Interpreter = 'latex';
    set(clab,'FontSize', Param.LegendFontSz, 'FontName', Param.LegendFontNm, 'Interpreter', 'latex');

    str_x = ['T [K]'];
    xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    xlab.Interpreter = 'latex';
    %xlim(XLimPlot);

    str_y = ['CE Coefficients'];
    ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    ylab.Interpreter = 'latex';
    %ylim(YLimPlot);

%     if Input.SaveFigsFlgInt > 0
%         [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
%         FolderPath = strcat(Input.Paths.SaveFigsFldr, '/');
%         [status,msg,msgID] = mkdir(FolderPath);
%         if Input.SaveFigsFlgInt == 1
%             FileName   = strcat(FolderPath, 'CDepletion_QSS');
%             export_fig(FileName, '-pdf')
%         elseif Input.SaveFigsFlgInt == 2
%             FileName   = strcat(FolderPath, 'CDepletion_QSS.fig');
%             savefig(FileName)
%         end
%         %close
%     end
%     Input.iFig = Input.iFig + 1;
%     
    if Input.SaveFigsFlgInt > 0
        [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
        FolderPath = strcat(Input.Paths.SaveFigsFldr, '/');
        [status,msg,msgID] = mkdir(FolderPath);
        if Input.SaveFigsFlgInt == 1
            FileName   = strcat(FolderPath, 'CDepletion');
            export_fig(FileName, '-pdf')
        elseif Input.SaveFigsFlgInt == 2
            FileName   = strcat(FolderPath, 'CDepletion.fig');
            savefig(FileName)
        end
        %close
    end
    Input.iFig = Input.iFig + 1;


    fprintf('====================================================\n\n')
    
end