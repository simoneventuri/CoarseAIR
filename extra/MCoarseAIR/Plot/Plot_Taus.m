%% The Function reads and plots the Vectors of Vibrational and Rotational Relaxation Times as Function of Temperatures
%
function Plot_Taus()

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

    fprintf('= Plot_Taus ========================================\n')
    fprintf('====================================================\n')


    if (~ strcmp(Input.Paths.TausWOExch, ''))
        opts = delimitedTextImportOptions("NumVariables", 5);
        opts.DataLines = [2, Inf];
        opts.Delimiter = ",";
        opts.VariableNames = ["T", "P", "tau_Int", "tau_Rot", "tau_Vib"];
        opts.VariableTypes = ["double", "double", "double", "double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        tbl = readtable(Input.Paths.TausWOExch, opts);
        T_NoExch       = tbl.T;
        P_NoExch       = tbl.P;
        tau_Int_NoExch = tbl.tau_Int;
        tau_Rot_NoExch = tbl.tau_Rot;
        tau_Vib_NoExch = tbl.tau_Vib;
        clear opts tbl
    end

    if (~ strcmp(Input.Paths.TausWExch, ''))
        opts = delimitedTextImportOptions("NumVariables", 5);
        opts.DataLines = [2, Inf];
        opts.Delimiter = ",";
        opts.VariableNames = ["T", "P", "tau_Int", "tau_Rot", "tau_Vib"];
        opts.VariableTypes = ["double", "double", "double", "double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        tbl = readtable(Input.Paths.TausWExch, opts);
        T       = tbl.T;
        P       = tbl.P;
        tau_Int = tbl.tau_Int;
        tau_Rot = tbl.tau_Rot;
        tau_Vib = tbl.tau_Vib;
        clear opts tbl
    end



    figure(Input.iFig)
    fig = gcf;
    screensize   = get( groot, 'Screensize' );
    %fig.Position = screensize;
    %fig.Color='None';

    PlotNames = [];
    if (~ strcmp(Input.Paths.TausWOExch, ''))
        h1 = semilogy(T_NoExch.^(-1/3), tau_Vib_NoExch, 'Color', Param.CMat(1,:), 'linestyle', char(Param.linS(1)), 'LineWidth', Param.LineWidth)
        PlotNames = [PlotNames, {'w/o Exch.'}];
        hold on
    end
    if (~ strcmp(Input.Paths.TausWExch, ''))
        h2 = semilogy(T.^(-1/3),        tau_Vib,        'Color', Param.CMat(2,:), 'linestyle', char(Param.linS(2)), 'LineWidth', Param.LineWidth)
        PlotNames = [PlotNames, {'w/ Exch.'}];
        hold on
    end

    xt = get(gca, 'XTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

    clab             = legend(PlotNames, 'Location', 'Best');
    clab.Interpreter = 'latex';
    set(clab,'FontSize', Param.LegendFontSz, 'FontName', Param.LegendFontNm, 'Interpreter', 'latex');

    str_x = ['T${}^{-1/3}$ [K${}^{-1/3}]$'];
    xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    xlab.Interpreter = 'latex';
    %xlim(XLimPlot);

    str_y = ['$p\tau_v$ [atm$\cdot$s]'];
    ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    ylab.Interpreter = 'latex';
    %ylim(YLimPlot);

    if Input.SaveFigsFlgInt > 0
        [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
        FolderPath = strcat(Input.Paths.SaveFigsFldr, '/');
        [status,msg,msgID] = mkdir(FolderPath);
        if Input.SaveFigsFlgInt == 1
            FileName   = strcat(FolderPath, 'TauVibS');
            export_fig(FileName, '-pdf')
        elseif Input.SaveFigsFlgInt == 2
            FileName   = strcat(FolderPath, 'TauVibS.fig');
            savefig(FileName)
        end
        %close
    end
    Input.iFig = Input.iFig + 1;



    figure(Input.iFig)
    fig = gcf;
    screensize   = get( groot, 'Screensize' );
    %fig.Position = screensize;
    %fig.Color='None';

    PlotNames = [];
    if (~ strcmp(Input.Paths.TausWOExch, ''))
        h1 = semilogy(T_NoExch.^(-1/3), tau_Rot_NoExch, 'Color', Param.CMat(1,:), 'linestyle', char(Param.linS(1)), 'LineWidth', Param.LineWidth)
        PlotNames = [PlotNames, {'w/o Exch.'}];
        hold on
    end
    if (~ strcmp(Input.Paths.TausWExch, ''))
        h2 = semilogy(T.^(-1/3),        tau_Rot,        'Color', Param.CMat(2,:), 'linestyle', char(Param.linS(2)), 'LineWidth', Param.LineWidth)
        PlotNames = [PlotNames, {'w/ Exch.'}];
        hold on
    end

    xt = get(gca, 'XTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

    clab             = legend(PlotNames, 'Location', 'Best');
    clab.Interpreter = 'latex';
    set(clab,'FontSize', Param.LegendFontSz, 'FontName', Param.LegendFontNm, 'Interpreter', 'latex');

    str_x = ['T${}^{-1/3}$ [K${}^{-1/3}]$'];
    xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    xlab.Interpreter = 'latex';
    %xlim(XLimPlot);

    str_y = ['$p\tau_r$ [atm$\cdot$s]'];
    ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    ylab.Interpreter = 'latex';
    %ylim(YLimPlot);

    if Input.SaveFigsFlgInt > 0
        [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
        FolderPath = strcat(Input.Paths.SaveFigsFldr, '/');
        [status,msg,msgID] = mkdir(FolderPath);
        if Input.SaveFigsFlgInt == 1
            FileName   = strcat(FolderPath, 'TauRotS');
            export_fig(FileName, '-pdf')
        elseif Input.SaveFigsFlgInt == 2
            FileName   = strcat(FolderPath, 'TauRotS.fig');
            savefig(FileName)
        end
        %close
    end
    Input.iFig = Input.iFig + 1;


    fprintf('====================================================\n\n')
    
end