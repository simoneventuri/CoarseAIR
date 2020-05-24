close all
clear all
clc


global Param Input
iFig = 1;

Input.SystNameLong       = 'CO2_NASA'
Input.Paths.TausWOExch   = '/home/venturi/WORKSPACE/Mars_Paper/Data/CO2_NASA/CO2_NASA/Taus_CO_0_1_0_0.csv'
Input.Paths.TausWExch    = '/home/venturi/WORKSPACE/Mars_Paper/Data/CO2_NASA/CO2_NASA/Taus_CO_0_1_1_0.csv'

Input.FigureFormat       = 'PrePrint';

Input.iFig               = 101;
Input.SaveFigsFlgInt     = 0;
Input.Paths.SaveFigsFldr = '/home/venturi/WORKSPACE/Mars_Paper/Figures/';


Initialize_ChemicalSyst()
Initialize_Parameters()


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


fprintf('= Plot_Taus ========================================\n')
fprintf('====================================================\n')


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

str_x = ['T${}^{-1/3}$ [1/K${}^{1/3}]$'];
xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
xlab.Interpreter = 'latex';
%xlim(XLimPlot);

str_y = ['$\tau_v$ atm*s'];
ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
ylab.Interpreter = 'latex';
%ylim(YLimPlot);


if Input.SaveFigsFlgInt > 0
    [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
    FolderPath = strcat(Input.Paths.SaveFigsFldr, '/');
    [status,msg,msgID] = mkdir(FolderPath);
    if Input.SaveFigsFlgInt == 1
        FileName   = strcat(FolderPath, 'GlobalRates');
        export_fig(FileName, '-pdf')
    elseif Input.SaveFigsFlgInt == 2
        FileName   = strcat(FolderPath, 'GlobalRates.fig');
        savefig(FileName)
    end
    close
end
Input.iFig = Input.iFig + 1;


fprintf('====================================================\n\n')