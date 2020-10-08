close all
clear all
clc


global Param Input Syst
iFig = 1;

% Input.SystNameLong       = 'O3_UMN'
% NBins                    = 45
% Input.SystNameLong       = 'N3_NASA'
% NBins                    = 61
% Input.SystNameLong       = 'CO2_NASA'
% NBins                    = 83
% Input.SystNameLong       = 'O2C_NASA'
% NBins                    = 83
% Input.SystNameLong       = 'N2O_UMN'
% NBins                    = 54
Input.SystNameLong       = 'NON_UMN'
NBins                    = 54
Input.FigureFormat       = 'PrePrint';

Input.iFig               = 101;
Input.SaveFigsFlgInt     = 2;
Input.Paths.SaveFigsFldr = strcat('/home/venturi/WORKSPACE/Air_Paper/Figures/', Input.SystNameLong);

Syst.NameLong = Input.SystNameLong
Syst = Initialize_ChemicalSyst(Syst)
Initialize_Parameters()


Input.NTaus         = 4
Input.Paths.Taus(1) = {strcat('/home/venturi/WORKSPACE/Air_Paper/Data/',Input.SystNameLong,'/Taus_',Syst.Molecule(1).Name,'_0_1_0_1.csv')}
Input.Paths.Taus(2) = {strcat('/home/venturi/WORKSPACE/Air_Paper/Data/',Input.SystNameLong,'_VSM/Taus_',Syst.Molecule(1).Name,'_0_1_0_1.csv')}
Input.Paths.Taus(3) = {strcat('/home/venturi/WORKSPACE/Air_Paper/Data/',Input.SystNameLong,'_RVE',num2str(NBins),'/Taus_',Syst.Molecule(1).Name,'_0_1_0_1.csv')}
Input.Paths.Taus(4) = {strcat('/home/venturi/WORKSPACE/Air_Paper/Data/',Input.SystNameLong,'_ADA',num2str(NBins),'/Taus_',Syst.Molecule(1).Name,'_0_1_0_1.csv')}
Input.TauNames      = [{'StS'}, {strcat('VS',num2str(NBins))}, {strcat('RVE',num2str(NBins))}, {strcat('ADA',num2str(NBins))}]

CMat      = [Param.KCVec; Param.RCVec; Param.GCVec;  Param.BCVec]; 
LineWidth = 3;

for iTau = 1:Input.NTaus
    
    opts = delimitedTextImportOptions("NumVariables", 5);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["T", "P", "tau_Int", "tau_Rot", "tau_Vib"];
    opts.VariableTypes = ["double", "double", "double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    tbl = readtable(char(Input.Paths.Taus(iTau)), opts);
    T_NoExch       = tbl.T;
    P_NoExch       = tbl.P;
    tau_Int_NoExch = tbl.tau_Int;
    tau_Rot_NoExch = tbl.tau_Rot;
    tau_Vib_NoExch = tbl.tau_Vib;
    clear opts tbl
    
    
    figure(Input.iFig)
    fig = gcf;
    screensize   = get( groot, 'Screensize' );
    %fig.Position = screensize;
    %fig.Color='None';
    semilogy(T_NoExch.^(-1/3), tau_Vib_NoExch, 'Color', CMat(iTau,:), 'linestyle', char(Param.linS(iTau)), 'LineWidth', LineWidth)
    hold on
    
    
    figure(Input.iFig+1)
    fig = gcf;
    screensize   = get( groot, 'Screensize' );
    %fig.Position = screensize;
    %fig.Color='None';
    semilogy(T_NoExch.^(-1/3), tau_Rot_NoExch, 'Color', CMat(iTau,:), 'linestyle', char(Param.linS(iTau)), 'LineWidth', LineWidth)
    hold on

    
    figure(Input.iFig+2)
    fig = gcf;
    screensize   = get( groot, 'Screensize' );
    %fig.Position = screensize;
    %fig.Color='None';
    semilogy(T_NoExch.^(-1/3), tau_Int_NoExch, 'Color', CMat(iTau,:), 'linestyle', char(Param.linS(iTau)), 'LineWidth', LineWidth)
    hold on

end

figure(Input.iFig)
fig = gcf;
screensize   = get( groot, 'Screensize' );
%fig.Position = screensize;
%fig.Color='None';

xt = get(gca, 'XTick');
set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
yt = get(gca, 'YTick');
set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

clab             = legend(Input.TauNames, 'Location', 'Best');
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




figure(Input.iFig+1)
fig = gcf;
screensize   = get( groot, 'Screensize' );
%fig.Position = screensize;
%fig.Color='None';

xt = get(gca, 'XTick');
set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
yt = get(gca, 'YTick');
set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

clab             = legend(Input.TauNames, 'Location', 'Best');
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





figure(Input.iFig+2)
fig = gcf;
screensize   = get( groot, 'Screensize' );
%fig.Position = screensize;
%fig.Color='None';

xt = get(gca, 'XTick');
set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
yt = get(gca, 'YTick');
set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

clab             = legend(Input.TauNames, 'Location', 'Best');
clab.Interpreter = 'latex';
set(clab,'FontSize', Param.LegendFontSz, 'FontName', Param.LegendFontNm, 'Interpreter', 'latex');

str_x = ['T${}^{-1/3}$ [K${}^{-1/3}]$'];
xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
xlab.Interpreter = 'latex';
%xlim(XLimPlot);

str_y = ['$p\tau_I$ [atm$\cdot$s]'];
ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
ylab.Interpreter = 'latex';
%ylim(YLimPlot);

if Input.SaveFigsFlgInt > 0
    [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
    FolderPath = strcat(Input.Paths.SaveFigsFldr, '/');
    [status,msg,msgID] = mkdir(FolderPath);
    if Input.SaveFigsFlgInt == 1
        FileName   = strcat(FolderPath, 'TauIntS');
        export_fig(FileName, '-pdf')
    elseif Input.SaveFigsFlgInt == 2
        FileName   = strcat(FolderPath, 'TauIntS.fig');
        savefig(FileName)
    end
    %close
end