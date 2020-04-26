close all
clc

global MoleculeMu mu_MW theta_MW
global ATMToPa UKb AvN
global SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz AxisLabelSz AxisLabelNm LegendFontNm 
global RCVec BCVec GCVec KCVec OCVec PCVec WCVec JCVec YCVec CCVec MCVec

iFigure  = 1;
SaveFigs = 0;

linS     = {'--','-.',':','-'};
linST    = {'-','-.',':','--'};

AxisFontSz = 36;
AxisFontNm = 'Arial';

AxisLabelSz = 44;
AxisLabelNm = 'Arial';

LegendFontSz = 40;
LegendFontNm = 'Arial';

RCVec = [255  50  20] ./ 255;
BCVec = [  0  70 200] ./ 255;
GCVec = [  0 140  50] ./ 255;
KCVec = [  0   0   0] ./ 255;
OCVec = [255 105  45] ./ 255;
PCVec = [155  45 175] ./ 255;
WCVec = [  1   1   1] ./ 255;
JCVec = [100 100 100] ./ 255;
YCVec = [255 255   0] ./ 255;
CCVec = [205 205 205] ./ 255;
MCVec = [100  25  15] ./ 255;


T_vec        = [   5000.d0,   7500.d0,    10000.d0,   12500.d0,   15000.d0];

K_Exch     = [1.071e-16; 2.299e-14; 3.650e-13; 1.834e-12; 5.260e-12];
K_Exch_QSS = [1.071e-16; 2.299e-14; 3.351e-13; 1.581e-12; 4.199e-12];

K_Diss     = [1.278e-18; 3.563e-15; 1.854e-13; 1.861e-12; 8.232e-12];
K_Diss_QSS = [1.278e-18; 1.169e-15; 5.171e-14; 4.622e-13; 1.865e-12];

T_Park = [0.5, 1.785166240409207];
K_Park = [4.4173447031400823e-11, 1.e-17];

T_Appleton = [0.6649616368286446, 1.2480818414322252]; 
K_Appleton = [1e-11, 3.2819278725114844e-14];

T_Johnston = [0.5767263427109974, 1.9117647058823528]; 
K_Johnston = [1e-10, 1.1178591777554e-17];

T_Hanson = [0.8529411764705882, 1.819693094629156]; 
K_Hanson = [9.284145445194764e-13, 2.2638034095214513e-16];


figure(iFigure)
fig = gcf;
screensize = get( groot, 'Screensize' );
fig.Position=screensize;
fig.Color='None';

h55=semilogy(T_Park, K_Park,':','Color',RCVec,'LineWidth',3);
hold on
h66=semilogy(T_Appleton, K_Appleton,'-.','Color',GCVec,'LineWidth',3);
h77=semilogy(T_Hanson, K_Hanson,'-','Color',BCVec,'LineWidth',3);
h88=semilogy(T_Johnston, K_Johnston,'--','Color',CCVec,'LineWidth',3);
h11=semilogy(10000.d0./T_vec,K_Diss,'x-','Color',KCVec,'MarkerSize',12,'MarkerFaceColor',KCVec,'LineWidth',2);
h22=semilogy(10000.d0./T_vec,K_Diss_QSS,'x','Color',KCVec,'MarkerSize',12,'MarkerFaceColor',KCVec,'LineWidth',2);
h33=semilogy(10000.d0./T_vec,K_Exch+K_Diss,'o-','Color',KCVec,'MarkerSize',8,'MarkerFaceColor','w','LineWidth',2);
h44=semilogy(10000.d0./T_vec,K_Exch_QSS+K_Diss_QSS,'o','Color',KCVec,'MarkerSize',8,'MarkerFaceColor','w','LineWidth',2);

LegendText = {'Park94';'Appleton70 2-par';'Hanson74';'Johnston2014';'QCT, $K_{Diss}$';'QCT, $K_{Diss_{QSS}}$';'QCT, $K_{(Diss+Exch)}$';'QCT, $K_{(Diss+Exch)_{QSS}}$'};
clab = legend(LegendText,'Location','Best');
clab.Interpreter = 'latex';
set(clab,'FontSize',LegendFontSz,'FontName',LegendFontNm,'Interpreter','latex');

xt = get(gca, 'XTick');
set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
yt = get(gca, 'YTick');
set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

str_x = ['$10000/T [K^{-1}]$'];
xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
xlab.Interpreter = 'latex';
%xlim(XLimPlot);

str_y = ['K $[cm^3/s]$'];
ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
ylab.Interpreter = 'latex';
%ylim(YLimPlot);

if SaveFigs == 1
  FolderPath = strcat(FigDirPath, '/FlowQuantities/');
  [status,msg,msgID] = mkdir(FolderPath);
  FileName   = strcat(FolderPath, MoleculesName(1,:),'-Energies');
  export_fig(FileName, '-pdf')
  close
elseif SaveFigs == 2
  FolderPath = strcat(FigDirPath, '/FlowQuantities/');
  [status,msg,msgID] = mkdir(FolderPath);
  FileName   = strcat(FolderPath, MoleculesName(1,:),'-Energies.fig');
  savefig(FileName)
  close
end

iFigure=iFigure+1;   




T_vec      = [   5000.d0,   7500.d0,    10000.d0,   12500.d0,   15000.d0];
K_Exch     = [4.6e-10; 5.219e-10; 6.287e-10; 7.209e-10; 7.944e-10];
K_Diss     = [2.3e-13; 1.233e-11; 5.868e-11; 1.424e-10; 2.504e-10];

T_Park = [0.5, 1.785166240409207];
K_Park = [4.4173447031400823e-11, 1.e-17];

T_Appleton = [0.6649616368286446, 1.2480818414322252]; 
K_Appleton = [1e-11, 3.2819278725114844e-14];

T_Johnston = [0.5767263427109974, 1.9117647058823528]; 
K_Johnston = [1e-10, 1.1178591777554e-17];

T_Hanson = [0.8529411764705882, 1.819693094629156]; 
K_Hanson = [9.284145445194764e-13, 2.2638034095214513e-16];

T_Park = linspace(5000,15000,1000);
K_Park = 1.0d22 * (T_Park).^(-3.d0/2.d0) .* exp( - 59750.d0 ./ T_Park) ./ AvN;
%K_Park = 5.05d21 * (T_Park).^(-1.5) .* exp( - 41500.d0 ./ T_Park) ./ AvN;

figure(iFigure)
fig = gcf;
screensize = get( groot, 'Screensize' );
fig.Position=screensize;
fig.Color='None';

h55=semilogy(10000.d0./T_Park, K_Park,':','Color',RCVec,'LineWidth',3);
hold on
h11=semilogy(10000.d0./T_vec,K_Diss,'x-','Color',KCVec,'MarkerSize',12,'MarkerFaceColor',KCVec,'LineWidth',2);
h33=semilogy(10000.d0./T_vec,K_Exch+K_Diss,'o-','Color',KCVec,'MarkerSize',8,'MarkerFaceColor','w','LineWidth',2);

LegendText = {'Park94';'QCT, $K_{Diss}$';'QCT, $K_{(Diss+Exch)}$'};
clab = legend(LegendText,'Location','Best');
clab.Interpreter = 'latex';
set(clab,'FontSize',LegendFontSz,'FontName',LegendFontNm,'Interpreter','latex');

xt = get(gca, 'XTick');
set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
yt = get(gca, 'YTick');
set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

str_x = ['$10000/T [K^{-1}]$'];
xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
xlab.Interpreter = 'latex';
xlim([0.6 2.1]);

str_y = ['K $[cm^3/s]$'];
ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
ylab.Interpreter = 'latex';
%ylim(YLimPlot);

if SaveFigs == 1
  FolderPath = strcat(FigDirPath, '/FlowQuantities/');
  [status,msg,msgID] = mkdir(FolderPath);
  FileName   = strcat(FolderPath, MoleculesName(1,:),'-Energies');
  export_fig(FileName, '-pdf')
  close
elseif SaveFigs == 2
  FolderPath = strcat(FigDirPath, '/FlowQuantities/');
  [status,msg,msgID] = mkdir(FolderPath);
  FileName   = strcat(FolderPath, MoleculesName(1,:),'-Energies.fig');
  savefig(FileName)
  close
end

iFigure=iFigure+1;   