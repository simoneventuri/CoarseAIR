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

T_MW        = linspace(1000.d0,60000.d0,1000);
MoleculeMu  = [28.0104e-3, 31.9988d-3];
% mu_MW       = [];
% theta_MW    = [];
% A_MW_    = 1.16d-3 .* mu_MW(iMol).^(0.5d0) .* theta_MW(iMol).^(4.d0/3.d0);
A_MW     = [  47.7d0,  47.7d0];
% B_MW_    = 0.015d0 * mu_MW(iMol).^(0.25d0);
B_MW     = [ 0.050d0, 0.060d0];
SigmaC   = [  2.d-21, 2.d-22];
MinT     = [20000.d0, 15000.d0]
for iMol = 1:2
  Sigma(:,iMol)     = SigmaC(iMol) * (50000.d0 ./ min(T_MW(:), MinT(iMol))).^(2.d0);
  tauP_PC(:,iMol)   = (8.d0 * UKb * T_MW(:) ./ (pi .* MoleculeMu(iMol) ./ AvN) ).^(0.5d0) .* Sigma(:,iMol) .* (ATMToPa ./ (T_MW(:) .* UKb));
  tauP_MW(:,iMol)   = exp(A_MW(iMol) * (T_MW(:).^-(1.d0/3.d0) - B_MW(iMol)) - 18.42d0);
end
tauP_Corr = tauP_PC.^(-1) + tauP_MW;
% exp(A * (x.^-(1.d0/3.d0) - B) - 18.42d0) + ((8.d0 .* 1.380658e-23 .* x ./ (pi .* 28.0104e-3 ./ 6.0221409e+23) ).^(0.5d0) .* SigmaC .* (50000.d0 ./ min(x, 20000.d0)).^(2.d0) .* (101325 ./ (x .* 1.380658e-23))).^(-1)

T_vec        = [   5000.d0,   7500.d0,    10000.d0,   12500.d0,   15000.d0];

% tau_vib_wo   = [2.7795e-08/0.081756, 2.6380e-08/0.122634, 2.8686e-08/0.1635124, 3.0425e-08/0.20439, 3.3382e-08/0.245268]%./0.05.*0.95;
% tau_vib      = [1.8394e-08/0.081756, 1.8255e-08/0.122634, 1.8979e-08/0.1635124, 2.0032e-08/0.20439, 2.1230e-08/0.245268]%./0.05.*0.95;
% 
% tau_rot_wo   = [5.7701e-09/0.081756, 8.277e-09/0.122634, 1.0552e-08/0.1635124, 1.2675e-08/0.20439, 1.5210e-08/0.245268]%./0.05.*0.95;
% tau_rot      = [5.2758e-09/0.081756, 7.023e-09/0.122634, 8.9550e-09/0.1635124, 1.0702e-08/0.20439, 1.2342e-08/0.245268]%./0.05.*0.95;

CO_tau_int_wo   = [1.7587e-08, 1.8996e-08, 2.1387e-08, 2.3724e-08, 2.7226e-08];
CO_tau_int      = [1.2664e-08, 1.3613e-08, 1.5400e-08, 1.6915e-08, 1.8013e-08];

CO_tau_vib_wo   = [2.7795e-08, 2.6380e-08, 2.8686e-08, 3.0425e-08, 3.3382e-08];
CO_tau_vib      = [1.8394e-08, 1.8255e-08, 1.8979e-08, 2.0032e-08, 2.1230e-08];

CO_tau_rot_wo   = [5.7701e-09, 8.277e-09, 1.0552e-08, 1.2675e-08, 1.5210e-08];
CO_tau_rot      = [5.2758e-09, 7.023e-09, 8.9550e-09, 1.0702e-08, 1.2342e-08];

CO_T_Park       = [ 0.0516e00 0.0996e00 ];
CO_tau_vib_Park = [1.0889e-08 1.1396e-07];

T_Center    = [6.28125e-2 6.33398e-2 6.68555e-2 6.75000e-2 6.79687e-2 6.79687e-2 6.82617e-2 6.83789e-2 7.20703e-2 7.28906e-2 7.32422e-2 7.34766e-2 7.32422e-2 7.29492e-2 ...
               7.41797e-2 7.40625e-2 7.41211e-2 7.80469e-2 7.88672e-2 7.95117e-2 7.94531e-2 7.94531e-2 8.04492e-2];
tauP_Center = [1.63784e-8 1.82213e-8 2.44531e-8 1.96846e-8 1.88351e-8 1.88351e-8 2.89668e-8 2.43729e-8 3.41365e-8 3.79535e-8 3.95953e-8 3.12220e-8 2.46443e-8 2.07543e-8 ...
               3.70412e-8 3.32654e-8 2.79934e-8 4.86128e-8 5.40486e-8 5.51495e-8 4.64209e-8 3.50775e-8 3.98344e-8];
% T_Center    = linspace(1000.d0,60000.d0,1000); 
% tauP_Center = exp(54.d0 * T_Center.^(-1.d0/3.d0) - 7.3d0) * 1.d-6;


figure(iFigure)
fig = gcf;
screensize = get( groot, 'Screensize' );
fig.Position=screensize;
fig.Color='None';

%h1=semilogy(T_vec.^(-1/3), CO_tau_int_wo,'+','Color',KCVec,'LineWidth',3,'MarkerSize',8,'MarkerFaceColor',KCVec);
%hold on 
%h2=semilogy(T_vec.^(-1/3), CO_tau_int,'+-','Color',KCVec,'LineWidth',3,'MarkerSize',8,'MarkerFaceColor',KCVec);
h3=semilogy(T_vec.^(-1/3), CO_tau_vib_wo,'x','Color',KCVec,'LineWidth',4,'MarkerSize',12,'MarkerFaceColor',KCVec);
hold on
h4=semilogy(T_vec.^(-1/3), CO_tau_vib,'x-','Color',KCVec,'MarkerSize',12,'MarkerFaceColor',KCVec,'LineWidth',2);
%h5=semilogy(CO_T_Park, CO_tau_vib_Park,':','Color',KCVec,'LineWidth',3);
h9=semilogy(T_MW(:).^(-1/3), tauP_Corr(:,1),'.-','Color',BCVec,'MarkerSize',8,'MarkerFaceColor',BCVec,'LineWidth',3);
h8=semilogy(T_MW(:).^(-1/3), tauP_MW(:,1),':','Color',BCVec,'MarkerSize',8,'MarkerFaceColor',BCVec,'LineWidth',2);
h8=semilogy(T_Center, tauP_Center,'^','Color',GCVec,'MarkerSize',8,'MarkerFaceColor',GCVec,'LineWidth',2);
h6=semilogy(T_vec.^(-1/3), CO_tau_rot_wo,'o','Color',KCVec,'MarkerSize',10,'MarkerFaceColor',KCVec);
h7=semilogy(T_vec.^(-1/3), CO_tau_rot,'o-','Color',KCVec,'MarkerSize',8,'MarkerFaceColor',KCVec,'LineWidth',2);

%LegendText = {'CO, $\tau^{IT}$ w/o exchange'; 'CO, $\tau^{IT}$ with exchange';'CO, $\tau^{VT}$ w/o exchange'; 'CO, $\tau^{VT}$ with exchange'; 'CO, $\tau^{VT}$ with exchange by Park et al.';'CO, $\tau^{RT}$ w/o exchange';'CO, $\tau^{RT}$ with exchange'};
LegendText = {'$\tau^{VT}$ w/o Exchange'; '$\tau^{VT}$ with Exchange'; '$\tau^{VT}$, Park 1994 ($\sigma_v=2 \times 10^{-21} m^2$)'; '$\tau^{VT}$, Millikan and White 1963'; '$\tau^{VT}$, Center 1973 '; '$\tau^{RT}$ w/o Exchange';'$\tau^{RT}$ with Exchange'};
clab = legend(LegendText,'Location','Best');
clab.Interpreter = 'latex';
set(clab,'FontSize',LegendFontSz,'FontName',LegendFontNm,'Interpreter','latex');

xt = get(gca, 'XTick');
set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
yt = get(gca, 'YTick');
set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

str_x = ['$T^{-1/3} [K^{-1/3}]$'];
xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
xlab.Interpreter = 'latex';
xlim([0.02, 0.09]);

str_y = ['$\tau p_O$ [atm-s]'];
ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
ylab.Interpreter = 'latex';
ylim([1.d-9, 1.d-5]);

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



T_vec        = [   5000.d0,   7500.d0,    10000.d0,   12500.d0,   15000.d0];

% O2_tau_vib = [3.5645e-08*0.08175620528, 6.0851e-08*0.1226343079, 9.1425e-08*0.1635124106, 0.0, 1.217e-07*0.2452686158];
% 
% O2_tau_rot = [1.0693e-08*0.08175620528, 1.4813e-08*0.1226343079, 1.9750e-08*0.1635124106, 0.0, 3.0828e-08*0.2452686158];

O2_tau_int = [2.3945e-08, 3.8106e-08, 5.532e-08, 7.2677e-08, 8.298e-08];
%O2_tau_int = [3.5645e-08, 6.1708e-08, 8.955e-08, 0.0, 1.186e-07];

O2_tau_vib = [3.622e-08, 6.2983e-08, 9.1891e-08, 1.0929e-07, 1.1873e-07];
%O2_tau_vib = [3.5645e-08, 6.1708e-08, 8.955e-08, 0.0, 1.186e-07];

O2_tau_rot = [1.078e-08, 1.5242e-08, 2.0049e-08, 2.5062e-08, 3.0074e-08];
%O2_tau_rot = [1.0693e-08, 1.5020e-08, 2.0027e-08, 0.0, 3.0041e-08];

O2_T_Park       = [0.059489671931956256,  0.10068043742405833];
O2_tau_vib_Park = [9.886051135385972e-9, 7.316273646063192e-8];


%T_Adr  = [9.88539e+2 2.00573e+3 2.97994e+3 4.98567e+3 1.00000e+4 3.98281e+3 6.64756e+3 8.23782e+3];
%O2_Adr = [2.80826e-8 3.92419e-8 4.86968e-8 6.51722e-8 9.20582e-8 5.54307e-8 7.26001e-8 8.00064e-8];
T_Adr = [ 4.83077e+2, 6.97350e+2, 9.83233e+2, 1.48399e+3, 1.99924e+3, 2.98711e+3, 3.97517e+3, 4.96333e+3, 5.96586e+3, 1.00050e+4];
O2_Adr  = [ 7.05403e-9, 9.54306e-9, 1.30504e-8, 1.80406e-8, 2.31240e-8, 3.19659e-8, 4.00982e-8, 4.81740e-8, 5.66402e-8, 9.20582e-8];


figure(iFigure)
fig = gcf;
screensize = get( groot, 'Screensize' );
fig.Position=screensize;
fig.Color='None';
 
%h11=semilogy(T_vec.^(-1/3),O2_tau_int,'+-','Color',KCVec,'MarkerSize',8,'MarkerFaceColor',KCVec,'LineWidth',3);
%hold on
h12=semilogy(T_vec.^(-1/3),O2_tau_vib,'x-','Color',KCVec,'MarkerSize',12,'MarkerFaceColor',KCVec,'LineWidth',3);
hold on
h15=semilogy(T_MW(:).^(-1/3), tauP_Corr(:,2),'-','Color',BCVec,'MarkerSize',8,'MarkerFaceColor',BCVec,'LineWidth',3);
h16=semilogy(T_MW(:).^(-1/3), tauP_MW(:,2),':','Color',BCVec,'MarkerSize',8,'MarkerFaceColor',BCVec,'LineWidth',3);
h17=semilogy(T_Adr(:).^(-1/3), O2_Adr(:),'^','Color',RCVec,'MarkerSize',8,'MarkerFaceColor',RCVec,'LineWidth',2);
h14=semilogy(T_vec.^(-1/3), O2_tau_rot,'o-','Color',KCVec,'MarkerSize',8,'MarkerFaceColor',KCVec,'LineWidth',3);


LegendText = {'$\tau^{VT}$'; '$\tau^{VT}$, Park 1994 ($\sigma_v=2 \times 10^{-22} m^2$)'; '$\tau^{VT}$, Millikan and White 1963'; '$\tau^{VT}$ for $O_2+O$, Andrienko 2016'; '$\tau^{RT}$'};
clab = legend(LegendText,'Location','Best');
clab.Interpreter = 'latex';
set(clab,'FontSize',LegendFontSz,'FontName',LegendFontNm,'Interpreter','latex');

xt = get(gca, 'XTick');
set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
yt = get(gca, 'YTick');
set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

str_x = ['$T^{-1/3} [K^{-1/3}]$'];
xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
xlab.Interpreter = 'latex';
%xlim(XLimPlot);

str_y = ['$\tau p_O$ [atm-s]'];
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