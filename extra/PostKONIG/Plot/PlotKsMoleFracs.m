%% The Function plots the Mole Fractions of the chemical system's component of interest vs its overall rates  
%
%  Input Arguments:  - iT:                      Index for the current Translational Temperature
%                    - t:                       Vector of time instants
%                    - ProcessesRatesOverall:   Matrix of Processes Overall Rates (Levels X 4(i.e.: Diss,Pair1,Pair2,Pair3) X NTint)
%                    - MolFracs:                Matrix of Mole Fractions (Time-Instants X Components)
%
%  Input Global Var: - CompOI:                  Component of interest
%                    - XLimPlot:                Vector of [Min, Max] for x-axes plot
%                    - YLimPlot:                Vector of [Min, Max] for y-axes plot
%                    - PlotPairUser:            Vector of Flags 0/1; if =1 the atomic pair is plotted
%

function [iFigure] = PlotKsMoleFracs(iT, iFigure, t, ProcessesRatesOverall, MolFracs)    
    
  %%==============================================================================================================
  % 
  % Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
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

  global PlotPairUser XLimPlot YLimPlot CompOI 
  global NBinnedMol BinnedMolToComp MoleculesName CompNames CompColor 
  global linST SaveFigs FigDirPath linS AxisFontSz AxisFontNm LegendFontSz AxisLabelSz AxisLabelNm LegendFontNm 
  global RCVec BCVec GCVec KCVec OCVec PCVec WCVec JCVec YCVec CCVec MCVec

  iMol=1;

  figure(iFigure)
  fig = gcf;
  screensize = get( groot, 'Screensize' );
  fig.Position=screensize;
  fig.Color='None';

  yyaxis right
  loglog(t(:),ProcessesRatesOverall(:,1,iT),'LineWidth',3,'linestyle',linST{1});
  ProcessesStr(1) = {'$\bar{K}_{Diss}$'};
  hold on
  iPP=2;
  for iP = 1:3
    if PlotPairUser(iP) == 1
      loglog(t(:),ProcessesRatesOverall(:,iP+1,iT),'LineWidth',3,'linestyle',linST{iP+1});
      ProcessesStr(iPP) = cellstr(strcat('$\bar{K}_{Exch',num2str(iP),'}$'));
      iPP = iPP + 1;
    end
  end


  xt = get(gca, 'XTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  yt = get(gca, 'YTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

  str_x = ['Time [s]'];
  xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  xlab.Interpreter = 'latex';
  xlim(XLimPlot);

  str_y = ['$\bar{K}_{Process} [cm^3/s]$'];
  ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  ylab.Interpreter = 'latex';
  %ylim(YLimPlot);

  clab = legend(ProcessesStr,'Location','Best');
  clab.Interpreter = 'latex';
  set(clab,'FontSize',LegendFontSz,'FontName',LegendFontNm,'Interpreter','latex');

  pbaspect([1 1 1])


  yyaxis left
  for iComp=CompOI
    semilogx(t(:),MolFracs(:,iComp),'linestyle',linS{iComp},'LineWidth',5)
    hold on
  end
  hold on

  xt = get(gca, 'XTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  yt = get(gca, 'YTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

  str_x = ['Time [s]'];
  xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  xlab.Interpreter = 'latex';
  xlim(XLimPlot);

  str_y = ['Mole Fraction'];
  ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  ylab.Interpreter = 'latex';
  ylim(YLimPlot);

  pbaspect([1 1 1])


%   clab = legend(CompNames,'Location','Best');
%   clab.Interpreter = 'latex';
%   set(clab,'FontSize',LegendFontSz,'FontName',LegendFontNm,'Interpreter','latex');

  if SaveFigs == 1
     FolderPath = strcat(FigDirPath, '/FlowQuantities/');
     [status,msg,msgID] = mkdir(FolderPath);
     FileName   = strcat(FolderPath, MoleculesName(1,:),'-MoleFractions_vs_K');
     export_fig(FileName, '-pdf')
     close
  elseif SaveFigs == 2
    FolderPath = strcat(FigDirPath, '/FlowQuantities/');
    [status,msg,msgID] = mkdir(FolderPath);
    FileName   = strcat(FolderPath, MoleculesName(1,:),'-MoleFractions_vs_K.fig');
    savefig(FileName)
    close
  end

  iFigure=iFigure+1;

end