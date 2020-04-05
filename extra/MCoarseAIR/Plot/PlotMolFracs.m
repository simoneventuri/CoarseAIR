%% The Function plots the Mole Fractions of the chemical system's components 
%
%  Input Arguments:  - iT:               Index for the current Translational Temperature
%                    - t:                Vector of time instants
%                    - MolFracs:         Matrix of Mole Fractions (Time-Instants X Components)
%
%  Input Global Var: - XLimPlot:         Vector of [Min, Max] for x-axes plot
%                    - YLimPlot:         Vector of [Min, Max] for y-axes plot
%

function i_figure = PlotMolFracs(iT, i_figure, t, MolFracs)    
    
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

  global CompOI NComp CompNames CompColor SaveFigs FigDirPath linS AxisFontSz AxisFontNm LegendFontSz LegendFontNm AxisLabelSz AxisLabelNm XLimPlot YLimPlot 

  figure(i_figure)
  fig = gcf;
  screensize = get( groot, 'Screensize' );
  fig.Position=screensize;
  %fig.Color='None';

  %for iComp=1:NComp
      iComp=CompOI
      semilogx(t(:),MolFracs(:,iComp),'Color',CompColor(iComp,:),'linestyle',linS{iComp},'LineWidth',4)
      hold on
  %end
  hold on

  xt = get(gca, 'XTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  yt = get(gca, 'YTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

  clab = legend(CompNames,'Location','Best');
  clab.Interpreter = 'latex';
  set(clab,'FontSize',LegendFontSz,'FontName',LegendFontNm,'Interpreter','latex');

  str_x = ['Time [s]'];
  xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  xlab.Interpreter = 'latex';
  %xlim(XLimPlot);

  str_y = ['Mole Fraction'];
  ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  ylab.Interpreter = 'latex';
  %ylim(YLimPlot);

  pbaspect([1 1 1])

  if SaveFigs == 1
     FolderPath = strcat(FigDirPath, '/FlowQuantities/');
     [status,msg,msgID] = mkdir(FolderPath);
     FileName   = strcat(FolderPath, 'MoleFractions');
     export_fig(FileName, '-pdf')
     close
  elseif SaveFigs == 2
    FolderPath = strcat(FigDirPath, '/FlowQuantities/');
    [status,msg,msgID] = mkdir(FolderPath);
    FileName   = strcat(FolderPath, 'MoleFractions.fig');
    savefig(FileName)
    close
  end

  i_figure=i_figure+1;

end