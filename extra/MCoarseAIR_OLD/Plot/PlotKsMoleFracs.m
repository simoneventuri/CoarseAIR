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

function [iFigure] = PlotKsMoleFracs(iT, iFigure, t, ProcessesRatesOverall, MolFracs, tStart, iQSS, tEnd)  
    
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

  global PlotPairUser XLimPlot YLimPlot CompOI 
  global NBinnedMol BinnedMolToComp MoleculesName CompNames CompColor 
  global linST SaveFigs FigDirPath linS AxisFontSz AxisFontNm LegendFontSz AxisLabelSz AxisLabelNm LegendFontNm 
  global RCVec BCVec GCVec KCVec OCVec PCVec WCVec JCVec YCVec CCVec 
  global DissFlg InelFlg ExchFlg T0_Vec

  iMol=1;
  
 
  figure(iFigure)
  fig = gcf;
  %screensize = get( groot, 'Screensize' );
  %fig.Position=screensize;
  %fig.Color='None';

  yyaxis right
  loglog(t(:),ProcessesRatesOverall(:,1,iT),'LineWidth',3,'linestyle',linST{1});
  ProcessesStr(1) = {'$\bar{K}_{Diss}$'};
  hold on
%   iPP=2;
%   for iP = 1:3
%     if PlotPairUser(iP) == 1
%       loglog(t(:),ProcessesRatesOverall(:,iP+1,iT),'LineWidth',3,'linestyle',linST{iP+1});
%       ProcessesStr(iPP) = cellstr(strcat('$\bar{K}_{Exch',num2str(iP),'}$'));
%       iPP = iPP + 1;
%     end
%   end
  
  yy = ProcessesRatesOverall(:,1,iT);
  SStart = 1;
  while (yy(SStart) < yy(iQSS)/10)
    SStart = SStart + 1;
  end
  EEnd = iQSS;
  while (abs(yy(EEnd) - yy(iQSS))/yy(iQSS) > 3.0)
    EEnd = EEnd + 1;
  end
  EEnd=length(yy);
  yyy(SStart:EEnd)  = ProcessesRatesOverall(SStart:EEnd,1,iT);
  [xData, yData] = prepareCurveData( t(SStart:EEnd), yyy(SStart:EEnd) );
  ft = 'splineinterp';
  [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
  ProcessesRatesOverall_Start  = fitresult(tStart);
  ProcessesRatesOverall_End    = fitresult(tEnd);
  loglog(tStart, ProcessesRatesOverall_Start, 'ko','MarkerSize',8,'MarkerFaceColor','k');
  loglog(tEnd,   ProcessesRatesOverall_End,   'ko','MarkerSize',8,'MarkerFaceColor','k');

  clear yyy 
  yyy(SStart:EEnd) = MolFracs(SStart:EEnd,CompOI);
  [xData, yData]   = prepareCurveData( t(SStart:EEnd), yyy(SStart:EEnd) );
  ft = 'splineinterp';
  [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );  
  Mol_Initial = MolFracs(1,CompOI);
  Mol_Final   = MolFracs(end,CompOI);
  Mol_Start   = fitresult(tStart);
  Mol_End     = fitresult(tEnd);
  Delta_Mol   = Mol_Initial - Mol_Final;
  Delta_Pre   = Mol_Initial - Mol_Start;
  Delta_Post  = Mol_End     - Mol_Final;
  Perc_Pre    = Delta_Pre  / Delta_Mol * 100.0
  Perc_Post   = Delta_Post / Delta_Mol * 100.0
  
  
  
  FileName = strcat('./PercDissQSS_',num2str(DissFlg),'_',num2str(InelFlg),'_',num2str(ExchFlg),'_',num2str(ExchFlg),'.csv');
  if exist(FileName, 'file')
    fileID1  = fopen(FileName,'a');
  else
    fileID1  = fopen(FileName,'w');
    fprintf(fileID1,'# T [K], QSS_Value, Perc_Pre, Perc_Post\n');
  end
  fprintf(fileID1,'%e,%e,%e,%e\n', T0_Vec(1), ProcessesRatesOverall(iQSS,1,iT), Perc_Pre, Perc_Post )
  fclose(fileID1);
  
  
  xt = get(gca, 'XTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  yt = get(gca, 'YTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

  str_x = ['Time [s]'];
  xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  xlab.Interpreter = 'latex';
  xlim(XLimPlot);

  str_y = ['$\bar{K}^D$ $[cm^3/s]$'];
  ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  ylab.Interpreter = 'latex';
  %ylim(YLimPlot);

%   clab = legend(ProcessesStr,'Location','Best');
%   clab.Interpreter = 'latex';
%   set(clab,'FontSize',LegendFontSz,'FontName',LegendFontNm,'Interpreter','latex');

  pbaspect([1 1 1])


  yyaxis left
  for iComp=CompOI
    semilogx(t(:),MolFracs(:,iComp),'linestyle',linS{iComp},'LineWidth',5)
    hold on
  end
  hold on
  
  semilogx([tStart, tStart],   [0, 1], ':k', 'LineWidth',2)
  semilogx([tEnd,     tEnd],   [0, 1], ':k', 'LineWidth',2)
  
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