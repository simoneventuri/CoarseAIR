%% The Function plots the Mole Fractions of the chemical system's components 
%
%  Input Arguments:  - iT:               Index for the current Translational Temperature
%                    - t:                Vector of time instants
%                    - MolFracs:         Matrix of Mole Fractions (Time-Instants X Components)
%
%  Input Global Var: - XLimPlot:         Vector of [Min, Max] for x-axes plot
%                    - YLimPlot:         Vector of [Min, Max] for y-axes plot
%

function iFigure = PlotMolFracsStoch(iT, iFigure, t, MolFracsStoch)    
    
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

  global NComp CompNames CompColor SaveFigs FigDirPath linS AxisFontSz AxisFontNm LegendFontSz LegendFontNm AxisLabelSz AxisLabelNm ...
         XLimPlot YLimPlot PlotPColorFlg PlotSamplesFlg PlotMomentsFlg NPESs iPESStart iPESEnd PlotPostCut TimePostCutVec

  iComp    = 2;
  tStart   = 1.e-8;%;
  tEnd     = 1.e-3;%;
  NPESTemp = iPESEnd - iPESStart + 1; 

  
   it=1;
    while (t(it) < tStart)
      it=it+1;
    end
    itStart=it;
    fprintf(strcat('  Starting Time = ', num2str(t(itStart)), 's'))
    while (t(it) < tEnd)
      it=it+1;
    end
    itEnd=it;
    fprintf(strcat('; Final Time = ', num2str(t(itEnd)), 's \n'))
  
    
  %%
  if PlotPColorFlg == 1
    
    x(1:(itEnd-itStart+1))            = t(itStart:itEnd);
    Y(1:(itEnd-itStart+1),1:NPESTemp) = squeeze(MolFracsStoch(itStart:itEnd,iComp,iPESStart:iPESEnd));
    
    figure(iFigure)
    fig = gcf;
    screensize = get( groot, 'Screensize' );
    fig.Position=screensize;
    %fig.Color='None';
    
    [iFigure, yHist] = DensityPlotOnlyY(iFigure, x, Y, 300, -1);
  
    set(gca, 'XScale', 'log');
    xt = get(gca, 'XTick');
    set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

    str_x = ['Time [s]'];
    xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
    xlab.Interpreter = 'latex';
    xlim([tStart, tEnd]);

    str_y = ['Mole Fraction'];
    ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
    ylab.Interpreter = 'latex';
    ylim(YLimPlot);

    %pbaspect([1 1 1])

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

    iFigure=iFigure+1;
    
  end
  %%
  
  
  %%     
  if PlotSamplesFlg == 1 || PlotMomentsFlg == 1     
       
    figure(iFigure)
    fig = gcf;
    screensize = get( groot, 'Screensize' );
    fig.Position=screensize;
    %fig.Color='None';

    Summ    = zeros(itEnd-itStart+1,1);
    SummSqr = zeros(itEnd-itStart+1,1);
        
    for iPES=iPESStart:iPESEnd
      %for iComp=1:2%NComp
        
        if PlotMomentsFlg == 1
          Summ(1:itEnd-itStart+1)    = Summ(1:itEnd-itStart+1)    + MolFracsStoch(itStart:itEnd,iComp,iPES);
          SummSqr(1:itEnd-itStart+1) = SummSqr(1:itEnd-itStart+1) + MolFracsStoch(itStart:itEnd,iComp,iPES).^2;
        end
          
        if PlotSamplesFlg == 1
          semilogx(t(:),MolFracsStoch(:,iComp,iPES),'k','LineWidth',1)
          hold on
        end
        
      %end
    end
    
    if PlotMomentsFlg == 1
      Meann  = Summ ./ NPESTemp;
      StDevv = sqrt(  SummSqr./NPESTemp - Meann.^2);

      semilogy(t(itStart:itEnd),Meann(:),'ro');
      hold on
      semilogy([t(itStart:itEnd)'; t(itStart:itEnd)'],[max(Meann'-3.0.*StDevv',Meann'.*0.0+1.d-100); Meann'+3.0.*StDevv'],'r-');
    end

    xt = get(gca, 'XTick');
    set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

    str_x = ['Time [s]'];
    xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
    xlab.Interpreter = 'latex';
    xlim([tStart, tEnd]);

    str_y = ['Mole Fraction'];
    ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
    ylab.Interpreter = 'latex';
    ylim(YLimPlot);

    %pbaspect([1 1 1])

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

    iFigure=iFigure+1;
    
  end
  %%
  
  
  %%
  if PlotPostCut == 1
    
    for Time = TimePostCutVec
      Time
      
      iTime = 1;
      while (t(iTime) < Time)
        iTime=iTime+1;
      end
      
      fig1=figure(iFigure);
      fig = gcf;
      if SaveFigs > 0
        screensize = get( groot, 'Screensize' );
        fig.Position=screensize;
        fig.Color='None';
      end 

      MoleFracsVec = squeeze(MolFracsStoch(iTime,iComp,:));
      histfit(MoleFracsVec,30)

      xlab = xlabel('$Mole Fraction$');
      xlab.Interpreter = 'latex';
      ylab = ylabel('Posterior Samples');
      grid on
      title(strcat('Time ', num2str(Time), 's'))

      set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
      %daspect([1 1 1])

      if SaveFigs == 1
         FolderPath = strcat(FigDirPath, '/Rates/');
         [status,msg,msgID] = mkdir(FolderPath);
         FileName   = strcat(FolderPath, MoleculesName(1,:), '-Process', num2str(iP-1), '-PosteriorForBin', num2str(iBin));
         export_fig(FileName, '-pdf')
         close
      elseif SaveFigs == 2
         FolderPath = strcat(FigDirPath, '/Rates/');
         [status,msg,msgID] = mkdir(FolderPath);
         FileName   = strcat(FolderPath, MoleculesName(1,:), '-Process', num2str(iP-1), '-PosteriorForBin', num2str(iBin), '.fig');
         saveas(fig,FileName);
      end
      iFigure=iFigure+1;

    end

  end
  %%

end