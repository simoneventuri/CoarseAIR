function PlotFigure(iFigure, X, Y, Exts_x, Exts_y, str_x, str_y, LineStr, LineColor, LineSize, LegendFlg, LegendStr, FolderStr, FileStr )

  global AxisFontSz AxisFontNm AxisLabelSz AxisLabelNm LegendFontSz LegendFontNm SaveFigs FigDirPath 
  
  figure(iFigure)
  fig = gcf;
  screensize = get( groot, 'Screensize' );
  fig.Position=screensize;
  fig.Color='None';
  plot(X,Y,LineStr, 'Color',LineColor, 'LineWidth',LineSize )
  hold on
  
  xlim(Exts_x);
  ylim(Exts_y);
  
  xt = get(gca, 'XTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  yt = get(gca, 'YTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
   
  xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  xlab.Interpreter = 'latex';
     
  ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  ylab.Interpreter = 'latex';

  if LegendFlg == 1
    clab = legend(LegendStr,'Location','Best');
    clab.Interpreter = 'latex';
    set(clab,'FontSize',LegendFontSz,'FontName',LegendFontNm,'Interpreter','latex');
  end

  if SaveFigs == 1
    FolderPath = strcat(FigDirPath, FolderStr);
    [status,msg,msgID] = mkdir(FolderPath);
    FileName   = strcat(FolderPath, FileStr);
    export_fig(FileName, '-pdf')
    close
  elseif SaveFigs == 2
    FolderPath = strcat(FigDirPath, FolderStr);
    [status,msg,msgID] = mkdir(FolderPath);
    FileName   = strcat(FolderPath, FileStr, '.fig');
    savefig(FileName)
    close
  end
        
end
