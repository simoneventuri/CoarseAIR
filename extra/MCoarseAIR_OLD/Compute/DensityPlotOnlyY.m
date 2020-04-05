function [iFigure, yHist] = DensityPlotOnlyY(iFigure, x, Y, nBinsy, LogScaleFlg)
  
  [xSorted,xSortIdx] = sort(x);
  
  if LogScaleFlg     == 1
    ySpace = logspace(min(min(log10(Y)))*1.01,max(max(log10(Y)))*0.99,nBinsy+1)';
  elseif LogScaleFlg == 0
    ySpace = linspace(min(min(Y))*0.8,max(max(Y))*1.2,nBinsy+1)';
  elseif LogScaleFlg == -1
    ySpace = linspace(min(min(Y)),max(max(Y)),nBinsy+1)';
  end
  
  yHist = [];
  xTemp = [];
  yTemp = [];
  for i=1:length(xSorted)
    X  = xSorted(i);
    if (i==1)
      XM = X;
    else
      XM = xSorted(i-1);
    end
    if (i==length(xSorted))
      XP = X;
    else
      XP = xSorted(i+1);
    end
    xTemp = [xTemp, ySpace.*0.0 + (XM + X)/2, ySpace.*0.0 + (X + XP)/2];
    
    yTemp = [yTemp, ySpace, ySpace];
   
    y     = Y(xSortIdx(i),:);
    Histt = histogram(y,ySpace')';
    yHist = [yHist, [Histt.Values', Histt.Values';0,0] ];
    
  end
  
  H = pcolor(xTemp, yTemp, yHist);
  box on
  hold on
  shading interp
  set(H,'edgecolor','none');
  if LogScaleFlg     == 1
    set(gca, 'YScale', 'log');
  end
  colorbar
  colormap jet
  
end