%% The Function plots the Overall Rates of the chemical system's component of interest
%
%  Input Arguments:  - iT:                      Index for the current Translational Temperature
%                    - t:                       Vector of time instants
%                    - ProcessesRatesOverall:   Matrix of Processes Overall Rates (Levels X 4(i.e.: Diss,Pair1,Pair2,Pair3) X NTint)
%
%  Input Global Var: - CompOI:                  Component of interest
%                    - XLimPlot:                Vector of [Min, Max] for x-axes plot
%                    - YLimPlot:                Vector of [Min, Max] for y-axes plot
%                    - PlotPairUser:            Vector of Flags 0/1; if =1 the atomic pair is plotted
%

function [iFigure, KDiss, PDF_KDiss, KDissQSS, PDF_KDissQSS, KExch, PDF_KExch, KExchQSS, PDF_KExchQSS] = PlotKsStoch(iT, iFigure, t, ProcessesRatesOverallStoch)    
    
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

  global PlotPairUser XLimPlot YLimPlot 
  global linST SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz AxisLabelSz AxisLabelNm LegendFontNm 
  global iPESStart iPESEnd PlotSamplesFlg PlotMomentsFlg MergePairsFlg tQSS

  iBinnedMol=1;
  
  NPESTemp = iPESEnd - iPESStart + 1; 

  
  
  %ProcessesRatesOverallStochSigmaTemp = ProcessesRatesOverallStochSigma.^2;

  if MergePairsFlg == 1
    for iPES=1:NPESTemp
      ProcessesRatesOverallStochTemp(:,1,iPES) = ProcessesRatesOverallStoch(:,1,iPES+iPESStart-1);
      ProcessesRatesOverallStochTemp(:,3,iPES) = ProcessesRatesOverallStoch(:,3,iPES+iPESStart-1) + ProcessesRatesOverallStoch(:,4,iPES+iPESStart-1);
%       for iP = 3:-1:2
%         jP = 1;
%         while jP < iP
%           if Pair_To_Molecule(iP) == Pair_To_Molecule(jP)
%             iP;
%             jP;
%             ProcessesRatesOverallStochTemp(:,jP+1,iPES)      = ProcessesRatesOverallStochTemp(:,jP+1,iPES)         + ProcessesRatesOverallStochTemp(:,iP+1,iPES);
%             %ProcessesRatesOverallStochSigmaTemp(:,jP+1,iPES) = ProcessesRatesOverallStochSigmaTemp(:,jP+1,iPES)    + ProcessesRatesOverallStochSigmaTemp(:,iP+1,iPES);
%             jP = iP;
%           end
%           jP = jP+1;
%         end
%       end
%       %ProcessesRatesOverallStochSigmaTemp = ProcessesRatesOverallStochSigmaTemp.^0.5d0;
    end
  end
  
  tStart   = t(1);
  tEnd     = t(end);
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

  itQSS=1;
  while (t(itQSS) < tQSS)
    itQSS=itQSS+1;
  end
  

  if PlotSamplesFlg == 1 || PlotMomentsFlg == 1        
    
    figure(iFigure)
    fig = gcf;
    screensize = get( groot, 'Screensize' );
    fig.Position=screensize;
    fig.Color='None';
    LegendText = [];
    
    Summ1    = zeros(NPESTemp,1);
    SummSqr1 = zeros(NPESTemp,1);
    Summ3    = zeros(NPESTemp,3);
    SummSqr3 = zeros(NPESTemp,3);
    
    for iPES=1:NPESTemp
        iPES
      %for iComp=1:2%NComp
        
        if PlotMomentsFlg == 1
          Summ1(1:itEnd-itStart+1)    = Summ1(1:itEnd-itStart+1)    + ProcessesRatesOverallStochTemp(1:itEnd-itStart+1,1,iPES);
          SummSqr1(1:itEnd-itStart+1) = SummSqr1(1:itEnd-itStart+1) + ProcessesRatesOverallStochTemp(1:itEnd-itStart+1,1,iPES).^2;
          Summ3(1:itEnd-itStart+1)    = Summ3(1:itEnd-itStart+1)    + ProcessesRatesOverallStochTemp(1:itEnd-itStart+1,3,iPES);
          SummSqr3(1:itEnd-itStart+1) = SummSqr3(1:itEnd-itStart+1) + ProcessesRatesOverallStochTemp(1:itEnd-itStart+1,3,iPES).^2;
        end
        
        KDissQSS(iPES) = ProcessesRatesOverallStochTemp(itQSS,1,iPES);
        KExchQSS(iPES) = ProcessesRatesOverallStochTemp(itQSS,3,iPES);
        
        KDiss(iPES) = ProcessesRatesOverallStochTemp(end,1,iPES);
        KExch(iPES) = ProcessesRatesOverallStochTemp(end,3,iPES);

        if PlotSamplesFlg == 1
          loglog(t(:),ProcessesRatesOverallStochTemp(:,1,iPES),'k','LineWidth',3);
          hold on  
          loglog(t(:),ProcessesRatesOverallStochTemp(:,3,iPES),'r:','LineWidth',3);
        end

      %end
    end

    if PlotMomentsFlg == 1
      Meann1  = Summ1 ./ NPESTemp;
      StDevv1 = sqrt(  SummSqr1./NPESTemp - Meann1.^2);
      Meann3  = Summ3 ./ NPESTemp;
      StDevv3 = sqrt(  SummSqr3./NPESTemp - Meann3.^2);

      semilogy(t(itStart:itEnd),Meann1(:),'ro');
      hold on
      semilogy([t(itStart:itEnd)'; t(itStart:itEnd)'],[max(Meann1'-3.0.*StDevv1',Meann1'.*0.0+1.d-100); Meann1'+3.0.*StDevv1'],'r-');
      
      semilogy(t(itStart:itEnd),Meann3(:),'ro');
      hold on
      semilogy([t(itStart:itEnd)'; t(itStart:itEnd)'],[max(Meann3'-3.0.*StDevv3',Meann3'.*0.0+1.d-100); Meann3'+3.0.*StDevv3'],'r-');
    end

    xt = get(gca, 'XTick');
    set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

    str_x = ['Time [s]'];
    xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
    xlab.Interpreter = 'latex';
    xlim([tStart, tEnd]);

    str_y = ['$\bar{K}_{Process}$'];
    ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
    ylab.Interpreter = 'latex';
    %ylim(YLimPlot);

    if SaveFigs == 1
    FolderPath = strcat(FigDirPath, '/FlowQuantities/');
    [status,msg,msgID] = mkdir(FolderPath);
    FileName   = strcat(FolderPath, 'Ks');
    export_fig(FileName, '-pdf')
    close
    elseif SaveFigs == 2
    FolderPath = strcat(FigDirPath, '/FlowQuantities/');
    [status,msg,msgID] = mkdir(FolderPath);
    FileName   = strcat(FolderPath, 'Ks.fig');
    savefig(FileName)
    close
    end

    iFigure=iFigure+1;   
      
  end
  
  
  
  fig1=figure(iFigure);
  fig = gcf;
  if SaveFigs > 0
    screensize = get( groot, 'Screensize' );
    fig.Position=screensize;
    fig.Color='None';
  end 

  h         = histfit(KDiss,1000,'lognormal');
  hold on 
  PDF_KDiss = fitdist(KDiss','lognormal');
  %pd = makedist('Lognormal','mu',PDF_KDiss.mu,'sigma',PDF_KDiss.sigma);
  %x = logspace(-13,-11,1000);
  %y = pdf(pd,x);
  %semilogx(x,y);
  %r = lognrnd(PDF_KDiss.mu,PDF_KDiss.sigma,[1,50])
  %histogram(r,30)
  h            = histfit(KDissQSS,30,'lognormal');
  PDF_KDissQSS = fitdist(KDissQSS','lognormal');
  h(1).FaceAlpha=0.4;
  set(gca, 'XScale', 'log')
  xlab = xlabel('$K^{Diss} [cm^3/s]$');
  xlab.Interpreter = 'latex';
  ylab = ylabel('Posterior Samples');
  grid on

  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  %daspect([1 1 1])
  iFigure = iFigure+1;
  
  
  fig1=figure(iFigure);
  fig = gcf;
  if SaveFigs > 0
    screensize = get( groot, 'Screensize' );
    fig.Position=screensize;
    fig.Color='None';
  end 
  h            = histfit(KExch,30,'lognormal');
  hold on
  PDF_KExch    = fitdist(KExch','lognormal');

  h            = histfit(KExchQSS,30,'lognormal');
  PDF_KExchQSS = fitdist(KExchQSS','lognormal');
  h(1).FaceAlpha=0.4;
  set(gca, 'XScale', 'log')
  xlab = xlabel('$K^{Exch} [cm^3/s]$');
  xlab.Interpreter = 'latex';
  ylab = ylabel('Posterior Samples');
  grid on
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  %daspect([1 1 1])
  iFigure = iFigure+1;
  
end