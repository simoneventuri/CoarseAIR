%% The Function plots the Mole Fractions of the chemical system's component of interest vs its overall rates  
%
%  Input Arguments:  - iT:                      Index for the current Translational Temperature
%                    - t:                       Vector of time instants
%                    - ProcessesRatesOverall:   Matrix of Processes Overall Rates (Levels X 4(i.e.: Diss,Pair1,Pair2,Pair3) X NTint)
%                    - ...
%                    - Steps:                   Vector of Time Steps Indexes to Plot, used if MovieFlg=1
%                    - StpInstants:             Vector of Time Steps Indexes to Plot, used if MovieFlg=0
%
%  Input Global Var: - BinnedFlg:               If BinnedFlg=0, a Bar Plot is plotted (plot for levels' probabilities)
%                                               If BinnedFlg>1, a Histogram plot is plotted (plot for bins' probabilities, where each bin contains a number of levels equal to BinnedFlg)
%                    - MovieFlg                 Flag 0/1; if =0, a figure is plotted for each of StpInstants time instants; 
%                                               if =1, a movie is plotted for each of Steps time instants.
%

function [iFigure] = PlotProcProb(iT, iFigure, t, ProcessesRates, ProcessesRatesOverall, niRatio, EeV, DeltaEintDiss, vEeVVib, Steps, StpInstants, Leveljqn, Levelvqn)
    
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
  
  global BinnedFlg MovieFlg
  global NBins NComp NBinnedMol BinnedMolToComp BinnedMolName T0_Vec CompNames KinMthd QNSpace
  global SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz AxisLabelSz AxisLabelNm LegendFontNm XLimPlot YLimPlot PlotPairUservqnColor ...
         CompColor linS linST
  global RCVec BCVec GCVec KCVec OCVec PCVec WCVec JCVec YCVec CCVec MCVec

  ProcessName = ['Dissocition'; '    Elastic'; '  Inelastic'; '   Exchange'];
  iBinnedMol  = 1;

  if QNSpace == 1
    PDFMatrix = ones(max(Leveljqn(1:NBins(iBinnedMol),iBinnedMol))+1,max(Levelvqn(1:NBins(iBinnedMol),iBinnedMol))+1,4) .* 1.d-100;
    xxx     = [0:1:max(Levelvqn(1:NBins(iBinnedMol),iBinnedMol))];
    yyy     = [0:1:max(Leveljqn(1:NBins(iBinnedMol),iBinnedMol))];
    [X,Y] = meshgrid(xxx,yyy);
  end
  
  if MovieFlg == 1
    
    figure(iFigure)
    pause(20)
    StepsTemp = Steps;
    
  else
    
    StepsTemp = StpInstants;
    
  end
  
  for iSteps = StepsTemp
  
    if BinnedFlg > 1
      NLevBin         = BinnedFlg;
      EeVBinned       = zeros(ceil(NBins(iBinnedMol)/NLevBin),1);
      ProcessesBinned = zeros(ceil(NBins(iBinnedMol)/NLevBin),4);
      for iLevels = 1:NBins(iBinnedMol)
        EeVBinned(ceil(iLevels/NLevBin))             = EeVBinned(ceil(iLevels/NLevBin)) +            EeV(iLevels,iBinnedMol)           ./ NLevBin;
        %EeVBinned(ceil(iLevels/NLevBin))             = EeVBinned(ceil(iLevels/NLevBin)) +            DeltaEintDiss(iLevels,iBinnedMol) ./ NLevBin;
      end
    end 
      
    for iP = [1,2,3,4]
      
      ProcessesRatesPDF(1:NBins(iBinnedMol),iP,iT) = ProcessesRates(1:NBins(iBinnedMol),iP,iT) .* niRatio(iSteps,1:NBins(iBinnedMol),iBinnedMol)' ./ ProcessesRatesOverall(iSteps,iP,iT);
      Sum = 0.d0;
      for iBins = 1:size(ProcessesRates,1)
        Sum = Sum + ProcessesRatesPDF(iBins,iP,iT);
        ProcessesRatesCDF(iBins,iP,iT) = Sum;
      end
      
      if QNSpace == 1
        PDFMatrix = 1.d-100;
        for iLevels = 1:NBins(iBinnedMol)
          PDFMatrix(Leveljqn(iLevels,iBinnedMol)+1,Levelvqn(iLevels,iBinnedMol)+1,iP) = ProcessesRatesPDF(iLevels,iP,iT);
        end
      end
      
      if BinnedFlg > 1
        for iLevels = 1:NBins(iBinnedMol)
          ProcessesBinned(ceil(iLevels/NLevBin),iP,iT) = ProcessesBinned(ceil(iLevels/NLevBin),iP,iT) + ProcessesRatesPDF(iLevels,iP,iT) ./ NLevBin;
        end 
      end
      
    end
    
    if MovieFlg == 1
      
      subplot(1,3,[1:2])
      if sum(KinMthd(iBinnedMol,:) == 'CGM') == 3 || sum(KinMthd(iBinnedMol,:) == 'STS') == 3
        if BinnedFlg <= 1
          %h1=bar(EeV(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesPDF(1:NBins(iBinnedMol),1,iT),'FaceColor',RCVec,'EdgeColor',RCVec,'LineWidth',1.5);
          h1=bar(DeltaEintDiss(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesPDF(1:NBins(iBinnedMol),1,iT),'FaceColor',RCVec,'EdgeColor',RCVec,'LineWidth',1.5);
        else
          h1=bar(EeVBinned(:),ProcessesBinned(:,1),'FaceColor',RCVec,'EdgeColor',RCVec,'LineWidth',1.5);
        end
        h1.FaceAlpha = 0.3;
        h1.EdgeAlpha = 0.2;
        hold on
        if BinnedFlg <= 1
          %h2=bar(EeV(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesPDF(1:NBins(iBinnedMol),4,iT),'FaceColor',GCVec,'EdgeColor',GCVec,'LineWidth',1.5);
          h2=bar(DeltaEintDiss(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesPDF(1:NBins(iBinnedMol),4,iT),'FaceColor',GCVec,'EdgeColor',RCVec,'LineWidth',1.5);
        else
          h2=bar(EeVBinned(:),ProcessesBinned(:,4),'FaceColor',GCVec,'EdgeColor',GCVec,'LineWidth',1.5);
        end
        h2.FaceAlpha = 0.3;
        h2.EdgeAlpha = 0.2;
        xlabel(['Ro-Vibrational Energy [eV]']);
        xlim([min(EeV(1:NBins(iBinnedMol))), max(EeV(1:NBins(iBinnedMol)))]);
        %xlim([min(EeV(1:NBins(iBinnedMol))), 0.d0]);
        hold off
      elseif KinMthd(iBinnedMol,:) == 'VIB'
        bar(vEeVVib(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesPDF(1:NBins(iBinnedMol),1,iT),'FaceColor',RCVec,'EdgeColor',RCVec,'LineWidth',1.5);
        hold on
        bar(vEeVVib(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesPDF(1:NBins(iBinnedMol),4,iT),'FaceColor',GCVec,'EdgeColor',GCVec,'LineWidth',1.5);
        xlabel(['Vibrational Energy [eV]']);
      end
      ylabel(['PDF']);
      legend('Dissociation','Exchange');
      set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
      set(gcf, 'PaperPositionMode', 'auto');

      subplot(1,3,3)
      %subplot(2,3,[3,6])
      if sum(KinMthd(iBinnedMol,:) == 'CGM') == 3 || sum(KinMthd(iBinnedMol,:) == 'STS') == 3
        plot( EeV(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesCDF(1:NBins(iBinnedMol),1,iT),'-','Color',RCVec,'LineWidth',3 );
        hold on
        plot( EeV(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesCDF(1:NBins(iBinnedMol),4,iT),'-','Color',GCVec,'LineWidth',3 );
        xlabel(['Ro-Vibrational Energy [eV]']);
        xlim([min(EeV(1:NBins(iBinnedMol),iBinnedMol)), max(EeV(1:NBins(iBinnedMol),iBinnedMol))]);
        hold off
      elseif KinMthd(iBinnedMol,:) == 'VIB'
        plot( vEeVVib(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesCDF(1:NBins(iBinnedMol),1,iT),'-','Color',RCVec,'LineWidth',3 );
        hold on
        plot( vEeVVib(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesCDF(1:NBins(iBinnedMol),4,iT),'-','Color',GCVec,'LineWidth',3 );
        xlabel(['Vibrational Energy [eV]']);
        xlim([min(vEeVVib(1:NBins(iBinnedMol),iBinnedMol)), max(vEeVVib(1:NBins(iBinnedMol),iBinnedMol))]);
      end  
      ylabel(['CDF']);
      ylim([0, 1]);
      
      title([strcat('t = ',num2str(t(iSteps)),'s')]);
      set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
      set(gcf, 'PaperPositionMode', 'auto');
    
      pause(0.001)
    
    else
      
      figure(iFigure)
      fig = gcf;
      screensize = get( groot, 'Screensize' );
      fig.Position=screensize;
      fig.Color='None';

      left_color  = [0 0 0];
      right_color = [0 0 0];
      set(fig,'defaultAxesColorOrder',[left_color; right_color]);
      
      yyaxis left
      if sum(KinMthd(iBinnedMol,:) == 'CGM') == 3 || sum(KinMthd(iBinnedMol,:) == 'STS') == 3
        
        if BinnedFlg <= 1
          %h1=bar(EeV(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesPDF(1:NBins(iBinnedMol),1,iT),'FaceColor',RCVec,'EdgeColor',RCVec,'LineWidth',1.5);
          h1=bar(DeltaEintDiss(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesPDF(1:NBins(iBinnedMol),1,iT),'FaceColor',RCVec,'EdgeColor',RCVec,'LineWidth',1.5);
        else
          h1=bar(EeVBinned(:),ProcessesBinned(:,1,iT),'FaceColor',RCVec,'EdgeColor',RCVec,'LineWidth',1.5);
        end
        h1.FaceAlpha = 0.4;
        h1.EdgeAlpha = 0.2;
        hold on
        
        if BinnedFlg <= 1
          %h2=bar(EeV(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesPDF(1:NBins(iBinnedMol),4,iT),'FaceColor',GCVec,'EdgeColor',GCVec,'LineWidth',1.5);
          h2=bar(DeltaEintDiss(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesPDF(1:NBins(iBinnedMol),4,iT),'FaceColor',GCVec,'EdgeColor',GCVec,'LineWidth',1.5);   
        else
          h2=bar(EeVBinned(:),ProcessesBinned(:,4,iT),'FaceColor',GCVec,'EdgeColor',GCVec,'LineWidth',1.5);
        end       
        h2.FaceAlpha = 0.4;
        h2.EdgeAlpha = 0.2;
        
        str_x = ['Energy [eV]'];
        %xlim([min(EeV(1:NBins(iBinnedMol))), max(EeV(1:NBins(iBinnedMol)))]);
        %xlim([min(EeV(1:NBins(iBinnedMol))), 0.d0]);
        hold off
      elseif KinMthd(iBinnedMol,:) == 'VIB'
        bar(vEeVVib(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesPDF(1:NBins(iBinnedMol),1,iT),'FaceColor',RCVec,'EdgeColor',RCVec,'LineWidth',1.5);
        hold on
        bar(vEeVVib(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesPDF(1:NBins(iBinnedMol),4,iT),'FaceColor',GCVec,'EdgeColor',GCVec,'LineWidth',1.5);
        str_x = ['Vibrational Energy [eV]'];
      end    
      
      str_y = ['PDF'];
      ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
      ylab.Interpreter = 'latex';
      
      yt = get(gca, 'YTick');
      set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
      
      
      yyaxis right
      plot( EeV(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesCDF(1:NBins(iBinnedMol),1,iT),'-','Color',RCVec,'LineWidth',3 );
      %plot( DeltaEintDiss(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesCDF(1:NBins(iBinnedMol),1,iT),'-','Color',RCVec,'LineWidth',3 );
      hold on
      plot( EeV(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesCDF(1:NBins(iBinnedMol),4,iT),'-','Color',GCVec,'LineWidth',3 );
      %plot( DeltaEintDiss(1:NBins(iBinnedMol),iBinnedMol), ProcessesRatesCDF(1:NBins(iBinnedMol),4,iT),'-','Color',GCVec,'LineWidth',3 );

      str_y = ['PDF'];
      ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
      ylab.Interpreter = 'latex';
      
      xt = get(gca, 'XTick');
      set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
      yt = get(gca, 'YTick');
      set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

      xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
      xlab.Interpreter = 'latex';

      str_y = ['CDF'];
      ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
      ylab.Interpreter = 'latex';
      
      TimeStr = strcat('Dissociation, Time = ',num2str(t(iSteps),'%7.2e\n'),' s');
      LegendText = {TimeStr; 'Exchange'};
      %LegendText = [TimeStr];
      clab = legend(LegendText,'Location','Best');
      clab.Interpreter = 'latex';
      set(clab,'FontSize',LegendFontSz,'FontName',LegendFontNm,'Interpreter','latex');

      if SaveFigs == 1
        FolderPath = strcat(FigDirPath, '/FlowQuantities/');
        [status,msg,msgID] = mkdir(FolderPath);
        FileName   = strcat(FolderPath, MoleculesName(1,:),'-Pops');
        export_fig(FileName, '-pdf')
        close
      elseif SaveFigs == 2
        FolderPath = strcat(FigDirPath, '/FlowQuantities/');
        [status,msg,msgID] = mkdir(FolderPath);
        FileName   = strcat(FolderPath, MoleculesName(1,:),'-Pops.fig');
        savefig(FileName)
        close
      end

      iFigure=iFigure+1;   
      
    end

  end
  
  if MovieFlg == 1
    iFigure=iFigure+1;   
  end
    
end