%% The Function plots the Levels' Populations at different time instants  
%
%  Input Arguments:  - iT:                      Index for the current Translational Temperature
%                    - t:                       Vector of time instants
%                    - ProcessesRatesOverall:   Matrix of Processes Overall Rates (Levels X 4(i.e.: Diss,Pair1,Pair2,Pair3) X NTint)
%                    - MolFracs:                Matrix of Mole Fractions (Time-Instants X Components)
%                    ...
%                    - LevToBin:                (Optional) Mapping Levels-Bins
%
%  Input Global Var: - StepsOverlappingSteps:   Flag 0/1; if =1, different time steps are plotted in the same figure
%                    - MinEvPlot:               Min for x-axes plot
%                    - MaxEvPlot:               Max for x-axes plot
%                    - LevToBinFlg:             Flag 0/1; if =1, populations are colored based on the levels' bins (LevToBin is used)
%                    - vqnColor:                Flag 0/1; if =1, populations are colored based on the levels' vqn (Levelvqn is used)
%

function [iFigure] = PlotSteps(iT, iFigure, t, MolFracs, ProcessesRatesOverall, NLevels, LevelEeV, LevelEeV0, LevelQ, Levelg, LevelToBin, Pop, QBins, StpInstants, Levelvqn, Leveljqn, LevToBin, DeltaEintDiss, ProcessesRates, KRec, KRecOverg)    

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
  
  global MinEvPlot MaxEvPlot StepsOverlappingSteps LevToBinFlg vqnColor KinMthd
  global MoleculesName
  global SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz AxisLabelSz AxisLabelNm LegendFontNm XLimPlot YLimPlot PlotPairUservqnColor 
  global Ue T0_Vec UKb
  
  iMol = 1

  if StepsOverlappingSteps == 1
    figure(iFigure)
    fig = gcf;
    screensize = get( groot, 'Screensize' );
    fig.Position=screensize;
    fig.Color='None';
    LegendText = [];
  end 

  ii = 0;
  for iSteps=StpInstants(1:end-1)
    ii = ii+1;
    
    %cmap     = colormap(lines(max(LevToBin(:,iMol,ii))));

    if StepsOverlappingSteps == 0
      figure(iFigure)
      fig = gcf;
      screensize = get( groot, 'Screensize' );
      fig.Position=screensize;
      fig.Color='None';
      title(['Time = ', num2str(t(iSteps)), ' s'])
      hold on
    end

    if sum(KinMthd(1,:) == 'CGM') == 3
      LevelPop(1:NLevels(iMol),iMol) = Pop(iSteps,LevelToBin(1:NLevels(iMol),iMol),iMol)' ./ QBins(LevelToBin(1:NLevels(iMol),iMol),iMol) .* exp( - LevelEeV0(1:NLevels(iMol),iMol) .* Ue ./ (T0_Vec(1) .* UKb) );
      LegendPlot(iSteps) = LevelPop(1,iMol);
    else
      LevelPop(1:NLevels(iMol),iMol) = Pop(iSteps,LevelToBin(1:NLevels(iMol),iMol),iMol)' ./ QBins(LevelToBin(1:NLevels(iMol),iMol),iMol);
      LegendPlot(iSteps) = LevelPop(1,iMol);
    end 

    if vqnColor == 0

      if LevToBinFlg == 0

        scatter(LevelEeV(1:NLevels(iMol),iMol),LevelPop(1:NLevels(iMol),iMol),20,'Filled');
        %scatter(LevelEeV(1:NLevels(iMol),iMol),LevelPop(1:NLevels(iMol),iMol),20,log10(ProcessesRates(1:NLevels(iMol),1,1)),'Filled');
        %caxis([-10 -8]);
        %scatter(LevelEeV(1:NLevels(iMol),iMol),LevelPop(1:NLevels(iMol),iMol),20,log10(KRecOverg(1:NLevels(iMol),1,1)),'Filled');
        %caxis([-11 -9]);
        %colormap(lines);
        %h=colorbar;
        %ylab = ylabel(h, '$log_{10}(K_{Diss_{i}}) [cm^3/(mol*s)]$');
        
      else

        %scatter3(LevelEeV(1:NLevels(iMol),iMol),LevToBin(1:NLevels(iMol),iMol),LevelPop(1:NLevels(iMol),iMol),20,LevToBin(1:NLevels(iMol),iMol,ii),'Filled');
        scatter(LevelEeV(1:NLevels(iMol),iMol),LevelPop(1:NLevels(iMol),iMol),20,LevToBin(1:NLevels(iMol),iMol,ii),'Filled');
        %scatter(LevelEeV(1:NLevels(iMol),iMol),LevelPop(1:NLevels(iMol),iMol),20,'Filled');

%         colorMap = [140 0, 0; 20, 120, 45]./256;
%         colormap(colorMap);
%  
        %colormap(lines(max(max(max(LevToBin)))))
        colormap(distinguishable_colors(max(max(max(LevToBin)))))
        cb=colorbar;
        %cb.Ticks = [1, 1.5]; %Create 8 ticks from zero to 1
        %cb.TickLabels = {'1','2'}
        ylab = ylabel(cb, '$Group$');
        ylab.Interpreter = 'latex';
        set(cb,'FontSize',LegendFontSz,'FontName',LegendFontNm,'TickLabelInterpreter','latex');

      end

    else
      
      ColorMat = distinguishable_colors(max(max(Levelvqn))+1);
      
      for iv = 0:max(Levelvqn(1:NLevels(iMol),iMol));
        jj = 0;
        for iLevel=1:NLevels(iMol)
          if Levelvqn(iLevel,iMol) == iv
            jj = jj + 1;
            LevelEeVTemp(jj) = LevelEeV(iLevel,iMol);
            LevelPopTemp(jj) = LevelPop(iLevel,iMol);
          end
        end
        semilogy(LevelEeVTemp',LevelPopTemp','Color',ColorMat(iv+1,:),'LineWidth',3)
        hold on
        clear LevelEeVTemp LevelPopTemp
      end
        
      
      %scatter(LevelEeV(1:NLevels(iMol),iMol),LevelPop(1:NLevels(iMol),iMol),10,DeltaEintDiss(1:NLevels(iMol),iMol),'Filled');
      %scatter(LevelEeV(1:NLevels(iMol),iMol),LevelPop(1:NLevels(iMol),iMol),50,Levelvqn(1:NLevels(iMol),iMol),'Filled');
      %scatter3(LevelEeV(1:NLevels(iMol),iMol),Levelvqn(1:NLevels(iMol),iMol),LevelPop(1:NLevels(iMol),iMol),50,Levelvqn(1:NLevels(iMol),iMol),'Filled');
      colormap(distinguishable_colors(max(max(Levelvqn))+1))
      %colormap(lines(max(max(Levelvqn)))+1);
      cb=colorbar;
      %ylab = ylabel(cb, '$E_i - E_i^{Diss}$');
      ylab = ylabel(cb, 'v');
      ylab.Interpreter = 'latex';
      set(cb,'FontSize',LegendFontSz,'FontName',LegendFontNm,'TickLabelInterpreter','latex');
      %ylim(cb,[-8, 0]); caxis([-8 0]);
      %ylim(cb,[-6.7, 0]); caxis([-6.7 0]);
      %ylim(cb,[-4.3, 0]); caxis([-4.3 0]);
      %ylim(cb,[-3.6, 0]); caxis([-3.6 0]);
      %ylim(cb,[-3, 0]); caxis([-3 0]);
      %ylim(cb,[-2, 0]); caxis([-2 0]);
      %ylim(cb,[-1.1, 0]); caxis([-1.1 0]);
      %ylim(cb,[-0.8, 0]); caxis([-0.8 0]);
      %ylim(cb,[-0.8, 0]); caxis([-0.8 0]);
      %colormap(lines)
      %colormap(jet)
      
    end
    hold on

    xt = get(gca, 'XTick');
    set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');

%     str_x = ['Energy [eV]'];
%     xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
%     xlab.Interpreter = 'latex';
%     xlim([max(min(LevelEeV(:,iMol)),MinEvPlot(iMol)), min(max(LevelEeV(:,iMol)),MaxEvPlot(iMol))]);

    str_y = ['$N_{i} / g_{i} [m^{-3}]$'];
    ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
    ylab.Interpreter = 'latex';
    ylim([1.d5, 1.d23]);
    set(gca, 'YScale', 'log')
    
%     str_z = ['$N_{i} / g_{i} [m^{-3}]$'];
%     zlab = zlabel(str_z,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
%     zlab.Interpreter = 'latex';
%     zlim([1.d0, 1.d20]);
%     set(gca, 'ZScale', 'log')

    if SaveFigs == 1
      FolderPath = strcat(FigDirPath, '/FlowQuantities/');
      [status,msg,msgID] = mkdir(FolderPath);
      FileName   = strcat(FolderPath, MoleculesName(1,:),'-Pop_',num2str(iSteps),'_t',num2str(t(iSteps),'%7.2e'));
      export_fig(FileName, '-pdf')
      close
    elseif SaveFigs == 2
      FolderPath = strcat(FigDirPath, '/FlowQuantities/');
      [status,msg,msgID] = mkdir(FolderPath);
      FileName   = strcat(FolderPath, MoleculesName(1,:),'-Pop_',num2str(iSteps),'_t',num2str(t(iSteps),'%7.2e'),'.fig');
      savefig(FileName)
      close
    end

    if StepsOverlappingSteps == 0
      iFigure=iFigure+1;
    else
      LegendText = [LegendText; strcat('Time = ',num2str(t(iSteps),'%7.2e\n'),' s')]
    end

  end

  if StepsOverlappingSteps == 1
    
    scatter(0.d0.*LegendPlot+LevelEeV(1,iMol),LegendPlot,20,'Filled');
    
    clab = legend(LegendText,'Location','Best');
    clab.Interpreter = 'latex';
    set(clab,'FontSize',LegendFontSz,'FontName',LegendFontNm,'Interpreter','latex');

    if SaveFigs == 1
      FolderPath = strcat(FigDirPath, '/FlowQuantities/');
      [status,msg,msgID] = mkdir(FolderPath);
      FileName   = strcat(FolderPath, MoleculesName(1,:),'-Pop_',num2str(iSteps),'_t',num2str(t(iSteps),'%7.2e'));
      export_fig(FileName, '-pdf')
      close
    elseif SaveFigs == 2
      FolderPath = strcat(FigDirPath, '/FlowQuantities/');
      [status,msg,msgID] = mkdir(FolderPath);
      FileName   = strcat(FolderPath, MoleculesName(1,:),'-Pop_',num2str(iSteps),'_t',num2str(t(iSteps),'%7.2e'),'.fig');
      savefig(FileName)
      close
    end
    
    iFigure=iFigure+1;
    
  end

end