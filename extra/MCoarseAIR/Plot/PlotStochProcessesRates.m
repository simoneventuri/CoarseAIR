%% The Function plots the processes rates in the Energy Space: K(iLevel,iProcess) = Mol1(iLevel) -> Mol2(iProcess)
%
%  Input Arguments:  - iFigure              Index for the first figure to plot
%                    - ProcessesRatesStoch       Matrix of Processes Rates
%                    - ProcessesRatesStochSigma  Matrix of Processes Rates St. Deviations
%                    - LevelEeV             Levels energy in eV, with the 0 corresponding to the dissociation energy
%                    - EeV                  Levels/Bins' energy in eV, with the 0 corresponding to the dissociation energy
%                    - LevelToBin           (Optional) Mapping Levels To Bins 
%   
%  Input Global Var: - MergePairsFlg        Flag 0/1; if 1, all the rates to same molecules but different pairs, are merged together
%                    - PlotPairUser         Vector of flags for the pairs to plot (e.g.: [1 0 0])
%                    - SigmaOn              Flag 0/1; if 1, error bars are plotted.
%                    - OverlapFlg           Flag 0/1; if 1, overlapping plots in the same figure
%                    - LevToBinFlg          Flag 0/1; if 1, rates are colored based on the mapping levels-to-bins
%                    - vqnColor             Flag 0/1; if 1, rates are colored based on their vqn
%

%function iFigure = PlotProcessesRatesStoch(iFigure, ProcessesRatesStoch, ProcessesRatesStochSigma, LevelEeV, EeV, LevToBin)   
function iFigure = PlotStochProcessesRates(iFigure, NLevels, ProcessesRatesStoch, ProcessesRatesStochMmnts, Levelg, LevelEeV, EeV, DeltaEintDiss, Levelvqn, Leveljqn, LevToBin, rOut, rIn, Egam, Tau, LevelQ, dVIn, ddVIn, dVOut, ddVOut, dVdJIn, dVdJOut, QBins)   

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
  
  global MergePairsFlg PlotPairUser SigmaOn OverlapFlg LevToBinFlg vqnColor
  global NPESs NBins ColorVec linS MoleculesName Pair_To_Molecule 
  global SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz LegendFontNm
  global Plnck UKb Ue KeV AvN AMUToKg EhToeV DSWtoKg ATMToPa
  global T0_Vec PlotPColorFlg PlotSamplesFlg PlotMomentsFlg PlotPostCut BinPostCutVec
  global iPESEnd iPESStart

  NPESTemp = iPESEnd - iPESStart + 1; 
  
  
  ProcessesRatesStochTemp      = ProcessesRatesStoch;
  %ProcessesRatesStochSigmaTemp = ProcessesRatesStochSigma.^2;

  %if MergePairsFlg == 1
    for iPES=1:NPESTemp
      ProcessesRatesStochTemp(:,3,iPES)      = ProcessesRatesStochTemp(:,3,iPES+iPESStart-1) + ProcessesRatesStochTemp(:,4,iPES+iPESStart-1);
%       for iP = 3:-1:2
%         jP = 1;
%         while jP < iP
%           if Pair_To_Molecule(iP) == Pair_To_Molecule(jP)
%             iP;
%             jP;
%             ProcessesRatesStochTemp(:,jP+1,iPES)      = ProcessesRatesStochTemp(:,jP+1,iPES)         + ProcessesRatesStochTemp(:,iP+1,iPES);
%             %ProcessesRatesStochSigmaTemp(:,jP+1,iPES) = ProcessesRatesStochSigmaTemp(:,jP+1,iPES)    + ProcessesRatesStochSigmaTemp(:,iP+1,iPES);
%             jP = iP;
%           end
%           jP = jP+1;
%         end
%       end
%       %ProcessesRatesStochSigmaTemp = ProcessesRatesStochSigmaTemp.^0.5d0;
    end
  %end
  
  
  jj = 1;
  for i = 1:size(LevelEeV,1)
%     if ddVOut(i,1) < 0.d0 
      BoundLevels(jj) =  i;
      jj = jj + 1;
%     end
  end

  
  for iP=[1,3,4]
    if (PlotPairUser(iP) == 1) 

      
      %%
      if PlotPColorFlg == 1

        x = EeV(1:NBins(1),1);
        Y = squeeze(ProcessesRatesStochTemp(1:NBins(1),iP,1:NPESTemp));

        iP
        fig1=figure(iFigure);
        fig = gcf;
        if SaveFigs > 0
          screensize = get( groot, 'Screensize' );
          fig.Position=screensize;
          fig.Color='None';
        end 

        [iFigure, yHist] = DensityPlotOnlyY(iFigure, x, Y, 300, 1);

        set(gca, 'YScale', 'log');
        xlab = xlabel('$E_i [eV]$');
        xlab.Interpreter = 'latex';
        if iP==1
          ylab = ylabel('$K_{i}^{Diss} [cm^3/s]$');
        else
          ylab = ylabel('$K_{i}^{Exch} [cm^3/s]$');
        end
        ylab.Interpreter = 'latex';
        grid on

        set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
        %daspect([1 1 1])

        if SaveFigs == 1
           FolderPath = strcat(FigDirPath, '/Rates/');
           [status,msg,msgID] = mkdir(FolderPath);
           FileName   = strcat(FolderPath, MoleculesName(1,:),'-Process', num2str(iP-1), 'PosteriorPColor');
           export_fig(FileName, '-pdf')
           close
        elseif SaveFigs == 2
           FolderPath = strcat(FigDirPath, '/Rates/');
           [status,msg,msgID] = mkdir(FolderPath);
           FileName   = strcat(FolderPath, MoleculesName(1,:),'-Process', num2str(iP-1), 'PosteriorPColor.fig');
           saveas(fig,FileName);
        end
        iFigure=iFigure+1;

      end
      %%
      
      
      %%
      if PlotSamplesFlg == 1 || PlotMomentsFlg == 1

        iP
        fig1=figure(iFigure);
        fig = gcf;
        if SaveFigs > 0
          screensize = get( groot, 'Screensize' );
          fig.Position=screensize;
          fig.Color='None';
        end 
        
        Summ    = zeros(NBins(1),1);
        SummSqr = zeros(NBins(1),1);
        
        for iPES=1:NPESTemp
          
          if PlotMomentsFlg == 1
            Summ(1:NBins(1))    = Summ(1:NBins(1))    + ProcessesRatesStochTemp(1:NBins(1),iP,iPES);
            SummSqr(1:NBins(1)) = SummSqr(1:NBins(1)) + ProcessesRatesStochTemp(1:NBins(1),iP,iPES).^2;
          end
          
          if PlotSamplesFlg == 1
            semilogy(EeV(1:NBins(1),1),ProcessesRatesStochTemp(1:NBins(1),iP,iPES),'ko','Marker','o','Markersize',3,'MarkerFaceColor','k');
          end
          hold on

        end
        
        if PlotMomentsFlg == 1
          Meann  = Summ ./ NPESTemp;
          StDevv = sqrt(  SummSqr./NPESTemp - Meann.^2);

          semilogy(EeV(1:NBins(1),1),Meann,'ro');
          hold on
          semilogy([EeV(1:NBins(1),1)'; EeV(1:NBins(1),1)'],[max(Meann'-3.0.*StDevv',Meann'.*0.0+1.d-100); Meann'+3.0.*StDevv'],'r-');
        end
        
        set(gca, 'YScale', 'log');
        xlab = xlabel('$E_i [eV]$');
        xlab.Interpreter = 'latex';
        if iP==1
          ylab = ylabel('$K_{i}^{Diss} [cm^3/s]$');
          ylim([1.d-15, 1.d-8])
        else
          ylab = ylabel('$K_{i}^{Exch} [cm^3/s]$');
        end
        ylab.Interpreter = 'latex';
        grid on

        set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
        %daspect([1 1 1])

        if SaveFigs == 1
           FolderPath = strcat(FigDirPath, '/Rates/');
           [status,msg,msgID] = mkdir(FolderPath);
           FileName   = strcat(FolderPath, MoleculesName(1,:),'-Process', num2str(iP-1), 'Posterior');
           export_fig(FileName, '-pdf')
           close
        elseif SaveFigs == 2
           FolderPath = strcat(FigDirPath, '/Rates/');
           [status,msg,msgID] = mkdir(FolderPath);
           FileName   = strcat(FolderPath, MoleculesName(1,:),'-OverallToPair', num2str(iP-1), 'Posterior.fig');
           saveas(fig,FileName);
        end
        iFigure=iFigure+1;

      end
      %%
     
      
      %% 
      if PlotPostCut == 1
        
        for iBin = BinPostCutVec

          fig1=figure(iFigure);
          fig = gcf;
          if SaveFigs > 0
            screensize = get( groot, 'Screensize' );
            fig.Position=screensize;
            fig.Color='None';
          end 
   
          RatesVec = squeeze(ProcessesRatesStochTemp(iBin,iP,:));
          histfit(RatesVec,30)
           
          if iP==1
            xlab = xlabel('$K_{i}^{Diss} [cm^3/s]$');
          else
            xlab = xlabel('$K_{i}^{Exch} [cm^3/s]$');
          end
          xlab.Interpreter = 'latex';
          ylab = ylabel('Posterior Samples');
          grid on
          title(strcat('Bin ', num2str(iBin)))

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
  end

  
  iBinnedMol = 1;
  ExpVec(:,1)  = QBins(:,1);
  ExpVec       = ExpVec ./ sum(ExpVec);
  
  
%   ThermalRatesStoch = zeros(NPESs,4);
%   for iPES=iPESStart:iPESEnd
%     for iP = 1:4
%       ThermalRatesStoch(iPES,iP) = sum(ProcessesRatesStoch(:,iP,iPES) .* ExpVec(:,1));
%     end
%   end
%   
%   for iP = 1:4
%     if (iP == 1)
%       fprintf('  Dissociation Thermal Rate Mean = %e cm^s/(mol*s); St. Deviation = %e\n', mean(ThermalRatesStoch(:,iP)), std(ThermalRatesStoch(:,iP)) ) 
%     elseif (iP == 2)
%       fprintf('  Inelastic    Thermal Rate Mean = %e cm^s/(mol*s); St. Deviation = %e\n', mean(ThermalRatesStoch(:,iP)), std(ThermalRatesStoch(:,iP)) )
%     elseif (iP == 3)
%       fprintf('  Exchange 1   Thermal Rate Mean = %e cm^s/(mol*s); St. Deviation = %e\n', mean(ThermalRatesStoch(:,iP)), std(ThermalRatesStoch(:,iP)) )
%     elseif (iP == 4)
%       fprintf('  Exchange 2   Thermal Rate Mean = %e cm^s/(mol*s); St. Deviation = %e\n', mean(ThermalRatesStoch(:,iP)), std(ThermalRatesStoch(:,iP)) )
%     end 
%   end
  

  ThermalRatesStoch = zeros(NPESTemp,4);
  for iPES=1:NPESTemp
    for iP = [1,3]
      ThermalRatesStoch(iPES,iP) = sum(ProcessesRatesStochTemp(:,iP,iPES) .* ExpVec(:,1));
    end
  end
  
  for iP = [1,3]
    if (iP == 1)
      fprintf('  Dissociation Thermal Rate Mean = %e cm^s/(mol*s); St. Deviation = %e\n', mean(ThermalRatesStoch(:,iP)), std(ThermalRatesStoch(:,iP)) ) 
    
      iP
      fig1=figure(iFigure);
      fig = gcf;
      if SaveFigs > 0
        screensize = get( groot, 'Screensize' );
        fig.Position=screensize;
        fig.Color='None';
      end 

      histfit(ThermalRatesStoch(:,iP),30)
      
      xlab = xlabel('$K_{Thermal}^{Diss} [cm^3/s]$');
      xlab.Interpreter = 'latex';
      ylab = ylabel('Posterior Samples');
      grid on

      set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
      %daspect([1 1 1])

      if SaveFigs == 1
         FolderPath = strcat(FigDirPath, '/Rates/');
         [status,msg,msgID] = mkdir(FolderPath);
         FileName   = strcat(FolderPath, MoleculesName(1,:),'-Diss_Thermal_Posterior');
         export_fig(FileName, '-pdf')
         close
      elseif SaveFigs == 2
         FolderPath = strcat(FigDirPath, '/Rates/');
         [status,msg,msgID] = mkdir(FolderPath);
         FileName   = strcat(FolderPath, MoleculesName(1,:),'-Diss_Thermal_Posterior.fig');
         saveas(fig,FileName);
      end
      iFigure=iFigure+1;

    elseif (iP == 3)
      fprintf('  Exchange     Thermal Rate Mean = %e cm^s/(mol*s); St. Deviation = %e\n', mean(ThermalRatesStoch(:,iP)), std(ThermalRatesStoch(:,iP)) )
    
      iP
      fig1=figure(iFigure);
      fig = gcf;
      if SaveFigs > 0
        screensize = get( groot, 'Screensize' );
        fig.Position=screensize;
        fig.Color='None';
      end 

      histfit(ThermalRatesStoch(:,iP),30)
      
      xlab = xlabel('$K_{Thermal}^{Exch} [cm^3/s]$');
      xlab.Interpreter = 'latex';
      ylab = ylabel('Posterior Samples');
      grid on

      set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
      %daspect([1 1 1])

      if SaveFigs == 1
         FolderPath = strcat(FigDirPath, '/Rates/');
         [status,msg,msgID] = mkdir(FolderPath);
         FileName   = strcat(FolderPath, MoleculesName(1,:),'-Exch_Thermal_Posterior');
         export_fig(FileName, '-pdf')
         close
      elseif SaveFigs == 2
         FolderPath = strcat(FigDirPath, '/Rates/');
         [status,msg,msgID] = mkdir(FolderPath);
         FileName   = strcat(FolderPath, MoleculesName(1,:),'-Exch_Thermal_Posterior.fig');
         saveas(fig,FileName);
      end
      iFigure=iFigure+1;
      
    end 
  end
  
end