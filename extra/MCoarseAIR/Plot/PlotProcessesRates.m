%% The Function plots the processes rates in the Energy Space: K(iLevel,iProcess) = Mol1(iLevel) -> Mol2(iProcess)
%
%  Input Arguments:  - iFigure              Index for the first figure to plot
%                    - ProcessesRates       Matrix of Processes Rates
%                    - ProcessesRatesSigma  Matrix of Processes Rates St. Deviations
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

%function iFigure = PlotProcessesRates(iFigure, ProcessesRates, ProcessesRatesSigma, LevelEeV, EeV, LevToBin)   
function iFigure = PlotProcessesRates(iFigure, NLevels, ProcessesRates, Levelg, LevelEeV, EeV, DeltaEintDiss, Levelvqn, Leveljqn, LevToBin, rOut, rIn, Egam, Tau, LevelQ, dVIn, ddVIn, dVOut, ddVOut, dVdJIn, dVdJOut, QBins)   

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
  
  global MergePairsFlg PlotPairUser SigmaOn OverlapFlg LevToBinFlg vqnColor
  global NTint NBins ColorVec linS MoleculesName Pair_To_Molecule 
  global SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz LegendFontNm
  global Plnck UKb Ue KeV AvN AMUToKg EhToeV DSWtoKg ATMToPa
  global T0_Vec
  
  ProcessesRatesTemp      = ProcessesRates;
  %ProcessesRatesSigmaTemp = ProcessesRatesSigma.^2;


  if MergePairsFlg == 1
    for iTint=1:NTint    
      for iP = 3:-1:2
        jP = 1;
        while jP < iP
          if Pair_To_Molecule(iP) == Pair_To_Molecule(jP)
            iP
            jP
            ProcessesRatesTemp(:,jP+1,iTint)      = ProcessesRatesTemp(:,jP+1,iTint)         + ProcessesRatesTemp(:,iP+1,iTint);
            %ProcessesRatesSigmaTemp(:,jP+1,iTint) = ProcessesRatesSigmaTemp(:,jP+1,iTint)    + ProcessesRatesSigmaTemp(:,iP+1,iTint);
            jP = iP
          end
          jP = jP+1;
        end
      end
      %ProcessesRatesSigmaTemp = ProcessesRatesSigmaTemp.^0.5d0;
    end
  end
  
  
  jj = 1;
  for i = 1:size(LevelEeV,1)
%     if ddVOut(i,1) < 0.d0 
      BoundLevels(jj) =  i;
      jj = jj + 1;
%     end
  end
  
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   jqnIn = Leveljqn(1:NBins(1),1);
%   vqnIn = Levelvqn(1:NBins(1),1);
%   
%   figure(iFigure);
%   fig = gcf;
%   screensize = get( groot, 'Screensize' );
%   fig.Position=screensize;
%   fig.Color='None';
%   x = [0:1:max(vqnIn)];
%   y = [0:1:max(jqnIn)];
%   [X,Y] = meshgrid(x,y);
%   
%   DissRate = 0.d0 * ones(max(jqnIn)+1,max(vqnIn)+1);
%   for iLevel = 1:NBins(1)
%     %DissRate(Leveljqn(iLevel)+1,Levelvqn(iLevel)+1) = ProcessesRatesTemp(iLevel,1,1);
%     DissRate(Leveljqn(iLevel)+1,Levelvqn(iLevel)+1) = LevelEeV(iLevel,1,1);
%   end
%   %contourf(X,Y,log10(DissRate),50);%,'ShowText','on');
%   %g=surf(X,Y,log10(DissRate));%,'ShowText','on');
%   g=surf(X,Y,DissRate);%,'ShowText','on');
%   h=colorbar;
%   %ylabel(h, '$log_{10}(K_{(v,J)}^{Diss})$')
%   ylabel(h, '$Energy [eV]$')
%   colormap(jet);
%   %caxis([log10(1.d-14), log10(2.d-8)]);
%   %zlim([-14, log10(2.d-8)])
%   caxis([-15, 5]);
%   zlim([-15, 5])
%   
%   view(-40,20)
% %  shading interp
%   lightangle(150,70)
%   g.FaceLighting = 'flat';
%   g.AmbientStrength = 0.7;
%   g.DiffuseStrength = 0.7;
%   g.SpecularStrength = 0.8;
%   g.SpecularExponent = 25;
%   g.BackFaceLighting = 'unlit';
%   
%   xlab = xlabel('v');
%   xlab.Interpreter = 'latex';
%   ylab = ylabel('J');
%   ylab.Interpreter = 'latex';
%   %zlab = zlabel('$log_{10}(K_{(v,J)}^{Diss})$');
%   %zlab.Interpreter = 'latex';
% %   if iP == 1
% %     ylab = ylabel('$K_{i}^{Diss} [cm^3/s]$');
% %   else
% %     if Pair_To_Molecule(iP-1) == Pair_To_Molecule(1)
% %       ystr = strcat('$K_{i,',MoleculesName(Pair_To_Molecule(iP-1),:),'}^{Exch} [cm^3/s]$')
% %       ylab = ylabel(ystr);
% %     else
% %       ystr = strcat('$K_{i,',MoleculesName(Pair_To_Molecule(iP-1),:),'}^{Exch} [cm^3/s]$')
% %       ylab = ylabel(ystr);
% %     end
% %   end
% %   ylab.Interpreter = 'latex';
% %   grid on
% 
%   set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
%   %daspect([1 1 1])
%   pause
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if OverlapFlg == 0
    
    for iP=[1,3,4]
      if (PlotPairUser(max(iP-1,1)) == 1) 

        for iTint=1:NTint

          fig1=figure(iFigure);
          fig = gcf;
          if SaveFigs > 0
            screensize = get( groot, 'Screensize' );
            fig.Position=screensize;
            fig.Color='None';
          end 

          if SigmaOn(1,:) == 'yes' 
            errorbar(EeV(1:NBins(1),1),ProcessesRatesTemp(:,iP,iTint),ProcessesRatesSigmaTemp(:,iP,iTint),'Color',ColorVec(1,:),'linestyle',linS{1},'LineWidth',3);%,'Marker','o','Markersize',6);
          else
            %plot(EeV(1:NBins(1),1),ProcessesRatesTemp(:,1,iTint),'o','Color',ColorVec(1),'Marker','o','Markersize',6);
            if LevToBinFlg == 1
              cmap = colormap(lines(max(LevToBin(1:NBins(1),1))));
              %scatter(DeltaEintDiss(1:NBins(1),1),ProcessesRatesTemp(:,iP,iTint),20,cmap(LevToBin(1:NBins(1),1)));
              scatter(EeV(1:NBins(1),1),ProcessesRatesTemp(:,iP,iTint),20,cmap(LevToBin(1:NBins(1),1)));         
              set(gca, 'YScale', 'log');
              h=colorbar;
            elseif vqnColor == 1
              %h=scatter(EeV(1:NBins(1),1)+LevelEeV(1,1),ProcessesRatesTemp(:,iP,iTint),20,Levelvqn(1:NBins(1),1),'filled');
              %h=scatter(Levelvqn(1:NBins(1),1),ProcessesRatesTemp(:,iP,iTint),20,Levelvqn(1:NBins(1),1));
              %h=scatter(DeltaEintDiss(BoundLevels(:),1),ProcessesRatesTemp(BoundLevels(:),iP,iTint),20,dVOut(BoundLevels(:),1));
              %h=scatter3(rIn(1:NBins(1),1),Leveljqn(1:NBins(1),1),EeV(1:NBins(1),1)+LevelEeV(1,1),20,log10(ProcessesRatesTemp(1:NBins(1),iP,iTint)),'filled');
              %hold on
              %h=scatter3(rOut(1:NBins(1),1),Leveljqn(1:NBins(1),1),EeV(1:NBins(1),1)+LevelEeV(1,1),20,log10(ProcessesRatesTemp(1:NBins(1),iP,iTint)),'filled');
              h=scatter(rOut(1:NBins(1),1)-rIn(1:NBins(1),1),ProcessesRatesTemp(:,iP,iTint),20,Levelvqn(1:NBins(1),1));
              %h=scatter(DeltaEintDiss(BoundLevels(:),1),ProcessesRatesTemp(BoundLevels(:),iP,iTint),20,dVOut(BoundLevels(:),1),'filled');
              %h=scatter(EeV(1:NBins(1),1)+LevelEeV(1,1),ProcessesRatesTemp(:,iP,iTint),20,Leveljqn(1:NBins(1),1));
              %h=scatter(DeltaEintDiss(1:NBins(1),1),ProcessesRatesTemp(:,iP,iTint),20,log10(Egam(1:NBins(1),1)));
              %h=scatter(DeltaEintDiss(1:NBins(1),1),ProcessesRatesTemp(:,iP,iTint),20,LevelEeV(1:NBins(1),1));
              set(gca, 'YScale', 'log');
              h=colorbar;
              ylabel(h, '$log_{10}(K_{i}^{Diss}) [cm^3/s]$')
              colormap(jet);
              %ylim(h,[0 0.1]);
              %caxis([2.4e-10 4.d-8]);
            else
              size(ProcessesRatesTemp)
              size(ColorVec)
              iP
              semilogy(EeV(1:NBins(1),1),ProcessesRatesTemp(1:NBins(1),iP,iTint),'o','Color',ColorVec(1,:),'Marker','o','Markersize',3,'MarkerFaceColor',ColorVec(1,:));
              %semilogy(DeltaEintDiss(1:NBins(1),1),ProcessesRatesTemp(:,iP,iTint),'o','Color',ColorVec(iP,:),'Marker','o','Markersize',3,'MarkerFaceColor',ColorVec(iP,:));
            end
          end
          hold on
          %ylim([1.d-14, 2.d-8]);
          %legend(TempChar2,'Location','NorthWest');
          %legend boxoff

          xlab = xlabel('$E_i [eV]$');
          xlab.Interpreter = 'latex';
          ylab = ylabel('$K_{i}^{Diss} [cm^3/s]$');
          ylab.Interpreter = 'latex';
          grid on

%           xlab = xlabel('R [a_0]');
%           xlab.Interpreter = 'latex';
%           ylab = ylabel('J');
%           ylab.Interpreter = 'latex';
%           zlab = zlabel('E [eV]');
%           zlab.Interpreter = 'latex';
%           grid on
          
%           SumContJ = zeros(1,max(Leveljqn(:,1))+1);
%           SumContV = zeros(1,max(Levelvqn(:,1))+1);
%           SumTemp  = sum(ProcessesRatesTemp(:,iP,iTint));
%           if iP == 1
%             FileName = strcat('./',MoleculesName(Pair_To_Molecule(1),:),'_DissRate.csv');
%             fileID = fopen(FileName,'w');
%             fprintf(fileID,'Variables = "R", "J", "E", "KDiss", "Q", "KDiss*Q", "DeltaEintDiss"\n');
%             for i=1:NBins(1)
%                 %if (DeltaEintDiss(i,1) < -0.0d0 && DeltaEintDiss(i,1) > -5.d0) 
%                 if ProcessesRatesTemp(i,iP,iTint) >= 1.d-13
%                   fprintf(fileID,'%15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e\n',rIn(i,1), Leveljqn(i,1), EeV(i,1)+LevelEeV(1,1), ProcessesRatesTemp(i,iP,iTint), LevelQ(i,1), ProcessesRatesTemp(i,iP,iTint)*LevelQ(i,1)./SumTemp, DeltaEintDiss(i,1));
%                 end
%             end
%             fclose(fileID);
%           end


%           xlab = xlabel('Energy [eV]');
%           xlab.Interpreter = 'latex';
%           if iP == 1
%             ylab = ylabel('$K_{i}^{Diss} [cm^3/s]$');
%           else
%             if Pair_To_Molecule(iP-1) == Pair_To_Molecule(1)
%               ystr = strcat('$K_{i,',MoleculesName(Pair_To_Molecule(iP-1),:),'}^{Exch} [cm^3/s]$')
%               ylab = ylabel(ystr);
%             else
%               ystr = strcat('$K_{i,',MoleculesName(Pair_To_Molecule(iP-1),:),'}^{Exch} [cm^3/s]$')
%               ylab = ylabel(ystr);
%             end
%           end
%           ylab.Interpreter = 'latex';
%           grid on

          set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
          %daspect([1 1 1])
          
          if SaveFigs == 1
             FolderPath = strcat(FigDirPath, '/Rates/');
             [status,msg,msgID] = mkdir(FolderPath);
             FileName   = strcat(FolderPath, MoleculesName(1,:),'-OverallToPair', num2str(iP-1));
             export_fig(FileName, '-pdf')
             close
          elseif SaveFigs == 2
             FolderPath = strcat(FigDirPath, '/Rates/');
             [status,msg,msgID] = mkdir(FolderPath);
             FileName   = strcat(FolderPath, MoleculesName(1,:),'-OverallToPair', num2str(iP-1), '.fig');
             saveas(fig,FileName);
          end
          iFigure=iFigure+1;
            
          
          
%           CDFContJ(1) = 0.d0
%           for i = 2:length(SumContJ)
%             CDFContJ(i) = CDFContJ(i-1) + SumContJ(i);
%           end
% 
%           fig1=figure(iFigure);
%           fig = gcf;
%           if SaveFigs > 0
%             screensize = get( groot, 'Screensize' );
%             fig.Position=screensize;
%             fig.Color='None';
%           end 
%           yyaxis left
%           bar(SumContJ)
%           ylab = ylabel('Contribution to Dissociation [%]');
%           yyaxis right
%           plot(CDFContJ)
%           xlab = xlabel('J');
%           xlab.Interpreter = 'latex';
%           ylab = ylabel('Cumulative Contribution to Dissociation [%]');
%           ylab.Interpreter = 'latex';
%           set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
%           iFigure=iFigure+1;
%           
%           
%           CDFContV(1) = 0.d0
%           for i = 2:length(SumContV)
%             CDFContV(i) = CDFContV(i-1) + SumContV(i);
%           end
% 
%           fig1=figure(iFigure);
%           fig = gcf;
%           if SaveFigs > 0
%             screensize = get( groot, 'Screensize' );
%             fig.Position=screensize;
%             fig.Color='None';
%           end 
%           yyaxis left
%           bar(SumContV)
%           ylab = ylabel('Contribution to Dissociation [%]');
%           ylab.Interpreter = 'latex';
%           yyaxis right
%           plot(CDFContV)
%           xlab = xlabel('v');
%           xlab.Interpreter = 'latex';
%           ylab = ylabel('Cumulative Contribution to Dissociation [%]');
%           ylab.Interpreter = 'latex';
%           set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
%           iFigure=iFigure+1;
          
          
        end

      end

    end

  
  else
    
    fig1=figure(iFigure);
    
    if SaveFigs > 0
      screensize = get( groot, 'Screensize' );
      fig.Position=screensize;
      fig.Color='None';
    end 
    
    for iP=[1,3,4]
    if (PlotPairUser(max(iP-1,1)) == 1) 

      for iTint=1:NTint
    
        if SigmaOn(1,:) == 'yes' 
          errorbar(EeV(1:NBins(1),1)+LevelEeV(1,1),ProcessesRatesTemp(:,iP,iTint),ProcessesRatesSigmaTemp(:,iP,iTint),'Color',ColorVec(1,:),'linestyle',linS{1},'LineWidth',3);%,'Marker','o','Markersize',6);
        else
          %plot(EeV(1:NBins(1),1),ProcessesRatesTemp(:,1,iTint),'o','Color',ColorVec(1),'Marker','o','Markersize',6);
          semilogy(EeV(1:NBins(1),1)+LevelEeV(1,1),ProcessesRatesTemp(:,iP,iTint),'o','Color',ColorVec(iP,:),'Marker','o','Markersize',3,'MarkerFaceColor',ColorVec(iP,:));
          %semilogy(DeltaEintDiss(1:NBins(1),1),ProcessesRatesTemp(:,iP,iTint),'o','Color',ColorVec(iP,:),'Marker','o','Markersize',3,'MarkerFaceColor',ColorVec(iP,:));
        end
        hold on
        ylim([1.d-14, 2.d-8]);
        if iP == 1
          legendstr = {'$K_{i}^{Diss}$'};
        else
          if Pair_To_Molecule(iP-1) == Pair_To_Molecule(1)
            ystr = strcat('$K_{i,',MoleculesName(Pair_To_Molecule(iP-1),:),'}^{Exch}$')
            legendstr = [legendstr,ystr];
          else
            ystr = strcat('$K_{i,',MoleculesName(Pair_To_Molecule(iP-1),:),'}^{Exch}$')
            legendstr = [legendstr,ystr];
          end
        end

      end

    end
    
    legendstr
    
    clab = legend(legendstr,'Location','Best','Interpreter', 'latex');
    clab.Interpreter = 'latex';
    set(clab,'FontSize',LegendFontSz,'FontName',LegendFontNm,'Interpreter','latex');
        
    %legend boxoff
    xlab = xlabel('Energy [eV]');
    xlab.Interpreter = 'latex';
    ylab = ylabel('$K_{i,...} [cm^3/s]$');
    ylab.Interpreter = 'latex';
    grid on

    set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
    
    if SaveFigs == 1
       FolderPath = strcat(FigDirPath, '/Rates/');
       [status,msg,msgID] = mkdir(FolderPath);
       FileName   = strcat(FolderPath, MoleculesName(1,:),'-OverallToPair', num2str(iP-1));
       export_fig(FileName, '-pdf')
       close
    elseif SaveFigs == 2
       FolderPath = strcat(FigDirPath, '/Rates/');
       [status,msg,msgID] = mkdir(FolderPath);
       FileName   = strcat(FolderPath, MoleculesName(1,:),'-OverallToPair', num2str(iP-1), '.fig');
       saveas(fig,FileName);
    end
    
    iFigure=iFigure+1;

    end
  
  end
  
  
  iBinnedMol = 1;
  ExpVec(:,1)  = QBins(:,1);
  ExpVec                           = ExpVec ./ sum(ExpVec);
  
  ThermalRates = zeros(4,0);
  for iP = 1:4
    ThermalRates(iP) = sum(ProcessesRates(:,iP,1) .* ExpVec(:,1));
    if (iP == 1)
      fprintf('  Dissociation Thermal Rate = %e cm^s/(mol*s)\n', ThermalRates(iP))
    elseif (iP == 2)
      fprintf('  Inelastic    Thermal Rate = %e cm^s/(mol*s)\n', ThermalRates(iP))
    elseif (iP == 3)
      fprintf('  Exchange 1   Thermal Rate = %e cm^s/(mol*s)\n', ThermalRates(iP))
    elseif (iP == 4)
      fprintf('  Exchange 2   Thermal Rate = %e cm^s/(mol*s)\n', ThermalRates(iP))
    end 
    
  end
  
  
end