% -- MATLAB --
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

%function iFigure = PlotLevelRates(iFigure, NLevels, iBinsStart, iBinsEnd, RatesMatrix, RatesSigmaMatrix, LevelEeV, Levelvqn)    
function iFigure = PlotLevelRates(iFigure, NLevels, iBinsStart, iBinsEnd, RatesMatrix, LevelEeV, Levelvqn)    

    global vqnColor SigmaOn NTint NBins ColorVec linS TempChar TempChar2 MoleculesName System SaveFigs FilePath PlotExchFlg
           
    NSigma     = 3.d0;
    iBinnedMol = 1;

    for iLevels = iBinsStart:iBinsEnd
      
      fig1=figure(iFigure)

      for iTint=1:NTint
        if SigmaOn(1,:) == 'yes' 
          %errorbar([1:size(RatesMatrix,2)],RatesMatrix(iLevels,:,iTint),RatesSigmaMatrix(iLevels,:,iTint));
          %errorbarlogy;
          semilogy([1:size(RatesMatrix,2)],RatesMatrix(iLevels,:,iTint),'o','Color',ColorVec(1,:),'Marker','o','Markersize',1,'MarkerFaceColor',ColorVec(1,:));
          hold on
          semilogy([[1:size(RatesMatrix,2)];[1:size(RatesMatrix,2)]],[RatesMatrix(iLevels,1:size(RatesMatrix,2),iTint)-NSigma.*RatesSigmaMatrix(iLevels,1:size(RatesMatrix,2),iTint); RatesMatrix(iLevels,1:size(RatesMatrix,2),iTint)+NSigma.*RatesSigmaMatrix(iLevels,1:size(RatesMatrix,2),iTint)],'-','Color',ColorVec(1,:),'LineWidth',0.5);
        else
          if vqnColor == 1
            scatter([1:size(RatesMatrix,2)]',RatesMatrix(iLevels,1:NLevels(iBinnedMol),iTint)',20,Levelvqn(1:NLevels(iBinnedMol),iBinnedMol)+1,'Filled');
          else
            scatter([1:size(RatesMatrix,2)]',RatesMatrix(iLevels,1:NLevels(iBinnedMol),iTint)',20,ColorVec(2,:),'Filled');            
          end
          cmap = colormap(lines(max(Levelvqn(1:NLevels(iBinnedMol),iBinnedMol))+1));
          set(gca, 'YScale', 'log');
          h = colorbar;
        end
        hold on
      end
      %ylim([1.d-18, 2.d-8]);
      %legend(TempChar2,'Location','NorthWest');
      legend boxoff
      xlabel('Processes')
      ylabel(strcat('K_{',num2str(iLevels),',j} [cm^3/s]'));
      grid on

      set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
      set(gcf, 'PaperPositionMode', 'auto');
      if SaveFigs == 1
          FileName = strcat(MoleculesName(1,:),'-Dissociation');
          FilePathNamePng = strcat(FilePath, FileName, '.png');
          FilePathNameFig = strcat(FilePath, FileName, '.fig');
          saveas(gcf,strcat(FilePath, FileName),'pdf')
          saveas(gcf,FilePathNamePng)
          savefig(FilePathNameFig);
      end

      iFigure=iFigure+1;
      
      
      fig2=figure(iFigure)

      for iTint=1:NTint
        if SigmaOn(1,:) == 'yes' 
          %errorbar([1:size(RatesMatrix,2)],RatesMatrix(iLevels,:,iTint),RatesSigmaMatrix(iLevels,:,iTint));
          %errorbarlogy;
          semilogy(LevelEeV(:,1),RatesMatrix(iLevels,:,iTint),'o','Color',ColorVec(1,:),'Marker','o','Markersize',1,'MarkerFaceColor',ColorVec(1,:));
          hold on
          semilogy([LevelEeV(:,1);LevelEeV(:,1)],[RatesMatrix(iLevels,1:size(RatesMatrix,2),iTint)-NSigma.*RatesSigmaMatrix(iLevels,1:size(RatesMatrix,2),iTint); RatesMatrix(iLevels,1:size(RatesMatrix,2),iTint)+NSigma.*RatesSigmaMatrix(iLevels,1:size(RatesMatrix,2),iTint)],'-','Color',ColorVec(1,:),'LineWidth',0.5);
        else
          if vqnColor == 1
            scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),RatesMatrix(iLevels,:,iTint),20,Levelvqn(1:NLevels(iBinnedMol),iBinnedMol)+1,'Filled');
            hold on
            cmap = colormap(lines(max(Levelvqn(1:NLevels(iBinnedMol),iBinnedMol))+1));
            h = colorbar;
            for vqn=0:max(Levelvqn(1:NLevels(iBinnedMol),iBinnedMol))
              xAbs = [];
              yAbs = [];
              for kLevel = 1:NLevels(iBinnedMol)
                 if Levelvqn(kLevel,iBinnedMol) == vqn
                   xAbs = [xAbs; LevelEeV(kLevel,iBinnedMol)];
                   yAbs = [yAbs; RatesMatrix(iLevels,kLevel,iTint)];
                 end 
              end
              plot(xAbs,yAbs,'-','LineWidth',1,'Color',cmap(vqn+1,:));
            end
          else
            scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),RatesMatrix(iLevels,:,iTint),20,ColorVec(2,:),'Filled');            
          end
          set(gca, 'YScale', 'log');
 
        end  
        hold on
      end
      %ylim([1.d-18, 2.d-8]);
      %legend(TempChar2,'Location','NorthWest');
      legend boxoff
      xlabel('Energy [eV]')
      ylabel(strcat('K_{',num2str(iLevels),',j} [cm^3/s]'));
      grid on

      set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
      set(gcf, 'PaperPositionMode', 'auto');
      if SaveFigs == 1
          FileName = strcat(MoleculesName(1,:),'-Dissociation');
          FilePathNamePng = strcat(FilePath, FileName, '.png');
          FilePathNameFig = strcat(FilePath, FileName, '.fig');
          saveas(gcf,strcat(FilePath, FileName),'pdf')
          saveas(gcf,FilePathNamePng)
          savefig(FilePathNameFig);
      end

      iFigure=iFigure+1;
      
    end
    
end
