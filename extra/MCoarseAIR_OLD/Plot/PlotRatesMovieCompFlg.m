% -- MATLAB --
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

function [iFigure] = PlotRatesMovieComp(iFigure, iP, EeV, RatesMatrix)

  global SigmaOn NTint NBins ColorVec linS TempChar TempChar2 MoleculesName System SaveFigs FilePath StartBin FinalBin Pair_Name

  fig = figure(iFigure);
  pause(1)
  for iBins = 1:NBins(1)

    left_color  = [0 0 0];
    right_color = [1 0 0];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);

    yyaxis left 
    for iT=1:NTint
        semilogy(EeV(1:iBins,1),RatesMatrix(iBins,1:iBins,1,iT),'o','Color','k','linestyle',linS{iT},'LineWidth',3,'Markersize',5,'MarkerFaceColor','auto')
        hold on
    end
    %xlim([EeV(1,1), EeV(NBins(1),1)]);
    %xlabel('CO(j) Bin Energy [eV]');
    ylim([1.d-14, 1.d-7]);
    ylabel(strcat('K_{CO(',num2str(iBins),') -> CO(j)}'));
    hold off

    yyaxis right
    for iT=1:NTint
        semilogy(EeV(:,2),RatesMatrix(iBins,:,3,iT)','o','Color',[1 0 0],'linestyle',linS{iT},'LineWidth',3,'Markersize',5,'MarkerFaceColor','auto')
        hold on
    end
    %xlim([EeV(1,2), EeV(NBins(1),2)]);
    %xlabel('O_2(k) Bin Energy [eV]');
    xlim([EeV(1,1), EeV(NBins(1),1)]);
    xlabel(strcat('Final Bin Energy [eV]  (CO(',num2str(iBins),') Bin Energy = ',num2str(EeV(iBins,1)),'eV)'));
    ylim([1.d-14, 1.d-7]);
    ylabel(strcat('K_{CO(',num2str(iBins),') -> O_2(k)}'));

    legend(TempChar2,'Location','NorthWest');
    legend boxoff
    set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
    set(gcf, 'PaperPositionMode', 'auto');
    pause(0.5)
    hold off
    
  end
  if SaveFigs == 1
    FileName = strcat(MoleculesName(1,:),'-',Pair_Name(iP,:),'-',TempChar(iT,5:end));
    FilePathNamePng = strcat(FilePath, FileName, '.png');
    FilePathNameFig = strcat(FilePath, FileName, '.fig');
    saveas(gcf,strcat(FilePath, FileName),'pdf')
    saveas(gcf,FilePathNamePng)
    savefig(FilePathNameFig);
  end
  iFigure=iFigure+1;
  
end
