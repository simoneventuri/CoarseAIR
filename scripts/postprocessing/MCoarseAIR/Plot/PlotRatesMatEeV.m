%% The Function plots the rates in the Energy Space: K(iLevel,jLevel) = Mol1(iLevel) -> Mol2(jLevel)
%
%  Input Arguments:  - iFigure         Index for the first figure to plot
%                    - RatesMatrix     Matrix of Rates
%                    - ProcessesRates  Matrix of Processes Rates
%                    - Levelg          Levels/Bins' degeneracy
%                    - LevelEeV        Levels energy in eV, with the 0 corresponding to the dissociation energy
%                    - EeV             Levels/Bins' energy in eV, with the 0 corresponding to the dissociation energy
%   
%  Input Global Var: - PlotPairUser    Vector of flags for the pairs to plot (e.g.: [1 0 0])
%

function [iFigure] = PlotRatesMatEeV(iFigure, RatesMatrix, ProcessesRates, Levelg, LevelEeV, EeV)

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
  
  global PlotPairUser EeVSpace
  global NTint NBins TempChar MoleculesName Pair_Name Pair_To_Molecule PlotPair
  global KeV
  global SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz LegendFontNm

  for iTint=1:NTint
    
    
    RatesMatrixTemp = RatesMatrix;
    
    for iP = 3:2
      jP = 1;
      while jP < iP
        if (Pair_To_Molecule(iP) == Pair_To_Molecule(jP))
          RatesMatrixTemp(1:NBins(Pair_To_Molecule(iP)),1:NBins(Pair_To_Molecule(iP))) = RatesMatrixTemp(1:NBins(Pair_To_Molecule(iP)),1:NBins(Pair_To_Molecule(iP))) + RatesMatrix(EeVOrder1(1:NBins(Pair_To_Molecule(iP))),EeVOrder1(1:NBins(Pair_To_Molecule(iP))),iP,iTint);
          jP = iP;
        end
        jP = jP + 1;
      end
    end
    
%     for iLev = 1:size(RatesMatrix,1)
%       for jLev = 1:size(RatesMatrix,2)
%         if iLev > jLev
%           RatesMatrixExo(jLev,iLev)  = RatesMatrixAll(jLev,iLev);
%         elseif iLev < jLev
%           RatesMatrixEndo(jLev,iLev) = RatesMatrixAll(jLev,iLev);
%           %RatesMatrixEndo(jLev,iLev) = RatesMatrixExo(iLev,jLev) * Levelg(iLev) / Levelg(jLev) * exp( - (LevelEeV(iLev) - LevelEeV(jLev)) / (KeV*T0_Vec(iTint)) ) ;
%         end
%       end
%     end
    
%     figure(iFigure)
%     h = imagesc(log10(RatesMatrixExo));
%     str_y = ['Energy of ', MoleculesName(1,:), '(j) [eV]'];
%     str_x = ['Energy of ', MoleculesName(1,:), '(i) [eV]'];
%     h=colorbar;
%     colormap(jet)
%     hold on
%     ylabel(h, 'log_{10}(K_{ij})')
%     ylim(h,[-15 -8]);
%     caxis([-15 -8]);
%     xlabel(str_x);
%     ylabel(str_y);
%     set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
%     set(gcf, 'PaperPositionMode', 'auto');
%     if SaveFigs == 1
%         FileName = strcat(MoleculesName(1,:),'-',Pair_Name(iP,:),'-',TempChar(iTint,5:end));
%         FilePathNamePng = strcat(FilePath, FileName, '.png');
%         FilePathNameFig = strcat(FilePath, FileName, '.fig');
%         saveas(gcf,strcat(FilePath, FileName),'pdf')
%         saveas(gcf,FilePathNamePng)
%         savefig(FilePathNameFig);
%     end
%     iFigure=iFigure+1;
%     
%     
%     figure(iFigure)
%     h = imagesc(log10(RatesMatrixEndo));
%     str_y = ['Energy of ', MoleculesName(1,:), '(j) [eV]'];
%     str_x = ['Energy of ', MoleculesName(1,:), '(i) [eV]'];
%     h=colorbar;
%     colormap(jet)
%     hold on
%     ylabel(h, 'log_{10}(K_{ij})')
%     ylim(h,[-15 -8]);
%     caxis([-15 -8]);
%     xlabel(str_x);
%     ylabel(str_y);
%     set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
%     set(gcf, 'PaperPositionMode', 'auto');
%     if SaveFigs == 1
%         FileName = strcat(MoleculesName(1,:),'-',Pair_Name(iP,:),'-',TempChar(iTint,5:end));
%         FilePathNamePng = strcat(FilePath, FileName, '.png');
%         FilePathNameFig = strcat(FilePath, FileName, '.fig');
%         saveas(gcf,strcat(FilePath, FileName),'pdf')
%         saveas(gcf,FilePathNamePng)
%         savefig(FilePathNameFig);
%     end
%     iFigure=iFigure+1;

    for iP = 1:3      
      if PlotPairUser(iP) == 1
        
        figure(iFigure);
        fig = gcf;
        screensize = get( groot, 'Screensize' );
        fig.Position=screensize;
        fig.Color='None';
        
        iMol = Pair_To_Molecule(iP);
        [EeVSorted1, EeVOrder1] = sort(EeV(1:NBins(1),1));
        EeVSorted1 = EeVSorted1 + LevelEeV(1,1);
        [EeVSorted2, EeVOrder2] = sort(EeV(1:NBins(iMol),iMol));
        EeVSorted2 = EeVSorted2 + LevelEeV(1,iMol);
%         if EeVSpace == 1
%           imagesc('YData',EeVSorted1(1:NBins(1)),'XData',EeVSorted2(1:NBins(iMol)),'CData',log10(RatesMatrixTemp(EeVOrder1(1:NBins(1)),EeVOrder2(1:NBins(iMol)),iP,iTint)))
%           ylim([min(EeVSorted1), max(EeVSorted1)]);
%           xlim([min(EeVSorted2), max(EeVSorted2)]);
%         else
          imagesc(log10(RatesMatrixTemp(EeVOrder1(1:NBins(1)),EeVOrder2(1:NBins(iMol)),iP,iTint)))
%         end
        hold on

        str_y = [MoleculesName(1,:),'(i)'];
        str_x = [MoleculesName(iMol,:),'(k)'];
        
        cb=colorbar;
        ylim(cb,[-15 -10.5]);
        ylab = ylabel(cb, '$log_{10}(K_{ik})$');
        ylab.Interpreter = 'latex';
        set(cb,'FontSize',LegendFontSz,'FontName',LegendFontNm,'TickLabelInterpreter','latex');
        caxis([-15 -10.5]);
        colormap(jet)

        xlab = xlabel(str_x);
        xlab.Interpreter = 'latex';
        ylab = ylabel(str_y);
        ylab.Interpreter = 'latex';
        
%         tit = title(str_title);
%         tit.Interpreter = 'latex';
            
        set(gca,'FontSize',AxisFontSz,'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
        daspect([1 1 1])
        
        if SaveFigs == 1
           FolderPath = strcat(FigDirPath, '/Rates/');
           [status,msg,msgID] = mkdir(FolderPath);
           FileName   = strcat(FolderPath,'RatesMatrix-',MoleculesName(1,:),'-To-',MoleculesName(Pair_To_Molecule(iP),:));
           export_fig(FileName, '-pdf')
           close
        elseif SaveFigs == 2
          FolderPath = strcat(FigDirPath, '/Rates/');
          [status,msg,msgID] = mkdir(FolderPath);
          FileName   = strcat(FolderPath,'RatesMatrix-',MoleculesName(1,:),'-To-',MoleculesName(Pair_To_Molecule(iP),:),'.fig');
          savefig(FileName)
          close
        end
        iFigure=iFigure+1;
        
      end
    end

    for iP = 1:3      
      if PlotPair(iP) == 1
        
        figure(iFigure);
        fig = gcf;
        screensize = get( groot, 'Screensize' );
        fig.Position=screensize;
        fig.Color='None';
        
        iMol = Pair_To_Molecule(iP);
        [EeVSorted1, EeVOrder1] = sort(EeV(1:NBins(1),1));
        EeVSorted1 = EeVSorted1 + LevelEeV(1,1);
        [EeVSorted2, EeVOrder2] = sort(EeV(1:NBins(iMol),iMol));
        EeVSorted2 = EeVSorted2 + LevelEeV(1,iMol);
%         if EeVSpace == 1
%           imagesc('YData',EeVSorted1(1:NBins(1)),'XData',EeVSorted2(1:NBins(iMol)),'CData',log10(RatesMatrix(EeVOrder1(1:NBins(1)),EeVOrder2(1:NBins(iMol)),iP,iTint)))
%           ylim([min(EeVSorted1), max(EeVSorted1)]);
%           xlim([min(EeVSorted2), max(EeVSorted2)]);
%         else
          imagesc(log10(RatesMatrix(EeVOrder1(1:NBins(1)),EeVOrder2(1:NBins(iMol)),iP,iTint)))
%         end 
        hold on

        str_y = [MoleculesName(1,:),'(i)'];
        str_x = [MoleculesName(iMol,:),'(k)'];
        
        cb=colorbar;
        ylim(cb,[-15 -10.5]);
        ylab = ylabel(cb, '$log_{10}(K_{ik})$');
        ylab.Interpreter = 'latex';
        set(cb,'FontSize',LegendFontSz,'FontName',LegendFontNm,'TickLabelInterpreter','latex');
        caxis([-15 -10.5]);
        colormap(jet)

        xlab = xlabel(str_x);
        xlab.Interpreter = 'latex';
        ylab = ylabel(str_y);
        ylab.Interpreter = 'latex';
        
%         tit = title(str_title);
%         tit.Interpreter = 'latex';
            
        set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
        daspect([1 1 1])
        
        if SaveFigs == 1
           FolderPath = strcat(FigDirPath, '/Rates/');
           [status,msg,msgID] = mkdir(FolderPath);
           FileName   = strcat(FolderPath,'RatesMatrix-',MoleculesName(1,:),'-To-Pair',num2str(iP));
           export_fig(FileName, '-pdf')
           close
        elseif SaveFigs == 2
          FolderPath = strcat(FigDirPath, '/Rates/');
          [status,msg,msgID] = mkdir(FolderPath);
           FileName   = strcat(FolderPath,'RatesMatrix-',MoleculesName(1,:),'-To-Pair',num2str(iP),'.fig');
          savefig(FileName)
          close
        end
        iFigure=iFigure+1;
        
      end 
    end
    
    
%     figure(iFigure)
%     [EeVSorted1, EeVOrder1] = sort(EeV(1:NBins(1),1)-abs(LevelEeV(1,1)));
%     DissMatrix = diag(ProcessesRates(EeVOrder1(1:NBins(1)),1,1));
%     imagesc(log10(DissMatrix.'.*6.022e23));
%     str_y = ['Energy of ', MoleculesName(1,:), '(j) [eV]'];
%     str_x = ['Energy of ', MoleculesName(1,:), '(i) [eV]'];
%     h=colorbar;
%     colormap(jet)
%     hold on
%     ylabel(h, 'log_{10}(K_{ij})')
%     ylim(h,[0 14.3]);
%     caxis([0  14.3]);
%     xlabel(str_x);
%     ylabel(str_y);
%     set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
%     set(gcf, 'PaperPositionMode', 'auto');
%     if SaveFigs == 1
%         FileName = strcat(MoleculesName(1,:),'-',Pair_Name(iP,:),'-',TempChar(iTint,5:end));
%         FilePathNamePng = strcat(FilePath, FileName, '.png');
%         FilePathNameFig = strcat(FilePath, FileName, '.fig');
%         saveas(gcf,strcat(FilePath, FileName),'pdf')
%         saveas(gcf,FilePathNamePng)
%         savefig(FilePathNameFig);
%     end
%     iFigure=iFigure+1;
    
    
  end
  
end
