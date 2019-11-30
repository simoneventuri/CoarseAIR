%% The Function plots the rates in the Q.N.s' Space: K(iLevel,(v,J)) = Mol1(iLevel) -> Mol2(v,J), where iLevel is specified by iLevelsVec 
%
%  Input Arguments:  - iFigure:      Index for the first figure to plot
%                    - LevelEeV:     Levels/Bins' energy in eV, with the 0 corresponding to the dissociation energy
%                    - Leveljqn:     Levels/Bins' rotational Q.N.s
%                    - Levelvqn:     Levels/Bins' vibrational Q.N.s
%                    - RatesMatrix:  Matrix of Rates
%
%  Input Global Var: - iLevelsVec    Vector of levels/bins indexes for which we want to plot rates 
%                    - PlotPairUser  Vector of flags for the pairs to plot (e.g.: [1 0 0])
%                    - SubplotsFlg   Flag 0/1 for having multiple plots in the same figure 
%

function [iFigure] = PlotStSRatesQN(iFigure, EeV, LevelEeV, Leveljqn, Levelvqn, RatesMatrix, rIn, rOut)

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

  global SubplotsFlg PlotPairUser iLevelsVec
  global SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz LegendFontNm
  global NTint NBins MoleculesName Pair_To_Molecule TempChar
       
  jqnIn = Leveljqn(1:NBins(1),1);
  vqnIn = Levelvqn(1:NBins(1),1);
  
  figure(iFigure);
  fig = gcf;
  screensize = get( groot, 'Screensize' );
  fig.Position=screensize;
  fig.Color='None';
  x = [0:1:max(vqnIn)];
  y = [0:1:max(jqnIn)];
  [X,Y] = meshgrid(x,y);
  LevelToBin = max(LevelEeV(:,1)) * ones(max(jqnIn)+1,max(vqnIn)+1);
  for iLevel = 1:NBins(1)
    LevelToBin(Leveljqn(iLevel)+1,Levelvqn(iLevel)+1)=LevelEeV(iLevel,1);
  end
  contourf(X,Y,LevelToBin,100);%,'ShowText','on');
  %surf(X,Y,LevelToBin);
  clear x y X Y 
  cb=colorbar;
  ylabel(cb, 'Energy [eV]')
  set(cb,'FontSize',LegendFontSz,'FontName',LegendFontNm,'TickLabelInterpreter','latex');
  cb.Label.Interpreter = 'latex';
  ylim(cb,[min(LevelEeV(:,1)) max(LevelEeV(:,1))]);
  caxis([min(LevelEeV(:,1)) max(LevelEeV(:,1))]);
  colormap(jet)
  hold on
  str_x = ['v'];
  xlab = xlabel(str_x);
  xlab.Interpreter = 'latex';
  str_y = ['J'];
  ylab = ylabel(str_y);
  ylab.Interpreter = 'latex';
  str_title = strcat(MoleculesName(1,:),{' '},'Energy [eV]');
  %title(str_title);
  daspect([1 1 1])
  %pbaspect([1 1 1])
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  if SaveFigs == 1
     FolderPath = strcat(FigDirPath, '/Rates/');
     [status,msg,msgID] = mkdir(FolderPath);
     FileName   = strcat(FolderPath, MoleculesName(1,:),'-QNs-Energy');
     export_fig(FileName, '-pdf')
     close
  elseif SaveFigs == 2
    FolderPath = strcat(FigDirPath, '/Rates/');
    [status,msg,msgID] = mkdir(FolderPath);
    FileName   = strcat(FolderPath, MoleculesName(1,:),'-QNs-Energy.fig');
    savefig(FileName)
    close
  end
  iFigure=iFigure+1;
  
  
  filename = strcat('./DiatPot.csv');
  fileID = fopen(filename,'w');
  fprintf(fileID,'Variables = "R", "J", "E"\n');
  for jBins = 1:NBins(1)
    fprintf(fileID,'%15.6e, %15.6e, %15.6e\n',rIn(jBins,1),Leveljqn(jBins,1),EeV(jBins,1)+LevelEeV(1,1));
    fprintf(fileID,'%15.6e, %15.6e, %15.6e\n',rOut(jBins,1),Leveljqn(jBins,1),EeV(jBins,1)+LevelEeV(1,1));
  end
  fclose(fileID);
  
  
  for iTint = 1:NTint
    
    if SubplotsFlg == 1

      for iBins = iLevelsVec

        jqnIn = Leveljqn(1:NBins(1),1);
        vqnIn = Levelvqn(1:NBins(1),1);
        
        figure(iFigure);
        fig = gcf;
        screensize = get( groot, 'Screensize' );
        fig.Position=[1 1 900 900];
        fig.Color='None';
        
        jP = 0;
        for iP = 1:3    
          if PlotPairUser(iP) == 1
            
            jP = jP + 1;
            subplot(2,sum(PlotPairUser),[jP, jP+sum(PlotPairUser)])
            iMol = Pair_To_Molecule(iP);

            jqn = Leveljqn(1:NBins(iMol),iMol);
            vqn = Levelvqn(1:NBins(iMol),iMol);
            RatesMatrix_Temp2 = zeros(max(jqn)+1,max(vqn)+1);

            for jBins = 1:NBins(iMol)
                RatesMatrix_Temp2(jqn(jBins)+1,vqn(jBins)+1) = RatesMatrix(iBins,jBins,iP,iTint);
            end

            x = [0:1:max(vqn)];
            y = [0:1:max(jqn)];
            [X,Y] = meshgrid(x,y);
            %contourf(X,Y,log10(RatesMatrix_Temp2),50);
            %surf(X,Y,log10(RatesMatrix_Temp2));
            %az = 60;
            %el = 20;
            %view(az, el);
            pcolor(X,Y,log10(RatesMatrix_Temp2));
                   
            if iP == 1
              filename = strcat('./InelRate_v',num2str(Levelvqn(iBins,1)),'_j',num2str(Leveljqn(iBins,1)),'.csv')
              fileID = fopen(filename,'w');
              fprintf(fileID,'Variables = "R", "J", "E", "KInel"\n');
              for jBins = 1:NBins(1)
                if RatesMatrix(iBins,jBins,1,iTint) >= 5.d-13
                  fprintf(fileID,'%15.6e, %15.6e, %15.6e, %15.6e\n',rIn(jBins,1),Leveljqn(jBins,1),EeV(jBins,1)+LevelEeV(1,1),RatesMatrix(iBins,jBins,1,iTint));
                end
              end
              fclose(fileID);
              
              filename = strcat('./ExchRate_v',num2str(Levelvqn(iBins,1)),'_j',num2str(Leveljqn(iBins,1)),'.csv');
              fileID = fopen(filename,'w');
              fprintf(fileID,'Variables = "R", "J", "E", "KExch"\n');
              for jBins = 1:NBins(1)
                if RatesMatrix(iBins,jBins,2,iTint) >= 5.d-13
                  fprintf(fileID,'%15.6e, %15.6e, %15.6e, %15.6e\n',rIn(jBins,1),Leveljqn(jBins,1),EeV(jBins,1)+LevelEeV(1,1),RatesMatrix(iBins,jBins,2,iTint));
                end
              end
              fclose(fileID);
            end

%             clear x y X Y 
%             if jP==sum(PlotPairUser)
%               %title('With Exchange');
%               cb=colorbar;
%               ylim(cb,[-15 log10(5.e-10)]);
%               cbstr = num2str(iBins);
%               cbstr = ['$log_{10}(K_{',cbstr,',(v,J)}) [cm^3/s]$'];
%               ylabel(cb, cbstr);
%               set(cb,'FontSize',LegendFontSz,'FontName',LegendFontNm,'TickLabelInterpreter','latex');
%               cb.Label.Interpreter = 'latex';
%             end
%             %ylim(h,[-15 log10(5.e-10)]);
%             colormap(jet);
%             caxis([-15 log10(5.e-10)]);
%             hold on
%             plot(vqn(iBins,1),jqn(iBins,1),'ko','MarkerSize',10,'MarkerFaceColor','k')   
%                         
%             str_x = ['v'];
%             xlab  = xlabel(str_x);
%             xlab.Interpreter = 'latex';
%             
%             str_y = ['J'];
%             ylab  = ylabel(str_y);
%             ylab.Interpreter = 'latex';
% 
%             str_title = strcat('$',MoleculesName(1,:),'(',num2str(vqnIn(iBins)),',',num2str(jqnIn(iBins)),')',{' '},'->',{' '},MoleculesName(iMol,:),'(v,j)$');
%             %tit = title(str_title);
%             %tit.Interpreter = 'latex';
%             
%             %shading(gca,'interp')
%             set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
%             daspect([1 1 1])
%             %pbaspect([1 1 1])
            
          end
        end
        
        if SaveFigs == 1
          FolderPath = strcat(FigDirPath, '/Rates/');
          [status,msg,msgID] = mkdir(FolderPath);
          FileName   = strcat(FolderPath, 'Level',num2str(iBins),'-QNs');
          export_fig(FileName, '-pdf')
          close
        elseif SaveFigs == 2
          FolderPath = strcat(FigDirPath, '/Rates/');
          [status,msg,msgID] = mkdir(FolderPath);
          FileName   = strcat(FolderPath, 'Level',num2str(iBins),'-QNs.fig');
          savefig(fig,FileName)
          close
        end
        iFigure=iFigure+1;

      end
    
    else

      for iBins = iLevelsVec

        jqnIn = Leveljqn(1:NBins(1),1);
        vqnIn = Levelvqn(1:NBins(1),1);

        for iP = 1:3    
          if PlotPairUser(iP) == 1

            figure(iFigure);
            fig = gcf;
            screensize = get( groot, 'Screensize' );
            fig.Position=screensize;
            fig.Color='None';

            %subplot(2,1,[iP, iP+2])
            iMol = Pair_To_Molecule(iP);

            jqn = Leveljqn(1:NBins(iMol),iMol);
            vqn = Levelvqn(1:NBins(iMol),iMol);
            RatesMatrix_Temp2 = zeros(max(jqn)+1,max(vqn)+1);

            for jBins = 1:NBins(iMol)
                RatesMatrix_Temp2(jqn(jBins)+1,vqn(jBins)+1) = RatesMatrix(iBins,jBins,iP,iTint);
            end

            x = [0:1:max(vqn)];
            y = [0:1:max(jqn)];
            [X,Y] = meshgrid(x,y);
            %contourf(X,Y,log10(RatesMatrix_Temp2),50);
            %surf(X,Y,log10(RatesMatrix_Temp2));
            %az = 60;
            %el = 20;
            %view(az, el);
            pcolor(X,Y,log10(RatesMatrix_Temp2));
            str_y = ['J'];
            str_x = ['v'];
            clear x y X Y 
            cb=colorbar;
            ylim(cb,[-15 log10(5.e-10)]);
            cbstr = num2str(iBins);
            cbstr = ['$log_{10}(K_{',cbstr,',(v,J)}) [cm^3/s]$'];
            ylabel(cb, cbstr);
            set(cb,'FontSize',LegendFontSz,'FontName',LegendFontNm,'TickLabelInterpreter','latex');
            cb.Label.Interpreter = 'latex';
            hold on
            plot(vqn(iBins,1),jqn(iBins,1),'ko','MarkerSize',10,'MarkerFaceColor','k')   

            xlab = xlabel(str_x);
            xlab.Interpreter = 'latex';

            ylab = ylabel(str_y);
            ylab.Interpreter = 'latex';

            str_title = strcat('$',MoleculesName(1,:),'(',num2str(vqnIn(iBins)),',',num2str(jqnIn(iBins)),')',{' '},'->',{' '},MoleculesName(iMol,:),'(v,j)$');
            tit = title(str_title);
            tit.Interpreter = 'latex';
            
            shading(gca,'interp')
            set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
            daspect([1 1 1])
            %pbaspect([1 1 1])

            if SaveFigs == 1
              FolderPath = strcat(FigDirPath, '/Rates/');
              [status,msg,msgID] = mkdir(FolderPath);
              FileName   = strcat(FolderPath, 'Level',num2str(iBins),'-Pair',num2str(iP),'-QNs');
              export_fig(FileName, '-pdf')
              close
            elseif SaveFigs == 2
              FolderPath = strcat(FigDirPath, '/Rates/');
              [status,msg,msgID] = mkdir(FolderPath);
              FileName   = strcat(FolderPath, 'Level',num2str(iBins),'-Pair',num2str(iP),'-QNs.fig');
              savefig(FileName)
              close
            end
            iFigure=iFigure+1;

          end
        end

      end
      
    end
    
  end

end