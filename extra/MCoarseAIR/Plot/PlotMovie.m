function iFigure = PlotMovie(iT, iFigure, t, MolFracs, ProcessesRatesOverall, NLevels, LevelEeV, DeltaEintDiss, LevelQ, Levelg, Levelvqn, Leveljqn, LevelToBin, Pop, QBins, Steps, LevToBin, ProcessesRates)    
  
  global QNSpace SubplotsFlg MinEvPlot MaxEvPlot vqnColor 
  global NComp NBinnedMol BinnedMolToComp BinnedMolName T0_Vec CompNames CompColor linS linST FilePath SaveFigs Pair_To_Molecule NTint
    
  ProcessesRatesTemp      = ProcessesRates;
  
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
  
  if QNSpace == 0
    
    %cmap = colormap(lines(max(max(LevToBin))));
   
    figure(iFigure)
    pause(15)

    ii=0;
    for iSteps=Steps  

      ii=ii+1;
      for iBinnedMol=1:NBinnedMol
          LevelPop(1:NLevels(iBinnedMol),iBinnedMol) = Pop(iSteps,LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol)' ./ QBins(LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol);
          %LevelPop(1:NLevels(iBinnedMol),iBinnedMol) = LevelQ(1:NLevels(iBinnedMol),iBinnedMol) ./ Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* Pop(iSteps,LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol)' ./ QBins(LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol);
      end
      
      if SubplotsFlg == 1

        for iBinnedMol=1:NBinnedMol
            subplot(3,(NBinnedMol+1),[1,2+NBinnedMol,5+NBinnedMol]+(iBinnedMol-1))
            semilogy(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),LevelPop(1:NLevels(iBinnedMol),iBinnedMol),'o','Color',CompColor(BinnedMolToComp(iBinnedMol),:),'Markersize',3)
            xlim([min(LevelEeV(:,iBinnedMol)), max(LevelEeV(:,iBinnedMol))]);
            ylim([1.d5, 1.d25]);
            title([BinnedMolName(iBinnedMol,:)]);
            xlabel(['Energy [eV]']);
            ylabel(['N_{i} / g_{i}']);
            set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
            set(gcf, 'PaperPositionMode', 'auto');
        end

        subplot(3,(NBinnedMol+1),(NBinnedMol+1))
        for iComp=1:NComp
            semilogx(t(1:iSteps),MolFracs(1:iSteps,iComp),'Color',CompColor(iComp,:),'linestyle',linS{iComp},'LineWidth',5)
            hold on
            grid on
        end
        legend(CompNames,'location','West')
        set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
        set(gcf, 'PaperPositionMode', 'auto');
        title(['Time = ', num2str(t(iSteps)), ' s']);
        xlabel(['Time [s]']);
        ylabel(['Mole Fraction']);
        xlim([t(Steps(1)), t(Steps(end))]);
        ylim([0, 1]);

        iBinnedMol=1;
        subplot(3,(NBinnedMol+1),(NBinnedMol+1)*2)
        figure(iFigure)
        loglog(t(1:iSteps),ProcessesRatesOverall(1:iSteps,1,iT),'Color',CompColor(BinnedMolToComp(iBinnedMol),:),'LineWidth',5);
        hold on
        grid on
        %title([System, ' System Temperatures, T = ', num2str(T0_Vec(iT))]);
        xlabel(['Time [s]']);
        h=ylabel(['$\bar{K}_{Diss}$']);
        set(h,'interpreter','Latex','FontSize',18)
        set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
        set(gcf, 'PaperPositionMode', 'auto');
        xlim([t(Steps(1)), t(Steps(end))]);
        %ylim([ProcessesRatesOverall(Steps(1),1,iT)*0.9, ProcessesRatesOverall(Steps(end),1,iT)*3]);
        ylim([1.d-15, 1.d-12]);
        
        iBinnedMol=1;
        subplot(3,(NBinnedMol+1),(NBinnedMol+1)*3)
        figure(iFigure)
        loglog(t(1:iSteps),ProcessesRatesOverall(1:iSteps,4,iT),'Color',CompColor(BinnedMolToComp(iBinnedMol),:),'LineWidth',5);
        hold on
        grid on
        %title([System, ' System Temperatures, T = ', num2str(T0_Vec(iT))]);
        xlabel(['Time [s]']);
        h=ylabel(['$\bar{K}_{Exch}$']);
        set(h,'interpreter','Latex','FontSize',18)
        set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
        set(gcf, 'PaperPositionMode', 'auto');
        xlim([t(Steps(1)), t(Steps(end))]);
        %ylim([ProcessesRatesOverall(Steps(1),4,iT)*0.9, ProcessesRatesOverall(Steps(end),4,iT)*3]);
        ylim([1.d-15, 1.d-12]);
        
      else
        
        for iBinnedMol=1%:NBinnedMol
          %subplot(1,(NBinnedMol),iBinnedMol)
          %semilogy(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),LevelPop(1:NLevels(iBinnedMol),iBinnedMol),'o','Color',CompColor(BinnedMolToComp(iBinnedMol),:),'Markersize',3)
          if vqnColor == 0
%             scatter3(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),Leveljqn(1:NLevels(iBinnedMol),iBinnedMol),LevelPop(1:NLevels(iBinnedMol),iBinnedMol),15,LevToBin(1:NLevels(iBinnedMol)),'Filled','o');
%             set(gca, 'ZScale', 'log')
%             xlim([max(min(LevelEeV(:,iBinnedMol)),MinEvPlot(iBinnedMol)), min(max(LevelEeV(:,iBinnedMol)),MaxEvPlot(iBinnedMol))]);
%             zlim([1.d5, 1.d25]);
%             title([BinnedMolName(iBinnedMol,:)]);
%             xlabel(['Energy [eV]']);
%             ylabel(['J']);
%             zlabel(['N_{i} / g_{i}']);
%             grid on
%             view(55,15);
            
            scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),LevelPop(1:NLevels(iBinnedMol),iBinnedMol),15,LevToBin(1:NLevels(iBinnedMol),iBinnedMol),'Filled','o');
            %scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),LevelPop(1:NLevels(iBinnedMol),iBinnedMol),15,'Filled','o');
            %colormap(lines(max(max(LevToBin,1))));
            colormap(distinguishable_colors(max(max(max(LevToBin)))))
            %h=colorbar;
            set(gca, 'YScale', 'log')
            xlim([max(min(LevelEeV(:,iBinnedMol)),MinEvPlot(iBinnedMol)), min(max(LevelEeV(:,iBinnedMol)),MaxEvPlot(iBinnedMol))]);
            ylim([1.d5, 1.d25]);
            title([BinnedMolName(iBinnedMol,:)]);
            xlabel(['Energy [eV]']);
            ylabel(['N_{i} / g_{i}']);
            %caxis([-9.5 -8.5]);
            
          else
            %scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),LevelPop(1:NLevels(iBinnedMol),iBinnedMol),20,DeltaEintDiss(1:NLevels(iBinnedMol),iBinnedMol),'Filled');
            scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),LevelPop(1:NLevels(iBinnedMol),iBinnedMol),20,Levelvqn(1:NLevels(iBinnedMol),iBinnedMol)+1,'Filled');
            %colormap(lines(max(max(Levelvqn))+1));
            colormap(distinguishable_colors(max(max(Levelvqn))+1))
            cb=colorbar;
            %ylim(cb,[-8, 0]); caxis([-8 0]);
            %colormap(jet)
            set(gca, 'YScale', 'log')
            xlim([max(min(LevelEeV(:,iBinnedMol)),MinEvPlot(iBinnedMol)), min(max(LevelEeV(:,iBinnedMol)),MaxEvPlot(iBinnedMol))]);
            ylim([1.d0, 1.d23]);
            title([BinnedMolName(iBinnedMol,:)]);
            xlabel(['Energy [eV]']);
            ylabel(['N_{i} / g_{i}']);
          end
          set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
          set(gcf, 'PaperPositionMode', 'auto');
        end
        title(['Time = ', num2str(t(iSteps)), ' s']);  
          
      end

      hold off
      pause(0.1)

    end
    iFigure=iFigure+1;
    
  else
    
     for iBinnedMol=1:NBinnedMol

      QTot(iBinnedMol)         = sum( Levelg(:,iBinnedMol)  .* exp( - LevelEeV(:,iBinnedMol) .* Ue ./ (T0_Vec(iT) .* UKb) ) )
      LevelPopEq(:,iBinnedMol) = exp( - LevelEeV(:,iBinnedMol) .* Ue ./ (T0_Vec(iT) .* UKb) ) ./ QTot(iBinnedMol);

      figure(iFigure)
      pause(10)

      x     = [0:1:max(Levelvqn(1:NLevels(iBinnedMol),iBinnedMol))];
      y     = [0:1:max(Leveljqn(1:NLevels(iBinnedMol),iBinnedMol))];
      [X,Y] = meshgrid(x,y);

      ii=0;
      for iSteps=Steps

        ii=ii+1;
        LevelPop(1:NLevels(iBinnedMol),iBinnedMol) = LevelQ(1:NLevels(iBinnedMol),iBinnedMol) ./ Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* ...
            Pop(iSteps,LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol)' ./ max(1.d-100,QBins(LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol));  

        LevelPopMatrix = ones(max(Leveljqn(1:NLevels(iBinnedMol),iBinnedMol))+1,max(Levelvqn(1:NLevels(iBinnedMol),iBinnedMol))+1) .* 1.d0;
        for iLevels = 1:NLevels(iBinnedMol)
            LevelPopMatrix(Leveljqn(iLevels,iBinnedMol)+1,Levelvqn(iLevels,iBinnedMol)+1) = LevelPop(iLevels,iBinnedMol) ./ ( nd(iSteps) .* MolFracs(iSteps,3) .* LevelPopEq(iLevels,iBinnedMol) );
        end

%         surf(X,Y,log10(LevelPopMatrix))
%         zlim([8, 25]);
%         title([BinnedMolName(iBinnedMol,:)]);
%         xlabel(['v']);
%         ylabel(['J']);
%         zlabel(['N_{i} / g_{i}']);
%         h=colorbar;
%         caxis([8 25])
%         set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
%         set(gcf, 'PaperPositionMode', 'auto');
%         az = 85;
%         el = 20;
%         view(az, el);

        b = pcolor(log10(LevelPopMatrix));
        %b = contourf(log10(LevelPopMatrix),20);
        cmap = colormap(jet);
        caxis([-6 0.3])
        xlim([0, 20]);
        ylim([0, 150]);
        %b = contourf(log10(LevelPopMatrix),100);
%         for k = 1:length(b)
%           zdata = b(k).ZData;
%           b(k).CData = zdata;
%           b(k).FaceColor = 'interp';
%         end
        %zlim([8, 25]);
        title([BinnedMolName(iBinnedMol,:), ' Time = ', num2str(t(iSteps)), ' s']);
        xlabel(['v']);
        ylabel(['J']);
%         zlabel(['log_{10}(N_{i} / g_{i})']);
        h=colorbar;
        %caxis([10 22])
        set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
        set(gcf, 'PaperPositionMode', 'auto');
%           az = 70;
%           el = 10;
%           view(az, el);

        hold off
        pause(1)
          
      end

      iFigure=iFigure+1;
      
      LevelgMatrix = ones(max(Leveljqn(1:NLevels(iBinnedMol),iBinnedMol))+1,max(Levelvqn(1:NLevels(iBinnedMol),iBinnedMol))+1) .* 1.d0;
      for iLevels = 1:NLevels(iBinnedMol)
          LevelgMatrix(Leveljqn(iLevels,iBinnedMol)+1,Levelvqn(iLevels,iBinnedMol)+1) = Levelg(iLevels,iBinnedMol);
      end
      figure(iFigure)
      b = pcolor(LevelgMatrix);
      iFigure=iFigure+1;     

     end
    
  end

end