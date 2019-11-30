function [iFigure, Rcorr] = PlotIntPopVsTime(iT, iFigure, t, MolFracs, NLevels, LevelEeV, LevelQ, Levelg, LevelToBin, Pop, QBins, Steps)    

  global NComp NBinnedMol BinnedMolToComp BinnedMolName T0_Vec CompNames CompColor linS linST FilePath SaveFigs 

  for iBinnedMol=1:NBinnedMol

    for iSteps=Steps%1:size(Pop,1)%Steps
      LevelPop(1:NLevels(iBinnedMol),iBinnedMol,iSteps) = Pop(iSteps,LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol)' ./ Pop(1,LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol)';
    end

  end

  for iBinnedMol=1%:NBinnedMol
    
    figure(iFigure)
    A=[];
    for iLevels=1:NLevels(iBinnedMol)
      temp(:)=LevelPop(iLevels,iBinnedMol,Steps(:));
      A = [A, temp(:)];
      semilogx(t(Steps),temp','LineWidth',2);
      hold on
    end    
    Rcorr = corrcoef(A);
    xlim([min(LevelEeV(:,iBinnedMol)), max(LevelEeV(:,iBinnedMol))]);
    %ylim([1.d3, 1.d25]);
    title([BinnedMolName(iBinnedMol,:)]);
    xlabel(['time [s]']);
    ylabel(['(N_{i}(t) / N_{i}(0))']);
    set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
    set(gcf, 'PaperPositionMode', 'auto');
    iFigure=iFigure+1;
    
%     figure(iFigure)
%     for it = 1:length(temp)
%       ttemp(1)=LevelPop(1000,iBinnedMol,Steps(it));
%       ttemp(2)=LevelPop(2000,iBinnedMol,Steps(it));
%       loglog(ttemp(1),ttemp(2),'o')
%       hold on
%       pause(0.1)
%     end    
%     xlim([1.d-5, 1]);
%     ylim([1.d-5, 1]);
%     title([BinnedMolName(iBinnedMol,:)]);
%     xlabel(['N_{10}(t) / N_{10}(0)']);
%     ylabel(['N_{20}(t) / N_{20}(0)']);
%     set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
%     set(gcf, 'PaperPositionMode', 'auto');
%     iFigure=iFigure+1;
    
    figure(iFigure)
    for iLevels=1:NLevels(iBinnedMol)
      temp(:)=LevelPop(iLevels,iBinnedMol,Steps(:));
      for i = 2:length(temp)-1 
        dtemp(i-1) = ( temp(i+1) - temp(i-1) ) ./ ( t(Steps(i+1)) - t(Steps(i-1)) );
      end
      semilogx(t(Steps(2:end-1)),dtemp(:)','LineWidth',2);
      hold on
    end    
    xlim([min(LevelEeV(:,iBinnedMol)), max(LevelEeV(:,iBinnedMol))]);
    %ylim([1.d3, 1.d25]);
    title([BinnedMolName(iBinnedMol,:)]);
    xlabel(['time [s]']);
    ylabel(['dN_{i} / dt']);
    set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
    set(gcf, 'PaperPositionMode', 'auto');
    iFigure=iFigure+1;
    
  end

end