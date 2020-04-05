function iFigure = PlotTemperatures(iT, iFigure, t, Tint, Ttrans)    

    global NBinnedMol BinnedMolToComp BinnedMolName T0_Vec CompNames CompColor linS linST FilePath SaveFigs 

    figure(iFigure)
    for iBinnedMol=1:NBinnedMol
        semilogx(t(:),Tint(:,iBinnedMol),'Color',CompColor(BinnedMolToComp(iBinnedMol),:),'linestyle',linST{iBinnedMol},'LineWidth',5);
        TBinnedMol(iBinnedMol,:)=['T_{int} for ',BinnedMolName(1,:)];
        hold on
    end
    legend(TBinnedMol)
    hTrans=semilogx(t,Ttrans,'LineWidth',5);
    %title([System, ' System Temperatures, T = ', num2str(T0_Vec(iT))]);
    xlabel(['Time [s]']);
    ylabel(['Temperature [K]']);
    set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
    set(gcf, 'PaperPositionMode', 'auto');
    if SaveFigs == 1
        FileName = strcat('IntTemperature-',num2str(T0_Vec(iT)),'K');
        FilePathNamePng = strcat(FilePath, FileName, '.png');
        FilePathNameFig = strcat(FilePath, FileName, '.fig');
        saveas(gcf,strcat(FilePath, FileName),'pdf')
        saveas(gcf,FilePathNamePng)
        savefig(FilePathNameFig);
    end
    iFigure=iFigure+1;
    
end