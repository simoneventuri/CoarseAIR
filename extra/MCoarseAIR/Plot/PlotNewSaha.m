function [iFigure] = PlotNewSaha(iT, iFigure, iFigureStart, StpInstants, NLevels, Levelg, LevelEeV, DeltaEintDiss, LevToBin, BinEnergy, BinEnergyMax, t, MolFracs, nd, VMax, ProcessesRates, RatesMatrix, Tint)
  
  global Plnck UKb Ue KeV AvN AMUToKg ComponentMass T0_Vec RxLxIdx StepsOverlappingSteps MinEvPlot MaxEvPlot vqnColor NComp ComponentDeg DSWtoKg
 
  %RxLxIdx = 1 if Component is at the LHS; -1 if Component is at the RHS; 0 if Component is not present

  iBinnedMol = 1;

  
  KDiss                         = zeros(NLevels(iBinnedMol),1);
  ExpVec                        = zeros(NLevels(iBinnedMol),1);
  
  KDiss(1:NLevels(iBinnedMol))  = ProcessesRates(1:NLevels(iBinnedMol),1,1);
  ExpVec(1:NLevels(iBinnedMol)) = exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0_Vec(1) .* UKb) );
% 
% %   A      = RatesMatrix(:,:,1) - diag(diag(RatesMatrix(:,:,1)));
% %   B      = sum(A,2);
%   ExpMat = kron(ExpVec,1.d0./ExpVec');  
%   A      = RatesMatrix(:,:,1) - diag(diag(RatesMatrix(:,:,1)));
%   %A(A < 1.d-13) = 0.d0;
%   A      = A .* ExpMat - diag(sum(A,2)) - diag(KDiss);
%   tic
%   [L,D,R]  = eig(A);
%   toc
%   disp('Eig done')  
%   FileInelDissEigs = strcat('./FileInelDissEigs')
%   save(FileInelDissEigs, 'A', 'D', 'R', 'L', '-v7.3')
%   DInv   = diag(1.d0./diag(D));
%   tic
%   AInv = L * DInv * R;
%   toc
%   disp('AInv done')  
%   FileInelDissEigs = strcat('./FileInelDissEigs')
%   save(FileInelDissEigs, 'AInv', 'D', 'DInv', 'R', 'L', 'KDiss', '-v7.3')
  
  
  %cmap     = colormap(lines(max(LevToBin)));
  
  
  iFigure = iFigureStart;
  if StepsOverlappingSteps == 1
    figure(iFigure)
    hold on
    LegendText = [];
  end 
  
  
  for iSteps=StpInstants
    
    if StepsOverlappingSteps == 0
      figure(iFigure)
      hold on
      title(['Time = ', num2str(t(iSteps)), ' s']);
      hold on
    end
    
    iReac = 1;
    for iComp = 1:NComp
      QTran(iComp) = (2.d0 .* pi .* ComponentMass(iComp) .* DSWtoKg ./ AvN .* UKb .* T0_Vec(iT) ./ Plnck.^2).^(3/2 .* RxLxIdx(iComp));
      if RxLxIdx(iComp) == -1
        iSteps
        PopReac(iReac) = MolFracs(iSteps,iComp) .* nd(iSteps);
        iReac = iReac + 1;
      end
    end
    IntDeg(:)  = ComponentDeg(:).^(RxLxIdx(:));
    
    T0                 = T0_Vec(iT);
    Q0(iBinnedMol)     = sum( Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* exp( - LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) .* Ue ./ (T0 .* UKb) ) );
    Q(iBinnedMol)      = sum( Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* ExpVec );
    QDelta(iBinnedMol) = sum( Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* exp( - LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) .* Ue ./ (T0 .* UKb) ));

    %KDissTemp = KDiss + B;
    a   = - KDiss .* 1.d-6 .* PopReac(2);
    b   =   KDiss .* 1.d-6 .* PopReac(2) .* prod(PopReac) .* prod(QTran) .* prod(IntDeg) .* ExpVec;
    y0  =   MolFracs(1,3) .* nd(1) ./ Q0(iBinnedMol) .* exp( - LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) .* Ue ./ (T0 .* UKb) );
    %y0  =   MolFracs(1,3) .* nd(1) ./ Q(iBinnedMol) .* ExpVec;
    c   =   b ./ a + y0;
    Exp = exp( a .* t(iSteps) );
    y   = c .* Exp - b ./ a;
    
%      b = zeros(NLevels(iBinnedMol),1);
%     q = zeros(NLevels(iBinnedMol),1);
%     c = zeros(NLevels(iBinnedMol),1);
%     z = zeros(NLevels(iBinnedMol),NLevels(iBinnedMol));
%     y = zeros(NLevels(iBinnedMol),1);
%   
%     tic
%     b = KDiss .* prod(PopReac) .* prod(QTran) .* prod(IntDeg) .* ExpVec;
%     toc
%     disp('b done')
%     tic
%     y = - A\b;
%     toc
%     disp('y done')
%     tic
%     b = KDiss .* PopReac(2) .* prod(PopReac) .* prod(QTran) .* prod(IntDeg) .* ExpVec;
%     toc
%     tic
%     q = PopReac(2) .* AInv * b;
%     toc
%     disp('q done')
%     tic
%     c = MolFracs(1,3) .* nd(1) ./ Q(iBinnedMol) .* ExpVec + q;
%     toc
%     disp('c done')
%     tic
%     z = L * expm( D .* 1.d-6 .* MolFracs(iSteps,2) .* nd(iSteps) .* t(iSteps) ) * R;
%     toc
%     disp('z done')
%     tic
%     y = z * c - q;   
%     toc
%     disp('y done')  
    
    
%     KEq(1:NLevels(iBinnedMol)) = ( prod(QTran) .* prod(IntDeg) ).^(-1) .* exp( LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) .* Ue ./ (T0_Vec(iT) .* UKb) );
%     KRec(:) = KDiss(:) ./ KEq(:);

        
    LevelPopBoltz(1:NLevels(iBinnedMol),iBinnedMol) = MolFracs(iSteps,3) .* nd(iSteps) ./ Q(iBinnedMol) .* exp( - LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) .* Ue ./ (T0_Vec(iT) .* UKb) );
    LevelPop(1:NLevels(iBinnedMol),iBinnedMol)      = prod(QTran) .* prod(PopReac) .* prod(IntDeg) .* exp( - DeltaEintDiss(1:NLevels(iBinnedMol),iBinnedMol) .* Ue ./ (T0_Vec(iT) .* UKb) );
    LegendPlot(iSteps)                              = LevelPop(1,iBinnedMol);

     
    if StepsOverlappingSteps == 0  
      if vqnColor == 0
        %scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),LevelPop(1:NLevels(iBinnedMol),iBinnedMol),20,'Filled');
        hold on
        %scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),LevelPopBoltz(1:NLevels(iBinnedMol),iBinnedMol),20,'Filled');
        scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),y(1:NLevels(iBinnedMol)),20,log10(KDiss(1:NLevels(iBinnedMol))),'Filled');
        h=colorbar;
        ylim(h,[-10.5, -8]); caxis([-10.5, -8]);
      else
        scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),LevelPop(1:NLevels(iBinnedMol),iBinnedMol),20,KDiss(1:NLevels(iBinnedMol)),'Filled');
        colormap(lines(max(max(Levelvqn))+1));
        h=colorbar;
      end
    else 
      scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),LevelPop(1:NLevels(iBinnedMol),iBinnedMol),20,'Filled');
    end 
    hold on
    xlim([max(min(LevelEeV(:,iBinnedMol)),MinEvPlot(iBinnedMol)), min(max(LevelEeV(:,iBinnedMol)),MaxEvPlot(iBinnedMol))]);
    %ylim([1.d-10, 1.d22]);
    xlabel(['Energy [eV]']);
    ylabel(['N_{i} / g_{i}']);
    set(gca, 'YScale', 'log')
    set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
    set(gcf, 'PaperPositionMode', 'auto');

    if StepsOverlappingSteps == 0
      iFigure=iFigure+1;
    else
      LegendText = [LegendText; strcat('Time = ',num2str(t(iSteps),'%7.2e\n'),' s')]
    end
    
  end
  
  if StepsOverlappingSteps == 1
    scatter(0.d0.*LegendPlot+LevelEeV(1,iBinnedMol),LegendPlot,20,'Filled');
    legend(LegendText);
    iFigure=iFigure+1;
  end
  
end