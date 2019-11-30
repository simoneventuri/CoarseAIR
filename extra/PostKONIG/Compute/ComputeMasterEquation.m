function [] = ComputeMasterEquation(iT, NLevels, LevelEeV, LevelEeV0, Levelvqn, Leveljqn, Levelg, RatesMatrix, ProcessesRates, t, MolFracs, nd, Pop, StpInstants)

  global T0_Vec SaveFigs MoleculesName UKb Ue AvN FigDirPath AxisFontSz AxisFontNm AxisLabelSz AxisLabelNm LegendFontSz LegendFontNm MinEvPlot MaxEvPlot

  LoadFlg = 2

  P0 = 9940.737d0;
  X0 = [0.05d0, 0.95d0];
  T0 = 300;
  V0 = 1.d0;
  
  MinValA    = 1.d-15 .* 1.d-6;

  
  iSteps     = 1;
  iBinnedMol = 1;
  R          = 8.3144598;
  iFigure    = 2001
  
  NewFolder = strcat(FigDirPath,'/KONIG_vs_Analytics')
  [status,msg,msgID] = mkdir(NewFolder)
  
  
  ExpVec0(1:NLevels(iBinnedMol),1)     = exp( - LevelEeV0(1:NLevels(iBinnedMol),1) .* Ue ./ (T0 .* UKb) );
  Q0(1:NLevels(iBinnedMol),iBinnedMol) = Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* ExpVec0(1:NLevels(iBinnedMol),iBinnedMol);
  ExpVec0                              = ExpVec0 ./ sum(ExpVec0);

  ExpVec(1:NLevels(iBinnedMol),1)       = exp( - LevelEeV0(1:NLevels(iBinnedMol),1) .* Ue ./ (T0_Vec(1) .* UKb) );
  Q(1:NLevels(iBinnedMol),iBinnedMol)   = Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* ExpVec(1:NLevels(iBinnedMol),iBinnedMol);
  QEn(1:NLevels(iBinnedMol),iBinnedMol) = Q(1:NLevels(iBinnedMol),iBinnedMol) .* LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol);
  ExpVec                                = ExpVec ./ sum(ExpVec);
  ExpMat                                = kron(ExpVec,1.d0./ExpVec');
  QMat                                  = kron(Q,1.d0./Q');  
  
  
  n0      = (P0 * V0 / (R * T0)) * AvN;
  nO0     = n0 * X0(1);
  nCOTot0 = n0 * X0(2);
  nCO0    = nCOTot0 .* Q0(1:NLevels(iBinnedMol)) ./ sum(Q0) ./ Levelg(1:NLevels(iBinnedMol),iBinnedMol);
  
  
  if LoadFlg == 0
    
    Kij                                              = zeros(NLevels(iBinnedMol),NLevels(iBinnedMol));
    Kij(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol)) = (RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,1) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),2,1)) .* 1.d-6;
    %Kij(:,:)                        = RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,1) .* 1.d-6;
    Kij(Kij < MinValA)              = 0.d0;
    Kji                             = tril(Kij .* ExpMat,-1) + tril(Kij,-1)';%+ diag(diag(Kij));%Kij .* ExpMat;%  
    %DissRates                       = ProcessesRates(1:NLevels(iBinnedMol),1,1)' .* 1.d-6;

    AOrig                           = Kji - diag(diag(Kij));
    ExpA                            = expm_new(AOrig .* nO0 .* 1.d-14);
    clear Kij Kji
  
    filename = strcat('./AInelExch')
    save(filename, 'AOrig', 'ExpA', '-v7.3')
    clear AOrig
    
  else
    
    filename = strcat('./AInelExch')
    load(filename, 'AOrig', 'ExpA')
  
  end 
  clear RatesMatrix 
  
  tOld    = 1.d-14;
  ExpANew = ExpA;
  TempStpInstants = StpInstants;
  %TempStpInstants =[2:10:2000];
  for iSteps = TempStpInstants(1:end)
    iSteps
    
    filename = strcat('./AInelExchExp_Step_',num2str(iSteps),'.mat');
    if LoadFlg < 2 || exist(filename, 'file') ~= 2
      
      tic
      tNew    = t(iSteps)
      ExpANew = real(mpower(ExpANew, tNew ./ tOld));
      nCO     = ExpANew * nCO0;
      tOld    = tNew;
      toc
      
  %     C  = sum(DissRates .* nCO);
  %     nO = nO0 * exp(C * t(iSteps));
  % 
  %     A     = (AOrig - diag(DissRates)) .* nO .* t(iSteps);
  %     ExpA  = expm(A);
  %     nCO   = ExpA * nCO0;
      filename = strcat('./AInelExchExp_Step_',num2str(iSteps))
      save(filename, 'nCO', '-v7.3')
    
    else

      filename = strcat('./AInelExchExp_Step_',num2str(iSteps))
      load(filename, 'nCO')
      
    end

    
    nStSOverG(1:NLevels(iBinnedMol),iSteps) = Pop(iSteps,1:NLevels(iBinnedMol),iBinnedMol)' ./ Levelg(1:NLevels(iBinnedMol),iBinnedMol);
    nStSTot(iSteps)                         = sum( Pop(iSteps,1:NLevels(iBinnedMol),iBinnedMol) );
    eStSTot(iSteps)                         = sum( Pop(iSteps,1:NLevels(iBinnedMol),iBinnedMol)' .* LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) ) ./ nStSTot(iSteps);

    nAnalitic(1:NLevels(iBinnedMol),iSteps) = nCO(1:NLevels(iBinnedMol)) .* Levelg(1:NLevels(iBinnedMol),iBinnedMol);
    nAnaliticOverG(1,1:NLevels(iBinnedMol)) = nCO(1:NLevels(iBinnedMol));
    nAnaliticTot(iSteps)                    = sum( nAnalitic(1:NLevels(iBinnedMol),iSteps) );
    eAnaliticTot(iSteps)                    = sum( nAnalitic(1:NLevels(iBinnedMol),iSteps) .* LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) ) ./ nStSTot(iSteps);

    
    figure(iFigure)
    hold on
    scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),nStSOverG(1:NLevels(iBinnedMol),iSteps),20,'Filled');
    scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),nAnaliticOverG(1:NLevels(iBinnedMol)),20,'Filled');
    hold on
    %scatter(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),nKONIGOverG(1:NLevels(iBinnedMol),iSteps),20,LevToBin(1:NLevels(iBinnedMol),iBinnedMol,1),'Filled');
    %colormap(lines(max(max(max(LevToBin)))))
    %cb=colorbar;
    %cb.Ticks = [1, 1.5]; %Create 8 ticks from zero to 1
    %cb.TickLabels = {'1','2'}
    %ylab = ylabel(cb, '$Group$');
    %ylab.Interpreter = 'latex';
    %set(cb,'FontSize',LegendFontSz,'FontName',LegendFontNm,'TickLabelInterpreter','latex');
    hold on
    xt = get(gca, 'XTick');
    set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
    str_x = ['Energy [eV]'];
    xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
    xlab.Interpreter = 'latex';
    xlim([max(min(LevelEeV(:,iBinnedMol)),MinEvPlot(iBinnedMol)), min(max(LevelEeV(:,iBinnedMol)),MaxEvPlot(iBinnedMol))]); 
    str_y = ['$N_{i} / g_{i} [m^{-3}]$'];
    ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
    ylab.Interpreter = 'latex';
    ylim([1.d5, 1.d23]);
    set(gca, 'YScale', 'log')    
    if SaveFigs == 1
      FolderPath = strcat(FigDirPath, '/KONIG_vs_Analytics/');
      [status,msg,msgID] = mkdir(FolderPath);
      FileName   = strcat(FolderPath, MoleculesName(1,:),'-Pop_',num2str(iSteps),'_t',num2str(t(iSteps),'%7.2e'));
      export_fig(FileName, '-pdf')
      close
    elseif SaveFigs == 2
      FolderPath = strcat(FigDirPath, '/KONIG_vs_Analytics/');
      [status,msg,msgID] = mkdir(FolderPath);
      FileName   = strcat(FolderPath, MoleculesName(1,:),'-Pop_',num2str(iSteps),'_t',num2str(t(iSteps),'%7.2e'),'.fig');
      savefig(FileName)
      %close
    end
    iFigure = iFigure + 1;

  end


  figure(iFigure)
  hold on
  plot(t(TempStpInstants(:)),nStSTot(TempStpInstants(:)),'g');
  hold on
  plot(t(TempStpInstants(:)),nAnaliticTot(TempStpInstants(:)),'r');
  xt = get(gca, 'XTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  yt = get(gca, 'YTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  str_x = ['time [s]'];
  xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  xlab.Interpreter = 'latex';
  set(gca, 'XScale', 'log')
  %xlim([max(min(LevelEeV(:,iBinnedMol)),MinEvPlot(iBinnedMol)), min(max(LevelEeV(:,iBinnedMol)),MaxEvPlot(iBinnedMol))]); 
  str_y = ['$N_{CO}$'];
  ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  ylab.Interpreter = 'latex';
  %ylim([1.d5, 1.d20]);
  %set(gca, 'YScale', 'log')
  if SaveFigs == 1
    FolderPath = strcat(FigDirPath, '/KONIG_vs_Analytics/');
    [status,msg,msgID] = mkdir(FolderPath);
    FileName   = strcat(FolderPath, MoleculesName(1,:),'-Pop_',num2str(iSteps),'_t',num2str(t(iSteps),'%7.2e'));
    export_fig(FileName, '-pdf')
    close
  elseif SaveFigs == 2
    FolderPath = strcat(FigDirPath, '/KONIG_vs_Analytics/');
    [status,msg,msgID] = mkdir(FolderPath);
    FileName   = strcat(FolderPath, MoleculesName(1,:),'-Pop_',num2str(iSteps),'_t',num2str(t(iSteps),'%7.2e'),'.fig');
    savefig(FileName)
    close
  end
  iFigure = iFigure + 1;
  


  figure(iFigure)
  hold on
  plot(t(TempStpInstants(:)),eStSTot(TempStpInstants(:)),'g');
  hold on
  plot(t(TempStpInstants(:)),eAnaliticTot(TempStpInstants(:)),'r');
  xt = get(gca, 'XTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  yt = get(gca, 'YTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  str_x = ['time [s]'];
  xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  xlab.Interpreter = 'latex';
  set(gca, 'XScale', 'log')
  %xlim([max(min(LevelEeV(:,iBinnedMol)),MinEvPlot(iBinnedMol)), min(max(LevelEeV(:,iBinnedMol)),MaxEvPlot(iBinnedMol))]); 
  str_y = ['$Energy [eV]$'];
  ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  ylab.Interpreter = 'latex';
  %ylim([1.d5, 1.d20]);
  %set(gca, 'YScale', 'log')
  if SaveFigs == 1
    FolderPath = strcat(FigDirPath, '/KONIG_vs_Analytics/');
    [status,msg,msgID] = mkdir(FolderPath);
    FileName   = strcat(FolderPath, MoleculesName(1,:),'-Sts_vs_KONIG_Energy');
    export_fig(FileName, '-pdf')
    close
  elseif SaveFigs == 2
    FolderPath = strcat(FigDirPath, '/KONIG_vs_Analytics/');
    [status,msg,msgID] = mkdir(FolderPath);
    FileName   = strcat(FolderPath, MoleculesName(1,:),'-Sts_vs_KONIG_Pop');
    savefig(FileName)
    close
  end
  iFigure = iFigure + 1;


  % 
  % FileFormatName = strcat('/Users/sventuri/WORKSPACE/neqplasma_QCT/database_CO2Binned/thermo/CO_Format')
  % FileFormatID   = fopen(FileFormatName);
  % 
  % FileNewName    = strcat('/Users/sventuri/WORKSPACE/neqplasma_QCT/database_CO2Binned/thermo/CO')
  % fileNewID      = fopen(FileNewName,'w');
  % 
  % tline  = fgetl(FileFormatID);
  % ttline = strcat(tline,'\n');
  % while ischar(tline)
  %   fprintf(fileNewID,ttline);
  %   tline = fgetl(FileFormatID);
  %   ttline = strcat(tline,'\n');
  % end
  % fclose(FileFormatID);
  % 
  % fprintf(fileNewID,'NB_ENERGY_LEVELS = %g\n', NbBins);
  % for iBin = 1:NbBins
  %   fprintf(fileNewID,'%14.8E  %14.8E\n', QBin(iBin), 0.d0);%EeVBin(iBin))
  % end
  % fclose(fileNewID);  
 
end