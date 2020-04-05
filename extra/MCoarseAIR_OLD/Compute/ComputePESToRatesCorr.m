function [iFigure] = ComputePESToRatesCorr(iFigure, RatesMatrixStoch, ProcessesRatesStoch, RGrid , EGrid, EMeanGrid, ESDGrid, IniStateVec, FinStateVec, ProcVec)
    
  global AnglesVecGrid RMinGrid RMaxGrid NPointsGrid
  global iPESStart iPESEnd 
 
  NSigma   = 3.0;
  RRGrid   = RGrid';
  NPESTemp = iPESEnd - iPESStart + 1; 
  
  
  iBinnedMol           = 1;
  ExpVec(:,iBinnedMol) = QBins(:,iBinnedMol);
  ExpVec               = ExpVec ./ sum(ExpVec);
  
  ThermalRatesStoch = zeros(NPESTemp,2);
  for iPES=1:NPESTemp
    ProcessesRatesStochTemp(:,1,iPES) = ProcessesRatesStoch(:,1,iPES+iPESStart-1);
    ProcessesRatesStochTemp(:,2,iPES) = ProcessesRatesStoch(:,3,iPES+iPESStart-1) + ProcessesRatesStoch(:,4,iPES+iPESStart-1);
    for iP = [1,2]
      ThermalRatesStoch(iPES,iP) = sum(ProcessesRatesStochTemp(:,iP,iPES) .* ExpVec(:,1));
    end
  end
 
  
  A1 = squeeze(ThermalRatesStoch(:,1));
  A2 = squeeze(ThermalRatesStoch(:,2));
    
  iPoints = 1;
  for Ang=AnglesVecGrid
    File   = strcat('./ThermalRates_To_PES_CorrCoeff.csv.', num2str(Ang));
    fileID = fopen(File,'w');
    fprintf(fileID,'R1,R2,R3,EMean,EMinus,EPlus,KDissCorr,KExchCorr\n');
      
    for iR1=1:NPointsGrid
      for iR3=iR1:NPointsGrid
        R1 = RRGrid(1,iPoints);
        R2 = RRGrid(2,iPoints);
        R3 = RRGrid(3,iPoints);
        
        EMean         = EGridMean(iPoints);
        EMinus        = EGridMean(iPoints) - NSigma.*EGridSD(iPoints); 
        EPlus         = EGridMean(iPoints) + NSigma.*EGridSD(iPoints); 
        B(1:NPESTemp) = squeeze(EGrid(iPoints,iPESStart:iPESEnd));
        %B             = B ./ mean(B);
        
        KDissCorr = corrcoef(A1,B');
        KExchCorr = corrcoef(A2,B');
        
        fprintf(fileID,'%e,%e,%e,%e,%e,%e,%e,%e\n', R1, R2, R3, EMean, EMinus, EPlus, abs(KDissCorr(1,2)), abs(KExchCorr(1,2)) );
        iPoints=iPoints+1;
      end
    end
    
    fclose(fileID);
  end
  
  clear A1 A2 B
  
  
   
  for iCorr=1:length(IniStateVec)
    iIni = IniStateVec(iCorr);
    iFin = FinStateVec(iCorr);
    iP   = ProcVec(iCorr);
    
    A = squeeze(RatesMatrixStoch(iIni,iFin,iP,iPESStart:iPESEnd));
    
    iPoints = 1;
    for Ang=AnglesVecGrid
      File   = strcat('./i=',num2str(iIni),'_j=',num2str(iFin),'_iP=',num2str(iP),'_To_PES_CorrCoeff.csv.', num2str(Ang));
      fileID = fopen(File,'w');
      fprintf(fileID,'R1,R2,R3,EMean,EMinus,EPlus,KCorr\n');

      for iR1=1:NPointsGrid
        for iR3=iR1:NPointsGrid
          R1 = RRGrid(1,iPoints);
          R2 = RRGrid(2,iPoints);
          R3 = RRGrid(3,iPoints);

          EMean  = EGridMean(iPoints);
          EMinus = EGridMean(iPoints) - NSigma.*EGridSD(iPoints); 
          EPlus  = EGridMean(iPoints) + NSigma.*EGridSD(iPoints); 
          B      = squeeze(EGrid(iPoints,:));

          KCorr = corrcoef(A,B);

          fprintf(fileID,'%e,%e,%e,%e,%e,%e,%e\n', R1, R2, R3, EMean, EMinus, EPlus, abs(KCorr) );
          iPoints=iPoints+1;
        end
      end

      fclosef(fileID);
    end
    
  end
  

end