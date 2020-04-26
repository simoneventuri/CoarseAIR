function [iFigure, PESToMoleFracCorr] = ComputePESToMoleFracCorr(iFigure, iComp, MoleFracOI, tStoch, MolFracsStoch, RGrid , EGrid, EGridMean, EGridSD, TopPercentage, NbSampledPoints, DistThreshold)
   
  global AnglesVecGrid RMinGrid RMaxGrid NPointsGrid
  global iPESStart iPESEnd 
 
  NSigma   = 3.0;
  RRGrid   = RGrid';
  NPESTemp = iPESEnd - iPESStart + 1; 

  for iPES=1:NPESTemp
    iT=1;
    while MolFracsStoch(iT,iComp,iPES+iPESStart-1) > MoleFracOI
      iT=iT+1;
    end
    Deltat    = tStoch(iT) - tStoch(iT-1);
    DeltaMole = MolFracsStoch(iT,iComp,iPES+iPESStart-1) - MolFracsStoch(iT-1,iComp,iPES+iPESStart-1);
    SmallMole = MoleFracOI - MolFracsStoch(iT-1,iComp,iPES+iPESStart-1);
    tOI(iPES) = tStoch(iT-1) + (SmallMole/DeltaMole * Deltat)
  end

  
  A = tOI;
    
  iPoints = 1;
  for Ang=AnglesVecGrid
    File   = strcat('./MoleFrac_To_PES_CorrCoeff.csv.', num2str(Ang));
    fileID = fopen(File,'w');
    fprintf(fileID,'R1,R2,R3,EMean,EMinus,EPlus,CorrCoeff\n');
      
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
        
        Corr                       = corrcoef(A,B);
        PESToMoleFracCorr(iPoints) = Corr(1,2);
        
        
        fprintf(fileID,'%e,%e,%e,%e,%e,%e,%e\n', R1, R2, R3, EMean, EMinus, EPlus, abs(PESToMoleFracCorr(iPoints)) );
        iPoints=iPoints+1;
      end
    end
    
    fclose(fileID);
  end
 
  
  if (NbSampledPoints >0 )
    fprintf('Number of Points to Sample: %i\n', NbSampledPoints);

    %figure(iFigure)
    %histogram(PESToMoleFracCorr);
    %iFigure=iFigure+1;

    %figure(iFigure)
    CDFDist = histogram(abs(PESToMoleFracCorr),300,'Normalization','cdf');
    iBin=1;
    while CDFDist.Values(iBin) < (1.0-TopPercentage/100.0)
      iBin=iBin+1;
    end
    CorrThreshold = (CDFDist.BinEdges(iBin) + CDFDist.BinEdges(iBin-1)) / 2.0;
    %iFigure=iFigure+1;

    NPoints = sum( abs(PESToMoleFracCorr) > CorrThreshold );
    fprintf('Number of Points with Correlation Over Treshold Value %f: %i \n', CorrThreshold, NPoints);

    jPoints = 1;
    for iPoints=1:length(PESToMoleFracCorr)
      if abs(PESToMoleFracCorr(iPoints)) > CorrThreshold
        NewPointsVec(jPoints)    = iPoints;
        NewCorrCoeffVec(jPoints) = abs(PESToMoleFracCorr(iPoints)) - CorrThreshold;
        if jPoints == 1
          NewCorrCDFVec(jPoints)   = NewCorrCoeffVec(jPoints);
        else
          NewCorrCDFVec(jPoints)   = NewCorrCDFVec(jPoints-1) + NewCorrCoeffVec(jPoints);
        end
        jPoints                  = jPoints + 1;
      end
    end 
    NewCorrCDFVec = NewCorrCDFVec ./ max(NewCorrCDFVec);
    %figure(iFigure)
    %plot(NewCorrCDFVec)
    %hold on
    %histogram(NewCorrCoeffVec,300,'Normalization','cdf');
    %iFigure=iFigure+1;


    File2   = strcat('./MoleFrac_To_PES_Corr_RPoints.csv');
    fileID2 = fopen(File2,'w');
    fprintf(fileID2,'R1,R2,R3\n');
    clear iSample Temp iPoints RealPoint RTemp ETemp CorrTemp MinDist TempDist PointsData
    iSample=1; 
    NAttempts=1;
    while iSample <= NbSampledPoints

      Temp    = rand;
      iPoints = 1;
      while (NewCorrCDFVec(iPoints) < Temp)
        iPoints=iPoints+1;
      end
      RealPoint = NewPointsVec(iPoints);

      RTemp(1:3) = RRGrid(1:3,RealPoint);
      ETemp      = EGridMean(RealPoint);
      CorrTemp   = abs(PESToMoleFracCorr(RealPoint));

      MinDist = 1.e10;
      for jSample=1:iSample-1
        TempDist = sqrt( (RTemp(1) - PointsData(jSample,1)).^2 + (RTemp(2) - PointsData(jSample,2)).^2 + (RTemp(3) - PointsData(jSample,3)).^2 );
        if TempDist < MinDist
          MinDist = TempDist;
        end
      end
      fprintf('Number of Tot Attempts: %i; Number of Points Accepted so far: %i\n', NAttempts, iSample-1);
      NAttempts = NAttempts + 1;


      if MinDist > DistThreshold
        PointsData(iSample,:)  = [RTemp, ETemp, CorrTemp];
        iSample                = iSample+1;
        fprintf(fileID2,'%e,%e,%e\n',       RTemp(1), RTemp(2), RTemp(3) );
      end
    end
    %figure(iFigure)
    %histogram(PointsData(:,4),30);
    %iFigure=iFigure+1;
    %fclose(fileID2);

    for iAng=AnglesVecGrid
      File1   = strcat('./MoleFrac_To_PES_Corr_Points.csv.', num2str(iAng));
      fileID1 = fopen(File1,'w');
      fprintf(fileID1,'R1,R2,R3,EMean,CorrCoeff\n');
      for iPoints=1:size(PointsData,1)
        R1       = PointsData(iPoints,1);
        R2       = PointsData(iPoints,2);
        R3       = PointsData(iPoints,3);
        EMean    = PointsData(iPoints,4);
        CorrTemp = PointsData(iPoints,5);
        Ang3   = acos( (R1.^2 + R2.^2 - R3.^2) ./ (2.d0.*R1.*R2) ) .* 180 ./ pi;
        Ang1   = acos( (R2.^2 + R3.^2 - R1.^2) ./ (2.d0.*R2.*R3) ) .* 180 ./ pi;
        Ang2   = acos( (R1.^2 + R3.^2 - R2.^2) ./ (2.d0.*R1.*R3) ) .* 180 ./ pi;
        DeltaMaxAng = 0.05;
        if     ((Ang1 <= iAng + DeltaMaxAng) && (Ang1 >= iAng - DeltaMaxAng))
          fprintf(fileID1,'%f,%f,%f,%f,%f\n', R2,R1,R3,EMean,CorrTemp);
          %fprintf(fileID1,'%f,%f,%f,%f,%f\n', R3,R1,R2,EMean,CorrTemp);
        elseif ((Ang2 <= iAng + DeltaMaxAng) && (Ang2 >= iAng - DeltaMaxAng))
          fprintf(fileID1,'%f,%f,%f,%f,%f\n', R1,R2,R3,EMean,CorrTemp);
          %fprintf(fileID1,'%f,%f,%f,%f,%f\n', R3,R2,R1,EMean,CorrTemp);
        elseif ((Ang3 <= iAng + DeltaMaxAng) && (Ang3 >= iAng - DeltaMaxAng))
          fprintf(fileID1,'%f,%f,%f,%f,%f\n', R1,R3,R2,EMean,CorrTemp);
          %fprintf(fileID1,'%f,%f,%f,%f,%f\n', R2,R3,R1,EMean,CorrTemp);
        end
      end
      fclose(fileID1);
    end
  
  end
    
end