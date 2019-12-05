function [] = WriteRatesToParaviewCSV(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, rOut, Levelg, RatesMatrix, ProcessesRates)

  global PathToOutput System T0_Vec NBins KinMthd NTtra NTint TInt_Vec NMolecules StartBin FinalBin NBinnedMol BinnedMolName

  global NAtoms AtomsName MoleculesName DegeneracyFactor ColPartToComp BinnedMolToComp NComp CompNames CompColor AtomColor AtomSize AllMoleculesName ...
         PairColor AtomMass ComponentMass ColorVec ComponentDeg MoleculeMu

  global Plnck UKb Ue KeV AvN AMUToKg EhToeV DSWtoKg ATMToPa

  global SystemPath RatesPath DatabasePath RunKONIGDir OutputPath OutputPathSaveFigs FigDirPath linS linST AxisFontSz AxisFontNm AxisLabelSz AxisLabelNm ...
         LegendFontSz LegendFontNm XLimPlot YLimPlot 

  global PlotMolFracsFlg PlotTemperaturesFlg PlotKsFlg PlotKsMoleFracsFlg PlotMovieFlg PlotStepsFlg PlotEnergiesFlg PlotProcProbFlg PlotNewSahaFlg ...
         PlotTausFlg LevToBinFlg

  global NStepsMovie TimeMinMovie TimeMaxMovie NSubplots NInstantsPerSubplot tInstants MinEvPlot MaxEvPlot StepsOverlappingSteps vqnColor SubplotsFlg ...
         QNSpace RxLxIdx ReadAllRatesFlg SaveRatesFlg CompOI PlotPairUser BinnedFlg MovieFlg

  global xDim yDim TotDim PlotPair ProcToLevIP TempChar TempChar2 Pair_Name Pair_To_Molecule Pair_to_Atoms iPInternal iPExternal

  global RCVec BCVec GCVec KCVec OCVec PCVec WCVec JCVec YCVec CCVec MCVec
  
  global WriteRatesFlag WriteSrcTermsFlag WriteDissFlag WriteInelasticFlag WriteAllFlag ConsiderDegFlg
  
  MinValA  = 1.d-14 .* 1.d-6;
  
  WriteDissRatesFlg      = 1
  
  WriteInelasticRatesFlg = 1
%     vLevels  = [0,   0,   4,  10,  15,  19,  21, 24, 46];
%     JLevels  = [0, 285, 114,   0, 190, 230, 129, 40, 14];
    vLevels  = [0, 24,   4, 15, 10,  15,   5,   5, 15, 40, 40];
    JLevels  = [0, 40, 114,  0, 60, 140, 180, 210, 90,  0, 20];
  
  
  iBinnedMol         = 1;
  FloderName         = strcat('./',MoleculesName(iBinnedMol,:),'_Rates')
  [status,msg,msgID] = mkdir(FloderName)
  
  ExpVec(1:NLevels(iBinnedMol),1)  = Levelg(1:NLevels(iBinnedMol),1) .* exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0_Vec(1) .* UKb) );
  %ExpVec                           = ExpVec ./ sum(ExpVec);
  
  iReac = 1;
  for iComp = 1:NComp
    QTran(iComp) = (2.d0 .* pi .* ComponentMass(iComp) .* DSWtoKg ./ AvN .* UKb .* T0_Vec(iT) ./ Plnck.^2).^(3/2 .* RxLxIdx(iComp));
  end
  IntDeg(:) = ComponentDeg(:).^(RxLxIdx(:));
  EqVec     = prod(QTran) .* prod(IntDeg) .* ExpVec(1:NLevels(iBinnedMol),1);
  
  KDiss     = ProcessesRates(1:NLevels(iBinnedMol),1,iT) .* 1.d-6;
  Kij(:,:)  = (RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,iT)) .* 1.d-6;
  %Kij(:,:)  = (RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,iT) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),2,iT) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),3,iT)) .* 1.d-6;
  %Kij(:,:)  = (RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),2,iT) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),3,iT)) .* 1.d-6;
  
  KEq(1:NLevels(iBinnedMol),1)       = ( EqVec(1:NLevels(iBinnedMol)) ).^(-1);
  KRec(1:NLevels(iBinnedMol),1)      = KDiss(1:NLevels(iBinnedMol),1) ./ KEq(1:NLevels(iBinnedMol),1);
  KRecOverg(1:NLevels(iBinnedMol),1) = KRec(1:NLevels(iBinnedMol),1)  ./ Levelg(1:NLevels(iBinnedMol),1);
  
  
  if WriteDissRatesFlg == 1
    
    FileName = strcat(FloderName,'/DissRates.csv');
    fileID = fopen(FileName,'w');
    fprintf(fileID,'Id,rIn,rOut,v,J,EeV,KDissi,KReci,KRecOvergi\n');
    for iLevels = 1:NLevels(iBinnedMol)
      if KDiss(iLevels) > MinValA
        fprintf(fileID,'%i,%e,%e,%e,%e,%e,%e,%e,%e\n', iLevels, rIn(iLevels,iBinnedMol), rOut(iLevels,iBinnedMol), Levelvqn(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), LevelEeV(iLevels,iBinnedMol), KDiss(iLevels), KRec(iLevels), KRecOverg(iLevels) );
      end
    end
    fclose(fileID);   

  end
  
  
  if WriteInelasticRatesFlg == 1

    for iLevels = 1:size(vLevels,2)
      for jLevels = 1:NLevels(iBinnedMol)
        if Levelvqn(jLevels) == vLevels(iLevels) && Leveljqn(jLevels) == JLevels(iLevels)
          LevelsList(iLevels) = jLevels
        end
      end
    end


    for iLevels = LevelsList

      FileName = strcat(FloderName,'/InelRate_v',num2str(Levelvqn(iLevels)),'_j',num2str(Leveljqn(iLevels)),'.csv');
      fileID = fopen(FileName,'w');
      fprintf(fileID,'jLevel,rIn,rOut,v,J,EeV,Kij\n');
      
      fprintf(fileID,'%i,%e,%e,%e,%e,%e,%e\n', iLevels, rIn(iLevels,iBinnedMol), rOut(iLevels,iBinnedMol), Levelvqn(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), LevelEeV(iLevels,iBinnedMol), 1.d-8);

      for jLevels = 1:NLevels(iBinnedMol)
        if Kij(iLevels,jLevels) > MinValA
          fprintf(fileID,'%i,%e,%e,%e,%e,%e,%e\n', jLevels, rIn(jLevels,iBinnedMol), rOut(jLevels,iBinnedMol), Levelvqn(jLevels,iBinnedMol), Leveljqn(jLevels,iBinnedMol), LevelEeV(jLevels,iBinnedMol), Kij(iLevels,jLevels) .* 1.d6);
        end
      end

      fclose(fileID);   

    end
  
  end
  
  
end