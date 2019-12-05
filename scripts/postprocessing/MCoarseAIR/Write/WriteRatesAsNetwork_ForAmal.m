function [] = WriteRatesAsNetwork_ForAmal(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, Levelg, RatesMatrix, ProcessesRates, LevelEeVVib0, LevelEeVRot, DeltaEintDiss)

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
  
  global WriteRatesFlag WriteSrcTermsFlag WriteDissFlag WriteInelasticFlag WriteAllFlag

  iBinnedMol = 1;
  if FinalBin < 2
    FinalBin = NLevels(iBinnedMol)
  end
  fprintf('CALLING WriteRatesAsNetwork: Writing Rates for Levels %i-%i in a Network File.\n',StartBin,FinalBin)

  
  cm3_To_m3          = 1.d-6;
  MinValA            = 1.d-14 .* cm3_To_m3;         %%% For Cutting out Rates < MinValA
  fprintf('  Cutting Value for Rates = %e m3/s\n',MinValA)
  MinValW            = 1.d-4;
  MinValW_Bis        = 1.d4;

  
  %%% Translational Temperature Conditions:
  ExpVec(1:NLevels(iBinnedMol),1)  = Levelg(1:NLevels(iBinnedMol),1) .* exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0_Vec(1) .* UKb) );
  ExpVec                           = ExpVec ./ sum(ExpVec);
  ExpMat                           = kron(ExpVec,1.d0./ExpVec');  
  
  iReac = 1;
  for iComp = 1:NComp
    QTran(iComp) = (2.d0 .* pi .* ComponentMass(iComp) .* DSWtoKg ./ AvN .* UKb .* T0_Vec(iT) ./ Plnck.^2).^(3/2 .* RxLxIdx(iComp));
  end
  IntDeg(:)  = ComponentDeg(:).^(RxLxIdx(:));
  EqVec      = prod(QTran) .* prod(IntDeg) .* Levelg(1:NLevels(iBinnedMol),1) .* ExpVec(1:NLevels(iBinnedMol),1);

  Kij(:,:)                        = (RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,1) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),2,1)) .* 1.d-6;
  Kij(Kij < MinValA)              = 0.d0;
  Kji                             = tril(Kij) + tril(Kij.*ExpMat,-1)';
  
  KDiss(1:NLevels(iBinnedMol),1)  = ProcessesRates(1:NLevels(iBinnedMol),1,1)' .* 1.d-6;
  KEq(1:NLevels(iBinnedMol),1)    = ( EqVec(1:NLevels(iBinnedMol)) ).^(-1);
  KRecOverg(:,1)                  = KDiss(:,1) ./ KEq(:,1);
  KRec(:,1)                       = KRecOverg(:,1);   
 
  
  tic
  FileName = strcat('./DissRatesAmal.net');
  fileID   = fopen(FileName,'w');
  for i = StartBin:FinalBin
    if KDiss(i) > MinValA
      fprintf(fileID,'%i %e %e\n', i-1, KDiss(i), KRec(i) );
    end
  end
  
  
  tic
  FileName = strcat('./InelRatesAmal.net');
  fileID   = fopen(FileName,'w');
  for i = StartBin:FinalBin
    for j = i+1:FinalBin
      if Kji(i,j) > MinValA
        fprintf(fileID,'%i %i %e %e\n', i-1, j-1, Kji(i,j), Kji(j,i));
      end
    end
  end
  fclose(fileID);
  timee = toc;
  fprintf('  Inelastic Rates Matrix written in %e s\n', timee)


  FileName1 = strcat('./LevelsAmal.csv');
  fileID1   = fopen(FileName1,'w');
  for i = StartBin:FinalBin
    fprintf(fileID1,'%i %e %e %i\n', i-1, LevelEeV(i)-LevelEeV(1), Levelg(i), i-1);
  end
  fclose(fileID1);
  
  fprintf('DONE WITH WriteRatesAsNetwork')
  

end