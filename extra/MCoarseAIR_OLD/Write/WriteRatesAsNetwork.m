function [] = WriteRatesAsNetwork(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, Levelg, RatesMatrix, ProcessesRates, LevelEeVVib0, LevelEeVRot, DeltaEintDiss)

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
  MinValA            = 1.d-13 .* cm3_To_m3;         %%% For Cutting out Rates < MinValA
  fprintf('  Cutting Value for Rates = %e m3/s\n',MinValA)
  MinValW            = 1.d-4;
  MinValW_Bis        = 1.d4;
   
  DissEnLimit = 10.0d0; %0.d0;
  
  
  P0    = 9940.737d0;
  X0    = [0.05d0, 0.95d0];
  T0    = 300;

  
  %%% Initial Conditions:
  ExpVec0(1:NLevels(iBinnedMol),1) = exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0 .* UKb) );
  ExpVec0                          = ExpVec0 ./ sum(ExpVec0);
    
  %%% Translational Temperature Conditions:
  ExpVec(1:NLevels(iBinnedMol),1)  = exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0_Vec(1) .* UKb) );
  ExpVec                           = ExpVec ./ sum(ExpVec);
  ExpMat                           = kron(ExpVec,1.d0./ExpVec');  
  
  
  iReac = 1;
  for iComp = 1:NComp
    QTran(iComp) = (2.d0 .* pi .* ComponentMass(iComp) .* DSWtoKg ./ AvN .* UKb .* T0_Vec(iT) ./ Plnck.^2).^(3/2 .* RxLxIdx(iComp));
  end
  IntDeg(:)  = ComponentDeg(:).^(RxLxIdx(:));
  EqVec      = prod(QTran) .* prod(IntDeg) .* Levelg(1:NLevels(iBinnedMol),1) .* ExpVec(1:NLevels(iBinnedMol),1);

  Kij(:,:)                        = (RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,1) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),2,1)) .* 1.d-6;
  %Kij(:,:)                        = RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,1) .* 1.d-6;
  Kij(Kij < MinValA)              = 0.d0;
  Kji                             = tril(Kij .* ExpMat,-1) + tril(Kij,-1)';%+ diag(diag(Kij));%Kij .* ExpMat;%
  wOut                            = sum(Kji,1);
  wIn                             = sum(Kji,2);
  
  KDiss(1:NLevels(iBinnedMol),1)  = ProcessesRates(1:NLevels(iBinnedMol),1,1)' .* 1.d-6;
  KEq(1:NLevels(iBinnedMol),1)    = ( EqVec(1:NLevels(iBinnedMol)) ).^(-1);
  KRecOverg(:,1)                  = KDiss(:,1) ./ KEq(:,1);
  KRec(:,1)                       = KRecOverg(:,1) .* Levelg(:,1);   
      

  V0      = 1.d0;
  R       = 8.3144598;
  n0      = (P0 * V0 / (R * T0)) * AvN;
  nO0     = n0 * X0(1);
  nCOTot0 = n0 * X0(2);

  Temp = (Kji - diag(wOut)) .* nO0;      
  MinValA = 1.d-200;
  NEdges  = sum(sum(abs(Temp) > MinValA));

  
  NLevelTemp = 0;
  LevelToNew = [];
  for i = StartBin:FinalBin
    if (DeltaEintDiss(i) < DissEnLimit)
      NLevelTemp    = NLevelTemp + 1;
      LevelToNew(i) = NLevelTemp;
    end
  end
  
  
  tic
  FileName = strcat('./InelRates.net');
  fileID   = fopen(FileName,'w');
  fprintf(fileID,'# A network in Pajeks .net format\n');

  fprintf(fileID,'*Vertices %i\n', NLevelTemp);
  for i = StartBin:FinalBin
    if (DeltaEintDiss(i) < DissEnLimit)
      StrNb = ['"v=', num2str(Levelvqn(i)), ',J=', num2str(Leveljqn(i)), '"'];
      fprintf(fileID,'%i %s\n', LevelToNew(i), StrNb);
    end
  end
  

  fprintf(fileID,'*Arcs %i\n',NEdges);
  for i = StartBin:FinalBin
    for j = StartBin:FinalBin
      if (DeltaEintDiss(i) < DissEnLimit) && (DeltaEintDiss(j) < DissEnLimit)
      %if i == j
        if abs(Temp(j,i)) >= MinValA
          %fprintf(fileID,'%i,%i,%e,"direct"\n', j, i, Temp(i,j) );
          fprintf(fileID,'%i %i %e\n', LevelToNew(i), LevelToNew(j), abs(Temp(j,i)) );
        end
      end
    end
  end
  fclose(fileID);
  clear AInel;
  timee = toc;
  fprintf('  Inelastic Rates Matrix written in %e s\n', timee)


  FileName1 = strcat('./LevelsForInelastic.csv');
  fileID1   = fopen(FileName1,'w');
  FileName2 = strcat('./LevelsForDissociation.csv');
  fileID2   = fopen(FileName2,'w');
  fprintf(fileID1,'id,v,Longitude,Latitude,rIn,EeVVib,EeVRot,DeltaEintDiss,OriginalLevel\n');
  fprintf(fileID2,'id,v,Longitude,Latitude,rIn,EeVVib,EeVRot,DeltaEintDiss,OriginalLevel\n');
  for i = StartBin:FinalBin
    if (DeltaEintDiss(i) < DissEnLimit)
      fprintf(fileID1,'%i,%i,%e,%e,%e,%e,%e,%e,%i\n', LevelToNew(i), Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), LevelEeVVib0(i), LevelEeVRot(i), DeltaEintDiss(i), i);
    else
      fprintf(fileID2,'%i,%i,%e,%e,%e,%e,%e,%e,%i\n', i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), LevelEeVVib0(i), LevelEeVRot(i), DeltaEintDiss(i), i);
    end
  end
  fclose(fileID1);
  fclose(fileID2);  
  
  fprintf('DONE WITH WriteRatesAsNetwork')

end