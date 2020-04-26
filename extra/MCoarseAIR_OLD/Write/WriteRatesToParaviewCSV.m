function [] = WriteRatesToParaviewCSV(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, rOut, Levelg, LevelEeVVib0, LevelEeVRot, VMax, RatesMatrix, ProcessesRates)

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

  
  iBinnedMol = 1;
  if FinalBin < 2
    FinalBin = NLevels(iBinnedMol)
  end
  fprintf('CALLING WriteRatesAsNetwork: Writing Rates for Levels %i-%i in a Network File.\n',StartBin,FinalBin)

  
  cm3_To_m3          = 1.0;%1.d-6;
  MinValA            = 5.d-13 .* cm3_To_m3;         %%% For Cutting out Rates < MinValA
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

  %Kij(:,:)                        = RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,1) .* cm3_To_m3;
  Kij(:,:)                        = (RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,1) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),2,1) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),3,1) ) .* cm3_To_m3;
  %Kij(:,:)                        = (RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),2,1) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),3,1) ) .* cm3_To_m3;
  Kij(Kij < MinValA)              = 0.d0;
  Kji                             = tril(Kij) + tril(Kij.*ExpMat,-1)';
  
  KDiss(1:NLevels(iBinnedMol),1)  = ProcessesRates(1:NLevels(iBinnedMol),1,1)' .* cm3_To_m3;
  KEq(1:NLevels(iBinnedMol),1)    = ( EqVec(1:NLevels(iBinnedMol)) ).^(-1);
  KRecOverg(:,1)                  = KDiss(:,1) ./ KEq(:,1);
  KRec(:,1)                       = KRecOverg(:,1);   
 
  
  
    
  WriteDissRatesFlg      = 1
  
  WriteInelasticRatesFlg = 1
    vLevels  = [0,   0,   4,  10,  15,  19,  21, 24, 46, 40, 30,   5, 20, 12, 15, 20, 18,  5, 25,  10];
    JLevels  = [0, 220, 114,   0, 190, 180, 129, 40, 14, 10, 20, 200, 70, 70, 70,  0, 80, 45, 80, 155];
%     vLevels  = [0, 24,   4, 15, 10,  15,   5,   5, 15, 40, 40];
%     JLevels  = [0, 40, 114,  0, 60, 140, 180, 210, 90,  0, 20];  
%     vLevels  = [5,    5,   5,   5,   5, ...
%                 5,    5,   5,   5,   5, ...
%                 5,    5,  20,  20,  20, ...
%                 20,  20,  20,  20,  20, ...
%                 20,  20,  20,  20,  30, ...
%                 30,  30,  30,  30,  30, ...
%                 30,  30,  30,  30,  30,  30];
%     JLevels  = [  0,  20,  40,  60,  80, ...
%                 100, 120, 140, 160, 180, ...
%                 200, 220,   0,  20,  40, ...
%                  60,  80, 100, 120, 140, ...
%                 160, 180, 200, 220,   0, ...
%                  20,  40,  60,  80, 100, ...
%                 120, 140, 160, 180, 200, 220];  
  
  iBinnedMol         = 1;
  FloderName         = strcat('./',MoleculesName(iBinnedMol,:),'_Rates')
  [status,msg,msgID] = mkdir(FloderName)
  
%   ExpVec(1:NLevels(iBinnedMol),1)  = Levelg(1:NLevels(iBinnedMol),1) .* exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0_Vec(1) .* UKb) );
%   %ExpVec                           = ExpVec ./ sum(ExpVec);
%   
%   iReac = 1;
%   for iComp = 1:NComp
%     QTran(iComp) = (2.d0 .* pi .* ComponentMass(iComp) .* DSWtoKg ./ AvN .* UKb .* T0_Vec(iT) ./ Plnck.^2).^(3/2 .* RxLxIdx(iComp));
%   end
%   IntDeg(:) = ComponentDeg(:).^(RxLxIdx(:));
%   EqVec     = prod(QTran) .* prod(IntDeg) .* ExpVec(1:NLevels(iBinnedMol),1);
%   
%   KDiss     = ProcessesRates(1:NLevels(iBinnedMol),1,iT) .* 1.d-6;
%   Kij(:,:)  = (RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,iT)) .* 1.d-6;
%   %Kij(:,:)  = (RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,iT) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),2,iT) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),3,iT)) .* 1.d-6;
%   %Kij(:,:)  = (RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),2,iT) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),3,iT)) .* 1.d-6;
%   
%   KEq(1:NLevels(iBinnedMol),1)       = ( EqVec(1:NLevels(iBinnedMol)) ).^(-1);
%   KRec(1:NLevels(iBinnedMol),1)      = KDiss(1:NLevels(iBinnedMol),1) ./ KEq(1:NLevels(iBinnedMol),1);
%   KRecOverg(1:NLevels(iBinnedMol),1) = KRec(1:NLevels(iBinnedMol),1)  ./ Levelg(1:NLevels(iBinnedMol),1);
%   
  
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
  
  ii  = 0;%12;
  iii = 0;
  for iLevels = 1:size(vLevels,2)
    for jLevels = 1:NLevels(iBinnedMol)
      if Levelvqn(jLevels) == vLevels(iLevels) && Leveljqn(jLevels) == JLevels(iLevels)
        iii = iii + 1;
        LevelsList(iii) = jLevels;
      end
    end
  end
  
  if WriteInelasticRatesFlg == 1

    for iLevels = LevelsList

      FileName = strcat(FloderName,'/InelRate_v',num2str(Levelvqn(iLevels)),'_j',num2str(Leveljqn(iLevels)),'.csv');
      ii = ii + 1;
      %FileName = strcat(FloderName,'/InelRates.csv.',num2str(ii));
      %FileName = strcat(FloderName,'/InelRates_',num2str(ii),'.csv');
      fileID = fopen(FileName,'w');
      fprintf(fileID,'jLevel,rIn0,J0,EeV0,rIn,rOut,v,J,EeV,LevelEeVVib0,LevelEeVRot,Deltav,Deltaj,DeltaEv,DeltaEj,Kij\n');
      
      %fprintf(fileID,'%i,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n', iLevels, rIn(iLevels,iBinnedMol), rOut(iLevels,iBinnedMol), Levelvqn(iLevels,iBinnedMol), ...
      %                                                              Leveljqn(iLevels,iBinnedMol), LevelEeV(iLevels,iBinnedMol), LevelEeVVib0(iLevels,iBinnedMol), ...
      %                                                              LevelEeVRot(iLevels,iBinnedMol), 0.0, 0.0, ...
      %                                                              log(1.d-11));

      DeltaEj = zeros(max(Leveljqn)+1,1);
      DeltaEv = zeros(max(Levelvqn)+1,1);
      DeltaE  = zeros(NLevels,1);
      for jLevels = 1:NLevels(iBinnedMol)
        if Levelvqn(jLevels) == Levelvqn(iLevels)
          DeltaEj(Leveljqn(jLevels)+1,1) = Levelg(jLevels) ./ Levelg(iLevels) .* exp( - (LevelEeV(jLevels) - LevelEeV(iLevels)) .* Ue ./ (T0_Vec(1) .* UKb) );
        end 
        if Leveljqn(jLevels) == Leveljqn(iLevels)
          DeltaEv(Levelvqn(jLevels)+1,1) = exp( - (LevelEeV(jLevels) - LevelEeV(iLevels)) .* Ue ./ (T0_Vec(1) .* UKb) );
        end 
        DeltaE(jLevels) = LevelEeV(jLevels) - LevelEeV(iLevels);
        Deltav(jLevels) = Levelvqn(jLevels) - Levelvqn(iLevels);
        Deltaj(jLevels) = Leveljqn(jLevels) - Leveljqn(iLevels);
      end
      
      for jLevels = 1:NLevels(iBinnedMol)
        if (Kji(iLevels,jLevels) > MinValA) && (iLevels~=jLevels)%&& DeltaE(jLevels) < 0.0
          fprintf(fileID,'%i,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n', jLevels, rIn(iLevels,iBinnedMol),         Leveljqn(iLevels,iBinnedMol),     LevelEeV(iLevels,iBinnedMol),     ...
                                                                                       rIn(jLevels,iBinnedMol),         rOut(jLevels,iBinnedMol),         Levelvqn(jLevels,iBinnedMol),     ...
                                                                                       Leveljqn(jLevels,iBinnedMol),    LevelEeV(jLevels,iBinnedMol),     LevelEeVVib0(jLevels,iBinnedMol), ...
                                                                                       LevelEeVRot(jLevels,iBinnedMol), Deltav(jLevels), Deltaj(jLevels), DeltaEv(Levelvqn(jLevels)+1),     ...
                                                                                       DeltaEj(Leveljqn(jLevels)+1),    Kji(iLevels,jLevels) );
        end
      end

      fclose(fileID);   

    end
  
  end
  
  
  if WriteInelasticRatesFlg == 2
    
    load('/Users/sventuri/Downloads/Prob_convolution.mat')%FC_matrix.mat    
    for iLevels = LevelsList

      FileName = strcat(FloderName,'/Factor_v',num2str(Levelvqn(iLevels)),'_j',num2str(Leveljqn(iLevels)),'.csv');
      fileID = fopen(FileName,'w');
      fprintf(fileID,'jLevel,rIn,rOut,v,J,EeV,Factor\n');
      
      fprintf(fileID,'%i,%e,%e,%e,%e,%e,%e\n', iLevels, rIn(iLevels,iBinnedMol), rOut(iLevels,iBinnedMol), Levelvqn(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), LevelEeV(iLevels,iBinnedMol), 1.d-8);
      
      Daje(:) = exp( - (LevelEeV(iLevels,1) - LevelEeV(1:NLevels(iBinnedMol),1)) .* Ue ./ (T0_Vec(1) .* UKb) );
      for jLevels = 1:NLevels(iBinnedMol)
        %if Kij(iLevels,jLevels) > MinValA
          fprintf(fileID,'%i,%e,%e,%e,%e,%e,%e\n', jLevels, rIn(jLevels,iBinnedMol), rOut(jLevels,iBinnedMol), Levelvqn(jLevels,iBinnedMol), Leveljqn(jLevels,iBinnedMol), LevelEeV(jLevels,iBinnedMol), Prob_convolution(iLevels,jLevels).*Daje(jLevels) );
        %end
      end

      fclose(fileID);   

    end
  
  end
  
  
end