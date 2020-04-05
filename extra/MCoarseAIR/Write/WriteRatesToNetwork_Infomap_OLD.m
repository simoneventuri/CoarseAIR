function [] = WriteRatesToNetwork_Infomap(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, Levelg, RatesMatrix, ProcessesRates, LevelEeVVib0, LevelEeVRot, DeltaEintDiss)

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
  
  global ConsiderDegFlg WriteExchangeFlg WriteOrigInvFlg WritePajekFlg WriteDirUndFlg
  
  
  DeltaEMax    = -1.3d0; % 100.0 <---
  DeltaERanges = [DeltaEMax, -1.0, -0.8, -0.6, -0.4, -0.3, -0.2, -0.1, 0.0]
  NJGroups     = [0, 100, 300]
  
  MinValA            = 1.d-14 .* 1.d-6;

  iBinnedMol         = 1;  
  FolderName         = strcat(MoleculesName(iBinnedMol,:),'_Network')
  [status,msg,msgID] = mkdir(FolderName)
  
  
  ExpVec0   = zeros(NLevels(iBinnedMol),1);
  ExpVec    = zeros(NLevels(iBinnedMol),1);
  ExpMat    = zeros(NLevels(iBinnedMol),1);
  EqVec     = zeros(NLevels(iBinnedMol),1);
  LogEqVec  = zeros(NLevels(iBinnedMol),1);
  KDiss     = zeros(NLevels(iBinnedMol),1);
  KEq       = zeros(NLevels(iBinnedMol),1);
  KRecOverg = zeros(NLevels(iBinnedMol),1);
  KRec      = zeros(NLevels(iBinnedMol),1);
  Kij       = zeros(NLevels(iBinnedMol),NLevels(iBinnedMol));
  Kij       = zeros(NLevels(iBinnedMol),NLevels(iBinnedMol));
  BinDiss   = zeros(NLevels(iBinnedMol),1);
  
  
  LevelsToWrite    = [];
  LevelsToNOTWrite = [];
  for iLevel=1:NLevels(iBinnedMol) % 1000 <---
    if DeltaEintDiss(iLevel) < DeltaEMax
       LevelsToWrite    =  [LevelsToWrite; iLevel];
    else
       LevelsToNOTWrite = [LevelsToNOTWrite; iLevel];
       iJ=length(NJGroups);
       while (Leveljqn(iLevel) < NJGroups(iJ))
         iJ=iJ-1;
       end
       iE=length(DeltaERanges);
       while (DeltaEintDiss(iLevel) < DeltaERanges(iE))
         iE=iE-1;
       end
       BinDiss(iLevel)  = (iE-1)*(length(NJGroups)-1) + iJ;
    end
  end
  
  
  T0 = 1500;
  ExpVec0(LevelsToWrite(:),1) = exp( - LevelEeV(LevelsToWrite(:),1) .* Ue ./ (T0 .* UKb) );
  ExpVec0                     = ExpVec0 ./ sum(ExpVec0);
   
  TT = T0_Vec(1)
  if (ConsiderDegFlg == 1)
    ExpVec(LevelsToWrite(:),1)  = Levelg(LevelsToWrite(:),1) .* exp( - LevelEeV(LevelsToWrite(:),1) .* Ue ./ (TT .* UKb) );
  else
    ExpVec(LevelsToWrite(:),1)  = exp( - LevelEeV(LevelsToWrite(:),1) .* Ue ./ (TT .* UKb) );
  end
  ExpVec                           = ExpVec ./ sum(ExpVec);
  ExpMat                           = kron(ExpVec,1.d0./ExpVec');  
  
  iReac = 1;
  for iComp = 1:NComp
    QTran(iComp)             = (2.d0 .* pi .* ComponentMass(iComp) .* DSWtoKg ./ AvN .* UKb .* T0_Vec(iT) ./ Plnck.^2).^(3/2 .* RxLxIdx(iComp));
  end
  IntDeg(:)                  = ComponentDeg(:).^(RxLxIdx(:));
  EqVec(LevelsToWrite(:))    =      prod(QTran) .* prod(IntDeg)  .* Levelg(LevelsToWrite(:),1)   .* ExpVec(LevelsToWrite(:),1);
  LogEqVec(LevelsToWrite(:)) = -log(prod(QTran) .* prod(IntDeg)) .* LevelEeV(LevelsToWrite(:),1) .* Ue ./ UKb;
  
  KDiss(LevelsToWrite(:),1)     = ProcessesRates(LevelsToWrite(:),1,1)' .* 1.d-6;
  KEq(LevelsToWrite(:),1)       = ( EqVec(LevelsToWrite(:)) ).^(-1);
  KRecOverg(LevelsToWrite(:),1) = KDiss(LevelsToWrite(:),1) ./ KEq(LevelsToWrite(:),1);
  KRec(LevelsToWrite(:),1)      = KRecOverg(LevelsToWrite(:),1) .* Levelg(LevelsToWrite(:),1);
  
  if (WriteExchangeFlg == 1)
    Kij(LevelsToWrite(:),LevelsToWrite(:)) = (RatesMatrix(LevelsToWrite(:),LevelsToWrite(:),1,1) + RatesMatrix(LevelsToWrite(:),LevelsToWrite(:),2,1) + RatesMatrix(LevelsToWrite(:),LevelsToWrite(:),3,1)) .* 1.d-6;
  else
    Kij(LevelsToWrite(:),LevelsToWrite(:)) = RatesMatrix(LevelsToWrite(:),LevelsToWrite(:),1,1) .* 1.d-6;
  end
  Kij(Kij < MinValA)              = 0.d0;
  AInel                           = tril(Kij .* ExpMat,-1) + tril(Kij,-1)'; 
  wOut                            = sum(AInel,1);
  wIn                             = sum(AInel,2);
  AInel = (AInel - diag(wOut));
  clear Kij Kji ExpMat
  
  NLevelsTemp = size(LevelsToWrite,1);


  
  %% ROBYN's N3 CASE
  P0    = 9940.737d0;
  X0    = 0.05d0;
  T0    = 300;
  %% DANIL's O3 CASE
  %P0    = 1.38064852e+03;
  %X0    = 0.999d0;
  %T0    = 100;
  %%
  V0    = 1.d0;
  R     = 8.3144598;
  n0    = (P0 * V0 / (R * T0)) * AvN;
  nO0   = n0 * X0;

  AInel = AInel .* nO0;      
  clear Kij

  
  tic
  
  
  if WriteDirUndFlg(1) == 1
  %% WRITING DIRECTED NETWORK

      MinValA = 1.d-200;
      NEdges  = sum(sum(abs(AInel) > MinValA));
      if (WritePajekFlg == 0)

          if (WriteOrigInvFlg(1)==1)
              FileName = strcat(FolderName,'/InelEdges_Directed.csv');
              fileID1   = fopen(FileName,'w');
              fprintf(fileID1,'Source,Target,Weight\n');
          end
          if (WriteOrigInvFlg(2)==1)
              FileName = strcat(FolderName,'/InelEdges_Directed_Inverse.csv');
              fileID2   = fopen(FileName,'w');
              fprintf(fileID2,'Source,Target,Weight\n');
          end
          ii=1;
          for i = LevelsToWrite'
            jj=1;
            for j = LevelsToWrite'
              %if i == j
                if abs(AInel(j,i)) >= MinValA
                  if (WriteOrigInvFlg(1)==1)
                    fprintf(fileID1, '%i,%i,%e\n', ii, jj, abs(AInel(j,i)) );
                  end
                  if (WriteOrigInvFlg(2)==1)
                    fprintf(fileID2,'%i,%i,%e\n', ii, jj, abs(1.0./AInel(j,i)) );
                  end
                end
              %end
              jj=jj+1;
            end
            ii=ii+1;
          end
          if (WriteOrigInvFlg(1)==1)
            fclose(fileID1);
          end
          if (WriteOrigInvFlg(2)==1)
            fclose(fileID2);
          end

      else

          if (WriteOrigInvFlg(1)==1)
              FileName = strcat(FolderName,'/InelNetwork_Directed.net');
              fileID1   = fopen(FileName,'w');
              fprintf(fileID1,'# A network in Pajeks .net format\n');
              fprintf(fileID1,'*Vertices %i\n', NLevelsTemp);
              ii=1;
              for i = LevelsToWrite'
                StrNb = ['"v=', num2str(Levelvqn(i)), ',J=', num2str(Leveljqn(i)), '"'];
                fprintf(fileID1,'%i %s\n', ii, StrNb);
                ii=ii+1;
              end
              fprintf(fileID1,'*Arcs %i\n',NEdges);
          end
          if (WriteOrigInvFlg(2)==1)
              FileName = strcat(FolderName,'/InelNetwork_Directed_Inverse.net');
              fileID2   = fopen(FileName,'w');
              fprintf(fileID2,'# A network in Pajeks .net format\n');
              fprintf(fileID2,'*Vertices %i\n', NLevelsTemp);
              ii=1;
              for i = LevelsToWrite'
               StrNb = ['"v=', num2str(Levelvqn(i)), ',J=', num2str(Leveljqn(i)), '"'];
               fprintf(fileID2,'%i %s\n', ii, StrNb);
               ii=ii+1;
              end
              fprintf(fileID2,'*Arcs %i\n',NEdges);
          end
          ii=1;
          for i = LevelsToWrite'
            jj=1;
            for j = LevelsToWrite'
              %if i == j
                if abs(AInel(j,i)) >= MinValA
                  if (WriteOrigInvFlg(1)==1)
                    fprintf(fileID1, '%i %i %e\n', ii, jj, abs(AInel(j,i)) );
                  end
                  if (WriteOrigInvFlg(2)==1)
                    fprintf(fileID2,'%i %i %e\n', ii, jj, abs(1.0./AInel(j,i)) );
                  end
                end
              %end
              jj=jj+1;
            end
            ii=ii+1;
          end
          if (WriteOrigInvFlg(1)==1)
            fclose(fileID1);
          end
          if (WriteOrigInvFlg(2)==1)DeltaEMax
            fclose(fileID2);
          end

      end
      
  end
  
  
  if WriteDirUndFlg(2) == 1
  %% WRITING UNDIRECTED NETWORK

      AInelSimm = (AInel + AInel') ./ 2.d0;
      AInelSimm = AInelSimm - diag(diag(AInelSimm));
      MinValA   = 1.d-200;
      NEdges    = sum(sum(abs(AInelSimm) > MinValA));
      NEdges    = NEdges./2;
      if (WritePajekFlg==0)

          if (WriteOrigInvFlg(1)==1)
              FileName = strcat(FolderName,'/InelEdges_Undirected.csv');
              fileID1   = fopen(FileName,'w');
              fprintf(fileID1,'Source,Target,Weight\n');
          end
          if (WriteOrigInvFlg(2)==1)
              FileName = strcat(FolderName,'/InelEdges_Undirected_Inverse.csv');
              fileID2   = fopen(FileName,'w');
              fprintf(fileID2,'Source,Target,Weight\n');
          end
          ii=1;
          for i = LevelsToWrite'
            jj=1;
            for j = LevelsToWrite'
              if i < j
                if abs(AInelSimm(j,i)) >= MinValA
                  if (WriteOrigInvFlg(1)==1)
                    fprintf(fileID1, '%i,%i,%e\n', ii, jj, abs(AInelSimm(j,i)) );
                  end
                  if (WriteOrigInvFlg(2)==1)
                    fprintf(fileID2,'%i,%i,%e\n', ii, jj, abs(1.0./AInelSimm(j,i)) );
                  end
                end
              end
              jj=jj+1;
            end
            ii=ii+1;
          end
          if (WriteOrigInvFlg(1)==1)
            fclose(fileID1);
          end
          if (WriteOrigInvFlg(2)==1)
            fclose(fileID2);
          end

      else

          if (WriteOrigInvFlg(1)==1)
              FileName = strcat(FolderName,'/InelNetwork_Undirected.net');
              fileID1   = fopen(FileName,'w');
              fprintf(fileID1,'# A network in Pajeks .net format\n');
              fprintf(fileID1,'*Vertices %i\n', NLevelsTemp);
              ii = 1;
              for i = LevelsToWrite'
               StrNb = ['"v=', num2str(Levelvqn(i)), ',J=', num2str(Leveljqn(i)), '"'];
               fprintf(fileID1,'%i %s\n', ii, StrNb);
               ii=ii+1;
              end
              fprintf(fileID1,'*Arcs %i\n',NEdges);
          end
          if (WriteOrigInvFlg(2)==1)
              FileName = strcat(FolderName,'/InelNetwork_Undirected_Inverse.net');
              fileID2   = fopen(FileName,'w');
              fprintf(fileID2,'# A network in Pajeks .net format\n');
              fprintf(fileID2,'*Vertices %i\n', NLevelsTemp);
              ii = 1;
              for i = LevelsToWrite'
                StrNb = ['"v=', num2str(Levelvqn(i)), ',J=', num2str(Leveljqn(i)), '"'];
                fprintf(fileID2,'%i %s\n', ii, StrNb);
                ii=ii+1;
              end
              fprintf(fileID2,'*Arcs %i\n',NEdges);
          end
          ii=1;
          for i = LevelsToWrite'
            jj=1;
            for j = LevelsToWrite'
              if i < j
                if abs(AInelSimm(j,i)) >= MinValA
                  if (WriteOrigInvFlg(1)==1)
                    fprintf(fileID1, '%i %i %e\n', ii, jj, abs(AInelSimm(j,i)) );
                  end
                  if (WriteOrigInvFlg(2)==1)
                    fprintf(fileID2,'%i %i %e\n', ii, jj, abs(1.0./AInelSimm(j,i)) );
                  end
                end
              end
              jj=jj+1;
            end
            ii=ii+1;
          end
          if (WriteOrigInvFlg(1)==1)
            fclose(fileID1);
          end
          if (WriteOrigInvFlg(2)==1)
            fclose(fileID2);
          end

      end

  end
  
  clear AInel AInelSimm
  timee = toc;
  fprintf('Inelastic Rates Matrixes written in %e s\n', timee)
  

  %% WRITING LEVELS
  FileName = strcat(FolderName,'/InelLevels.csv');
  fileID = fopen(FileName,'w');
  fprintf(fileID,'id,Level,v,Longitude,Latitude,rIn,EeVVib,EeVRot,DeltaEBarrier,Bin\n');
  ii=1;
  for i = LevelsToWrite'
    fprintf(fileID,'%i,%i,%i,%e,%e,%e,%e,%e,%e,%i\n', ii, i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), LevelEeVVib0(i), LevelEeVRot(i), DeltaEintDiss(i), BinDiss(i));
    ii=ii+1;
  end
  fclose(fileID);
  
  if DeltaEMax < 0.d0
      
      FileName = strcat(FolderName,'/DissLevels.csv');
      fileID = fopen(FileName,'w');
      fprintf(fileID,'id,Level,v,Longitude,Latitude,rIn,EeVVib,EeVRot,DeltaEBarrier,Bin\n');
      ii=1;
      for i = LevelsToNOTWrite'
        fprintf(fileID,'%i,%i,%i,%e,%e,%e,%e,%e,%e,%i\n', ii, i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), LevelEeVVib0(i), LevelEeVRot(i), DeltaEintDiss(i), BinDiss(i));
        ii=ii+1;
      end
      fclose(fileID);
  
  end
  
end