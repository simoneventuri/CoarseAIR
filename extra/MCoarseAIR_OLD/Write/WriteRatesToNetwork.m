function [] = WriteRatesToNetwork(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, Levelg, RatesMatrix, ProcessesRates, LevelEeVVib0, LevelEeVRot, DeltaEintDiss)

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
  
  global ConsiderDegFlg WriteExchangeFlg WriteOrigInvFlg WritePajekFlg
  
  MinValA            = 1.d-14 .* 1.d-6;

  iBinnedMol         = 1;  
  FolderName         = strcat(MoleculesName(iBinnedMol,:),'_Network')
  [status,msg,msgID] = mkdir(FolderName)
  
  
  T0 = 300;
  ExpVec0(1:NLevels(iBinnedMol),1) = exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0 .* UKb) );
  ExpVec0                          = ExpVec0 ./ sum(ExpVec0);
    
%   if (ConsiderDegFlg == 1)
%     ExpVec(1:NLevels(iBinnedMol),1)  = Levelg(1:NLevels(iBinnedMol),1) .* exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0_Vec(1) .* UKb) );
%   else
%     ExpVec(1:NLevels(iBinnedMol),1)  = exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0_Vec(1) .* UKb) );
%   end
%   ExpVec                           = ExpVec ./ sum(ExpVec);
%   ExpMat                           = kron(ExpVec,1.d0./ExpVec');  %%% corect one
% 
%   iReac = 1;
%   for iComp = 1:NComp
%     QTran(iComp) = (2.d0 .* pi .* ComponentMass(iComp) .* DSWtoKg ./ AvN .* UKb .* T0_Vec(iT) ./ Plnck.^2).^(3/2 .* RxLxIdx(iComp));
%   end
%   IntDeg(:)  = ComponentDeg(:).^(RxLxIdx(:));
%   EqVec      =      prod(QTran) .* prod(IntDeg)  .* Levelg(1:NLevels(iBinnedMol),1) .* ExpVec(1:NLevels(iBinnedMol),1);
%   LogEqVec   = -log(prod(QTran) .* prod(IntDeg)) .* LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ UKb;
%   
%   KDiss(1:NLevels(iBinnedMol),1)     = ProcessesRates(1:NLevels(iBinnedMol),1,1)' .* 1.d-6;
%   KEq(1:NLevels(iBinnedMol),1)       = ( EqVec(1:NLevels(iBinnedMol)) ).^(-1);
%   KRecOverg(1:NLevels(iBinnedMol),1) = KDiss(1:NLevels(iBinnedMol),1) ./ KEq(1:NLevels(iBinnedMol),1);
%   KRec(1:NLevels(iBinnedMol),1)      = KRecOverg(1:NLevels(iBinnedMol),1) .* Levelg(1:NLevels(iBinnedMol),1);
%   
  if (WriteExchangeFlg == 1)
    Kij(:,:)                      = (RatesMatrix(1:NBins(iBinnedMol),1:NBins(iBinnedMol),1,1) + RatesMatrix(1:NBins(iBinnedMol),1:NBins(iBinnedMol),2,1) + RatesMatrix(1:NBins(iBinnedMol),1:NBins(iBinnedMol),3,1)) .* 1.d-6;
  else
    Kij(:,:)                      = RatesMatrix(1:NBins(iBinnedMol),1:NBins(iBinnedMol),1,1) .* 1.d-6;
  end
  Kij(Kij < MinValA)              = 0.d0;
  AInel                           = Kij;%tril(Kij .* ExpMat,-1) + tril(Kij,-1)'; %%% Kji = AInel
  wOut                            = sum(AInel,1);
  wIn                             = sum(AInel,2);
  AInel = (AInel - diag(wOut));
  clear Kij Kji ExpMat
  
  NLevelsTemp = NBins(iBinnedMol);%1000
  
  AInell = AInel(1:NLevelsTemp,1:NLevelsTemp);
  clear AInel
  AInel  = AInell;


  
  %% ROBYN's N3 CASE
  P0    = 9940.737d0;
  X0    = 0.05d0;
  T0    = 1400;
  %% DANIL's O3 CASE
  %P0    = 1.38064852e+03;  %%% these values are diff
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
      for i = 1:NLevelsTemp
        for j = 1:NLevelsTemp
          %if i == j
            if abs(AInel(j,i)) >= MinValA
              if (WriteOrigInvFlg(1)==1)
                fprintf(fileID1, '%i,%i,%e\n', i, j, abs(AInel(j,i)) );
              end
              if (WriteOrigInvFlg(2)==1)
                fprintf(fileID2,'%i,%i,%e\n', i, j, abs(1.0./AInel(j,i)) );
              end
            end
          %end
        end
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
          for i = 1:NLevelsTemp
           StrNb = ['"v=', num2str(Levelvqn(i)), ',J=', num2str(Leveljqn(i)), '"'];
           fprintf(fileID1,'%i %s\n', i, StrNb);
          end
          fprintf(fileID1,'*Arcs %i\n',NEdges);
      end
      if (WriteOrigInvFlg(2)==1)
          FileName = strcat(FolderName,'/InelNetwork_Directed_Inverse.net');
          fileID2   = fopen(FileName,'w');
          fprintf(fileID2,'# A network in Pajeks .net format\n');
          fprintf(fileID2,'*Vertices %i\n', NLevelsTemp);
          for i = 1:NLevelsTemp
           StrNb = ['"v=', num2str(Levelvqn(i)), ',J=', num2str(Leveljqn(i)), '"'];
           fprintf(fileID2,'%i %s\n', i, StrNb);
          end
          fprintf(fileID2,'*Arcs %i\n',NEdges);
      end
      for i = 1:NLevelsTemp
        for j = 1:NLevelsTemp
          %if i == j
            if abs(AInel(j,i)) >= MinValA
              if (WriteOrigInvFlg(1)==1)
                fprintf(fileID1, '%i %i %e\n', i, j, abs(AInel(j,i)) );
              end
              if (WriteOrigInvFlg(2)==1)
                fprintf(fileID2,'%i %i %e\n', i, j, abs(1.0./AInel(j,i)) );
              end
            end
          %end
        end
      end
      if (WriteOrigInvFlg(1)==1)
        fclose(fileID1);
      end
      if (WriteOrigInvFlg(2)==1)
        fclose(fileID2);
      end
      
  end
  
  
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
      for i = 1:NLevelsTemp
        for j = i+1:NLevelsTemp
          %if i == j
            if abs(AInelSimm(j,i)) >= MinValA
              if (WriteOrigInvFlg(1)==1)
                fprintf(fileID1, '%i,%i,%e\n', i, j, abs(AInelSimm(j,i)) );
              end
              if (WriteOrigInvFlg(2)==1)
                fprintf(fileID2,'%i,%i,%e\n', i, j, abs(1.0./AInelSimm(j,i)) );
              end
            end
          %end
        end
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
          for i = 1:NLevelsTemp
           StrNb = ['"v=', num2str(Levelvqn(i)), ',J=', num2str(Leveljqn(i)), '"'];
           fprintf(fileID1,'%i %s\n', i, StrNb);
          end
          fprintf(fileID1,'*Arcs %i\n',NEdges);
      end
      if (WriteOrigInvFlg(2)==1)
          FileName = strcat(FolderName,'/InelNetwork_Undirected_Inverse.net');
          fileID2   = fopen(FileName,'w');
          fprintf(fileID2,'# A network in Pajeks .net format\n');
          fprintf(fileID2,'*Vertices %i\n', NLevelsTemp);
          for i = 1:NLevelsTemp
           StrNb = ['"v=', num2str(Levelvqn(i)), ',J=', num2str(Leveljqn(i)), '"'];
           fprintf(fileID2,'%i %s\n', i, StrNb);
          end
          fprintf(fileID2,'*Arcs %i\n',NEdges);
      end
      for i = 1:NLevelsTemp
        for j = i+1:NLevelsTemp
          %if i == j
            if abs(AInelSimm(j,i)) >= MinValA
              if (WriteOrigInvFlg(1)==1)
                fprintf(fileID1, '%i %i %e\n', i, j, abs(AInelSimm(j,i)) );
              end
              if (WriteOrigInvFlg(2)==1)
                fprintf(fileID2,'%i %i %e\n', i, j, abs(1.0./AInelSimm(j,i)) );
              end
            end
          %end
        end
      end
      if (WriteOrigInvFlg(1)==1)
        fclose(fileID1);
      end
      if (WriteOrigInvFlg(2)==1)
        fclose(fileID2);
      end

  end
  
  clear AInel AInelSimm
  timee = toc;
  fprintf('Inelastic Rates Matrixes written in %e s\n', timee)
  

  %% WRITING LEVELS
  FileName = strcat(FolderName,'/InelLevels.csv');
  fileID = fopen(FileName,'w');
  fprintf(fileID,'id,v,Longitude,Latitude,rIn,EeVVib,EeVRot\n');
  for i = 1:NLevelsTemp
    fprintf(fileID,'%i,%i,%e,%e,%e,%e,%e\n', i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), LevelEeVVib0(i), LevelEeVRot(i));
  end
  fclose(fileID);
  
  
  %% WRITING LEVELS
  FileName = strcat(FolderName,'/InelBin.csv');
  fileID = fopen(FileName,'w');
  fprintf(fileID,'id,Bin\n');
  for i = 1:NLevelsTemp
    fprintf(fileID,'%i,%i\n', i, i);
  end
  fclose(fileID);
  
  
end