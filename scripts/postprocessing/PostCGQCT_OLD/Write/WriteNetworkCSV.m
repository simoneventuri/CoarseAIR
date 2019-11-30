function [] = WriteNetworkCSV(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, Levelg, RatesMatrix, ProcessesRates, t, MolFracs, nd, Pop, StpInstants, LevelEeVVib0, LevelEeVRot)

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

  
  MinValA            = 1.d-15 .* 1.d-6;
  MinValW            = 1.d-4;
  MinValW_Bis        = 1.d4;
  
  [status,msg,msgID] = mkdir('./matrix4binning')
  
  iBinnedMol         = 1;
  
  T0 = 300;
  ExpVec0(1:NLevels(iBinnedMol),1) = exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0 .* UKb) );
  ExpVec0                          = ExpVec0 ./ sum(ExpVec0);
    
  ExpVec(1:NLevels(iBinnedMol),1)  = exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0_Vec(1) .* UKb) );
  ExpVec                           = ExpVec ./ sum(ExpVec);
  ExpMat                           = kron(ExpVec,1.d0./ExpVec');  
  
%   LogVec(1:NLevels(iBinnedMol),1) = LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ UKb ;
%   LogMat = reshape(LogVec, [1, size(RatesMatrix,1)]) - reshape(LogVec', [size(RatesMatrix,1), 1]);
  
  iReac = 1;
  for iComp = 1:NComp
    QTran(iComp) = (2.d0 .* pi .* ComponentMass(iComp) .* DSWtoKg ./ AvN .* UKb .* T0_Vec(iT) ./ Plnck.^2).^(3/2 .* RxLxIdx(iComp));
  end
  IntDeg(:)  = ComponentDeg(:).^(RxLxIdx(:));
  EqVec                           = prod(QTran) .* prod(IntDeg) .* Levelg(1:NLevels(iBinnedMol),1) .* ExpVec(1:NLevels(iBinnedMol),1);
  LogEqVec                        = -log(prod(QTran) .* prod(IntDeg)) .* LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ UKb;

  %Kij(:,:)                        = (RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,1) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),2,1)) .* 1.d-6;
  Kij(:,:)                        = RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1,1) .* 1.d-6;
  Kij(Kij < MinValA)              = 0.d0;
  Kji                             = tril(Kij .* ExpMat,-1) + tril(Kij,-1)';%+ diag(diag(Kij));%Kij .* ExpMat;%
  wOut                            = sum(Kji,1);
  wIn                             = sum(Kji,2);
  
  if WriteDissFlag == 1 || WriteAllFlag == 1
    KDiss(1:NLevels(iBinnedMol),1)  = ProcessesRates(1:NLevels(iBinnedMol),1,1)' .* 1.d-6;
    KEq(1:NLevels(iBinnedMol),1)    = ( EqVec(1:NLevels(iBinnedMol)) ).^(-1);
    KRecOverg(:,1)                  = KDiss(:,1) ./ KEq(:,1);
    KRec(:,1)                       = KRecOverg(:,1) .* Levelg(:,1);
  end
   
  if WriteSrcTermsFlag == 1
    i = 0;
    for iSteps = StpInstants(1:end-1)
      i = i+1;
      PopMol = sum( Pop(iSteps,:,iBinnedMol) );
      for iLevels = 1:NLevels(iBinnedMol)
        n(i,iLevels) = Pop(iSteps,iLevels,iBinnedMol) ./ Levelg(iLevels,iBinnedMol);%./ PopMol .* MolFracs(iSteps,CompOI) .* nd(iSteps);
      end
      n(i,NLevels(iBinnedMol)+1) = MolFracs(iSteps,1) .* nd(iSteps);
      n(i,NLevels(iBinnedMol)+2) = MolFracs(iSteps,2) .* nd(iSteps);
    end
  end
  clear Kij ExpMat
  clear Pop RatesMatrix ProcessesRates
  
  
  if WriteRatesFlag == 1
    
    if WriteInelasticFlag == 1
      
      P0    = 9940.737d0;
      X0    = [1.d-20, 0.05d0];
      X0(3) = 1-sum(X0);
      T0    = 300;
      V0    = 1.d0;

      iBinnedMol = 1;
      iSteps     = 1;
      R          = 8.3144598;
      iFigure    = 2001;
      
      n0      = (P0 * V0 / (R * T0)) * AvN;
      nO0     = n0 * X0(2);
      nCOTot0 = n0 * X0(3);

      Temp = (Kji - diag(wOut)) .* nO0;      
%       AInel       = (Kji - diag(diag(Kij))) .* nO0;
%       [AInelExp]  = expm(AInel .* 1.d-12);
%       filenameSVD = strcat('./AInelExp')
%       save(filenameSVD, 'AInel', 'AInelExp', '-v7.3')
%       [U,S,V] = svd(AInel);
%       filenameSVD = strcat('./SVD')
%       save(filenameSVD, 'U', 'S', 'V', '-v7.3')
     
%       MinN    =   90.0d6;
%       MaxN    =   100.0d6;
%       a       =    20;
%       b       =   -20;
%       ExitFlg =  true;
%       while ExitFlg 
%         if sum(sum(abs(Temp) > 10.d0^(mean([a,b])))) < MinN
%           a = mean([a,b])
%         elseif sum(sum(abs(Temp) > 10.d0^(mean([a,b])))) > MaxN
%           b = mean([a,b])
%         else
%           MinValA = 10.d0^(mean([a,b]));
%           ExitFlg = false;
%         end
%       end
%       fprintf('Min Value for Inelastic Rates is %e\n', MinValA)
      MinValA = 1.d-200;
      NEdges  = sum(sum(abs(Temp) > MinValA));
      
             
%       Temp(Temp < MinValA) = 0.d0;
%       NodeNamesStr = num2str([1:1:NLevels(iBinnedMol)]');
%       NodeNames    = cellstr(NodeNamesStr);
%       Graph        = digraph(Temp,NodeNames);
%       figure; plot(Graph,'XData',Leveljqn(1:NLevels(iBinnedMol),iBinnedMol),'YData',LevelEeV(1:NLevels(iBinnedMol),iBinnedMol))


      tic
      FileName = strcat('./InelRates.net');
      %FileName = strcat('/Users/sventuri/Dropbox/TempRes/TransMatrix.net');
      fileID   = fopen(FileName,'w');
      %fprintf(fileID,'Source,Target,Weight,Type\n');
      fprintf(fileID,'# A network in Pajeks .net format\n');
      
      fprintf(fileID,'*Vertices %i\n', 300)%NLevels(iBinnedMol))
      for i = 1:300%NLevels(iBinnedMol)
       StrNb = ['"v=', num2str(Levelvqn(i)), ',J=', num2str(Leveljqn(i)), '"'];
       fprintf(fileID,'%i %s\n', i, StrNb);
      end

      fprintf(fileID,'*Arcs %i\n',NEdges);
      for i = 1:300%NLevels(iBinnedMol)
        for j = 1:300%NLevels(iBinnedMol)
          %if i == j
            if abs(Temp(j,i)) >= MinValA
              %fprintf(fileID,'%i,%i,%e,"direct"\n', j, i, Temp(i,j) );
              fprintf(fileID,'%i %i %e\n', i, j, abs(Temp(j,i)) );
            end
          %end
        end
      end
      fclose(fileID);
      clear AInel
      timee = toc;
      fprintf('Inelastic Rates Matrixes written in %e s\n', timee)
      
      
      FileName = strcat('./InelLevels.csv');
      fileID = fopen(FileName,'w');
      fprintf(fileID,'id,v,Longitude,Latitude,rIn,EeVVib,EeVRot\n');
      for i = 1:300%NLevels(iBinnedMol)
        fprintf(fileID,'%i,%i,%e,%e,%e,%e,%e\n', i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), LevelEeVVib0(i), LevelEeVRot(i));
      end
      fclose(fileID);

    end

  end

  

  if WriteSrcTermsFlag == 1

    if WriteDissFlag == 1
      
      for iSteps = 8%1:length(StpInstants)-1

        WDiss                                              = kron(n(iSteps,1:NLevels(iBinnedMol)),1.d0./n(iSteps,1:NLevels(iBinnedMol))')' ./ ExpMat;  
        WDiss(1:NLevels(iBinnedMol),NLevels(iBinnedMol)+1) = n(iSteps,1:NLevels(iBinnedMol))' ./ (n(iSteps,NLevels(iBinnedMol)+1) .* n(iSteps,NLevels(iBinnedMol)+2)) ./ EqVec(1:NLevels(iBinnedMol));
        WDiss(NLevels(iBinnedMol)+1,NLevels(iBinnedMol)+1) = 0;
        
        WDissNormal = triu(lognpdf(WDiss,1.d-10,1.d-5) ./ lognpdf(1,1.d-10,1.d-5), 1);
        clear WDiss
        
        MinN    = 1.8d6;
        MaxN    = 2.0d6;
        a       =   1;
        b       = -15;
        ExitFlg =  true;
        while ExitFlg 
          if sum(sum(WDissNormal > 10.d0^(mean([a,b])))) < MinN
            a = mean([a,b])
          elseif sum(sum(WDissNormal > 10.d0^(mean([a,b])))) > MaxN
            b = mean([a,b])
          else
            MinValW = 10.d0^mean([a,b]);
            ExitFlg = false;
          end
        end
        fprintf('Min Value for Sources at time %e s is %e\n', t(StpInstants(iSteps)), MinValW)

        
        tic
        FileName = strcat('/Users/sventuri/Dropbox/TempRes/DissRates_',num2str(StpInstants(iSteps)),'.csv');
        fileID = fopen(FileName,'w');
        fprintf(fileID,'Source,Target,Weight,Type\n');

        %       FileNameInverse = strcat('./matrix4binning/DissRatesInverse.csv');
        %       fileIDInverse = fopen(FileNameInverse,'w');
        %       fprintf(fileIDInverse,'Source,Target,Weight,Type\n');
        
        LevelWrtFlg = zeros(size(RatesMatrix,1)+2,1);
        for i = 1:size(RatesMatrix,1)
          for j = 1:size(RatesMatrix,1)
            if i ~= j
              if WDissNormal(j,i) >= MinValW
                fprintf(fileID,'%i,%i,%e,"direct"\n', i, j, (WDissNormal(j,i)-MinValW) / (1.d0-MinValW) + 2);
                LevelWrtFlg(i) = LevelWrtFlg(i)+1;
                LevelWrtFlg(j) = LevelWrtFlg(j)+1;
        %               fprintf(fileIDInverse,'%i,%i,%e,"direct"\n', i, j, -log10(AAll(i,j)) );
              end
            end
          end
        end
%         for i = 1:size(RatesMatrix,1)+2
%           for j = size(RatesMatrix,1)+1:size(RatesMatrix,1)+2
%             if i ~= j
%               if WDissNormal(j,i) >= MinValW
%                 LevelWrtFlg(i) = LevelWrtFlg(i)+1;
%                 LevelWrtFlg(j) = LevelWrtFlg(j)+1;
%                 fprintf(fileID,'%i,%i,%e,"direct"\n', i, j, log10( WDissNormal(j,i) ./ MinValW) );
%         %               fprintf(fileIDInverse,'%i,%i,%e,"direct"\n', i, j, -log10(AAll(i,j)) );
%               end
%             end
%           end
%         end
        fclose(fileID);
        clear WDissNormal
        %       fclose(fileIDInverse);
        timee = toc;
        fprintf('Diss Rates Matrixes at time %e s written in %e s\n', t(StpInstants(iSteps)), timee)

        
%         FileName = strcat('/Users/sventuri/Dropbox/TempRes/LevelsDiss_',num2str(StpInstants(iSteps)),'.csv');
%         fileID = fopen(FileName,'w');
%         fprintf(fileID,'id,v,Longitude,Latitude,rIn,Pop\n');
%         for i = 1:size(RatesMatrix,1)
%           if LevelWrtFlg(i) ~= 0
%             fprintf(fileID,'%i,%i,%i,%e,%e,%e\n', i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), log10(n(iSteps,i)));
%           end
%         end
%         fprintf(fileID,'%i,%i,%i,%e,%e,%e\n', size(RatesMatrix,1)+1, -1, -150, 30.d0, -1.d0, log10(n(iSteps,size(RatesMatrix,1)+1)));
%         fprintf(fileID,'%i,%i,%i,%e,%e,%e\n', size(RatesMatrix,1)+2, -2, -170, 30.d0, -2.d0, log10(n(iSteps,size(RatesMatrix,1)+2)));
%         fclose(fileID);   

        FileName = strcat('/Users/sventuri/Dropbox/TempRes/LevelsDiss_',num2str(StpInstants(iSteps)),'.csv');
        fileID = fopen(FileName,'w');
        fprintf(fileID,'id,v,Longitude,Latitude,rIn,Pop\n');
        for i = 1:size(RatesMatrix,1)
          %if LevelWrtFlg(i) ~= 0
            fprintf(fileID,'%i,%i,%i,%e,%e,%e\n', i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), log10(n(iSteps,i)));
          %end
        end
        fprintf(fileID,'%i,%i,%i,%e,%e,%e\n', size(RatesMatrix,1)+1, -1, -150, 30.d0, -1.d0, log10(n(iSteps,size(RatesMatrix,1)+1) .* n(iSteps,size(RatesMatrix,1)+2)) );
        fclose(fileID);   
      
      end

    end
   
    
    if WriteInelasticFlag == 1
      
      for iSteps = 1:length(StpInstants)-1
        
        WInel       = kron(n(iSteps,1:NLevels(iBinnedMol)),1.d0./n(iSteps,1:NLevels(iBinnedMol))')' ./ ExpMat;  
        WInelNormal = triu(lognpdf(WInel,1/100,1/10) ./ lognpdf(1,1/100,1/10), 1);
        clear WInel
        
%         WInel       = log( kron(n(iSteps,1:NLevels(iBinnedMol)),1.d0./n(iSteps,1:NLevels(iBinnedMol))')' );  
%         WInelNormal = triu(LogMat ./ WInel, 1);
%         clear WInel

%         WInel = (Kji - diag(sum(Kij,2))) .* n(iSteps,NLevels(iBinnedMol)+2);
%         %WInelNormal = WInel .* repmat(n(iSteps,1:size(RatesMatrix,1)),[size(RatesMatrix,1), 1]) ./ repmat(n(iSteps,1:size(RatesMatrix,1))',[1, size(RatesMatrix,1)]);
%         WInelNormal = WInel ./ repmat(n(iSteps,1:size(RatesMatrix,1))',[1, size(RatesMatrix,1)]);
%         clear WInel
        MinN    = 2.8d6;
        MaxN    = 3.0d6;
        a       =   10;
        b       = -15;
        ExitFlg =  true;
        while ExitFlg 
          if sum(sum(WInelNormal > 10.d0^(mean([a,b])))) < MinN
            a = mean([a,b])
          elseif sum(sum(WInelNormal > 10.d0^(mean([a,b])))) > MaxN
            b = mean([a,b])
          else
            MinValW = 10.d0^(mean([a,b]));
            ExitFlg = false;
          end
        end
        fprintf('Min Value for Inelastic Source Terms is %e\n', MinValW)
        %MinValW = 1.d-100;
        
        tic
        FileName = strcat('/Users/sventuri/Dropbox/TempRes/InelRates_',num2str(StpInstants(iSteps)),'.csv');
        fileID = fopen(FileName,'w');
        fprintf(fileID,'Source,Target,Weight,Type\n');

        %       FileNameInverse = strcat('./matrix4binning/DissRatesInverse.csv');
        %       fileIDInverse = fopen(FileNameInverse,'w');
        %       fprintf(fileIDInverse,'Source,Target,Weight,Type\n');
        
        LevelWrtFlg = zeros(size(RatesMatrix,1),1);
        for i = 1:size(RatesMatrix,1)
          for j = i+1:size(RatesMatrix,1)%1:size(RatesMatrix,1)
            if i ~= j
              if WInelNormal(i,j) >= MinValW
                fprintf(fileID,'%i,%i,%e,"direct"\n', j, i, WInelNormal(i,j) );
                LevelWrtFlg(i) = LevelWrtFlg(i)+1;
                LevelWrtFlg(j) = LevelWrtFlg(j)+1;
        %               fprintf(fileIDInverse,'%i,%i,%e,"direct"\n', i, j, -log10(AAll(i,j)) );
              end
            end
          end
        end
        fclose(fileID);
        clear WInelNormal
        %       fclose(fileIDInverse);
        timee = toc;
        fprintf('Inelastic Rates Matrixes at time %e s written in %e s\n', t(StpInstants(iSteps)), timee)

        
        FileName = strcat('/Users/sventuri/Dropbox/TempRes/LevelsInel_',num2str(StpInstants(iSteps)),'.csv');
        fileID = fopen(FileName,'w');
        fprintf(fileID,'id,v,Longitude,Latitude,rIn,Pop\n');
        for i = 1:size(RatesMatrix,1)
          %if LevelWrtFlg(i) ~= 0
            fprintf(fileID,'%i,%i,%i,%e,%e,%e\n', i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), log10(n(iSteps,i)));
          %end
        end
        fclose(fileID);
      
      end

    end
    

    
    if WriteAllFlag == 1

      for iSteps = 11%1:length(StpInstants)-1
        
    
        nO  = n(iSteps,NLevels(iBinnedMol)+2);
        nCO = n(iSteps,1:NLevels(iBinnedMol))';
        
        wInWeighted = sum(Kji .* repmat(nCO',[NLevels(iBinnedMol),1]), 2);

        Temp = (Kji - diag(wOut + KDiss)) .* nO;      
        Temp = [Temp; KDiss' .* nO];
        Temp = [Temp, [ -(wOut' + KDiss).*nCO + wInWeighted; sum(KDiss.*nCO) ]];
  %       AInel       = (Kji - diag(diag(Kij))) .* nO0;
  %       [AInelExp]  = expm(AInel .* 1.d-12);
  %       filenameSVD = strcat('./AInelExp')
  %       save(filenameSVD, 'AInel', 'AInelExp', '-v7.3')
  %       [U,S,V] = svd(AInel);
  %       filenameSVD = strcat('./SVD')
  %       save(filenameSVD, 'U', 'S', 'V', '-v7.3')

        MinN    =   90.0d6;
        MaxN    =   100.0d6;
        a       =    20;
        b       =   -20;
        ExitFlg =  true;
        while ExitFlg 
          if sum(sum(abs(Temp) > 10.d0^(mean([a,b])))) < MinN
            a = mean([a,b])
          elseif sum(sum(abs(Temp) > 10.d0^(mean([a,b])))) > MaxN
            b = mean([a,b])
          else
            MinValA = 10.d0^(mean([a,b]));
            ExitFlg = false;
          end
        end
        fprintf('Min Value for Inelastic Rates is %e\n', MinValA)
        NEdges = sum(sum(abs(Temp) > 10.d0^(mean([a,b]))));


  %       Temp(Temp < MinValA) = 0.d0;
  %       NodeNamesStr = num2str([1:1:NLevels(iBinnedMol)]');
  %       NodeNames    = cellstr(NodeNamesStr);
  %       Graph        = digraph(Temp,NodeNames);
  %       figure; plot(Graph,'XData',Leveljqn(1:NLevels(iBinnedMol),iBinnedMol),'YData',LevelEeV(1:NLevels(iBinnedMol),iBinnedMol))


        tic
        FileName = strcat('/Users/sventuri/Dropbox/TempRes/AllRates_', num2str(t(StpInstants(iSteps))), '.net');
        %FileName = strcat('/Users/sventuri/Dropbox/TempRes/TransMatrix.net');
        fileID   = fopen(FileName,'w');
        %fprintf(fileID,'Source,Target,Weight,Type\n');
        fprintf(fileID,'# A network in Pajeks .net format\n');

        fprintf(fileID,'*Vertices %i\n',NLevels(iBinnedMol)+1)
        for i = 1:NLevels(iBinnedMol)
         StrNb = ['"v=', num2str(Levelvqn(i)), ',J=', num2str(Leveljqn(i)), '"'];
         fprintf(fileID,'%i %s\n', i, StrNb);
        end
        StrNb = ['"O Atom"'];
        fprintf(fileID,'%i %s\n', NLevels(iBinnedMol)+1, StrNb);

        fprintf(fileID,'*Arcs %i\n',NEdges);
        for i = 1:NLevels(iBinnedMol)+1
          for j = 1:NLevels(iBinnedMol)+1
            %if i == j
              if abs(Temp(j,i)) >= MinValA
                %fprintf(fileID,'%i,%i,%e,"direct"\n', j, i, Temp(i,j) );
                fprintf(fileID,'%i %i %e\n', i, j, abs(Temp(j,i)) );
              end
            %end
          end
        end
        fclose(fileID);
        clear AInel
        timee = toc;
        fprintf('Inelastic Rates Matrixes written in %e s\n', timee)


        FileName = strcat('/Users/sventuri/Dropbox/TempRes/AllLevels.csv');
        fileID = fopen(FileName,'w');
        fprintf(fileID,'id,v,Longitude,Latitude,rIn,EeVVib,EeVRot\n');
        for i = 1:NLevels(iBinnedMol)
          fprintf(fileID,'%i,%i,%e,%e,%e,%e,%e\n', i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), LevelEeVVib0(i), LevelEeVRot(i));
        end
        fclose(fileID);
        fprintf(fileID,'%i,%i,%e,%e,%e,%e,%e\n', i, -1, -1.d0, 1.d0, 0.d0, 0.d0, 0.d0);

        
% %         WAll                            = ( Kji - diag(sum(Kij,2)) - diag(KDiss) ) .* n(iSteps,NLevels(iBinnedMol)+2);
% %         cVecCol                         = KRecOverg .* n(iSteps,NLevels(iBinnedMol)+2)^2 ;
% %         oVecCol                         = 2.d0 .* KRecOverg .* n(iSteps,NLevels(iBinnedMol)+1) .* n(iSteps,NLevels(iBinnedMol)+2) - KDiss .* n(iSteps,1:NLevels(iBinnedMol))';
% %         WAll                            = [WAll, cVecCol, oVecCol];
% %         cVecRow                         = [(KDiss(1:NLevels(iBinnedMol),iBinnedMol) .* Levelg(1:NLevels(iBinnedMol),iBinnedMol))' .* n(iSteps,NLevels(iBinnedMol)+2), -sum(KRec) .* n(iSteps,NLevels(iBinnedMol)+2)^2, -2.d0 .* sum(KRec) .* n(iSteps,NLevels(iBinnedMol)+1) .* n(iSteps,NLevels(iBinnedMol)+2) + sum(KDiss(1:NLevels(iBinnedMol)) .* Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* n(iSteps,1:NLevels(iBinnedMol))') ];
% %         oVecRow                         = [(KDiss(1:NLevels(iBinnedMol),iBinnedMol) .* Levelg(1:NLevels(iBinnedMol),iBinnedMol))' .* n(iSteps,NLevels(iBinnedMol)+2), -sum(KRec) .* n(iSteps,NLevels(iBinnedMol)+2)^2, -2.d0 .* sum(KRec) .* n(iSteps,NLevels(iBinnedMol)+1) .* n(iSteps,NLevels(iBinnedMol)+2) + sum(KDiss(1:NLevels(iBinnedMol)) .* Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* n(iSteps,1:NLevels(iBinnedMol))') ];
% %         WAll                            = [WAll; cVecRow; oVecRow];
% %         WAllNormal = WAll ./ repmat(n(iSteps,1:size(RatesMatrix,1)+2)',[1, size(RatesMatrix,1)+2]);
% % %        WAllNormal = WAll .* repmat(n(iSteps,1:size(RatesMatrix,1)+2),[size(RatesMatrix,1)+2, 1]) ./ repmat(n(iSteps,1:size(RatesMatrix,1)+2)',[1, size(RatesMatrix,1)+2]);
% %          clear WAll
%          
% %         MinN    = 2.8d6;
% %         MaxN    = 3.0d6;
% %         a       =   10;
% %         b       = -15;
% %         ExitFlg =  true;
% %         while ExitFlg 
% %           if sum(sum(WAllNormal > 10.d0^(mean([a,b])))) < MinN
% %             a = mean([a,b])
% %           elseif sum(sum(WAllNormal > 10.d0^(mean([a,b])))) > MaxN
% %             b = mean([a,b])
% %           else
% %             MinValW = 10.d0^mean([a,b]);
% %             ExitFlg = false;
% %           end
% %         end
% %         fprintf('Min Value for Sources at time %e s is %e\n', t(StpInstants(iSteps)), MinValW)
% % %        MinValW = min(WAllNormal(WAllNormal>0));
% % 
% %         
% %         tic
% %         FileName = strcat('/Users/sventuri/Dropbox/TempRes/AllRates_',num2str(StpInstants(iSteps)),'.csv');
% %         fileID = fopen(FileName,'w');
% %         fprintf(fileID,'Source,Target,Weight,Type\n');
% % 
% %         %       FileNameInverse = strcat('./matrix4binning/DissRatesInverse.csv');
% %         %       fileIDInverse = fopen(FileNameInverse,'w');
% %         %       fprintf(fileIDInverse,'Source,Target,Weight,Type\n');
% %         
% %         LevelWrtFlg = zeros(size(RatesMatrix,1)+1,1);
% %         for i = 1:size(RatesMatrix,1)+1
% %           
% %           for j = 1:size(RatesMatrix,2)+1
% %             if i ~= j
% %               if WAllNormal(j,i) >= MinValW
% %                 fprintf(fileID,'%i,%i,%e,"direct"\n', i, j, WAllNormal(j,i) );
% %                 LevelWrtFlg(i) = LevelWrtFlg(i)+1;
% %                 LevelWrtFlg(j) = LevelWrtFlg(j)+1;
% %         %               fprintf(fileIDInverse,'%i,%i,%e,"direct"\n', i, j, -log10(AAll(i,j)) );
% %               end
% %             end
% %           end
% %         end
% %         fclose(fileID);
% %         clear WAllNormal
% %         %       fclose(fileIDInverse);
% %         timee = toc;
% %         fprintf('All Rates Matrixes at time %e s written in %e s\n', t(StpInstants(iSteps)), timee)
%         
%         
%         tic
%         FileName = strcat('/Users/sventuri/Dropbox/TempRes/AllRates_',num2str(StpInstants(iSteps)),'.csv');
%         fileID = fopen(FileName,'w');
%         fprintf(fileID,'Source,Target,Weight,Type\n');
% 
%         %       FileNameInverse = strcat('./matrix4binning/DissRatesInverse.csv');
%         %       fileIDInverse = fopen(FileNameInverse,'w');
%         %       fprintf(fileIDInverse,'Source,Target,Weight,Type\n');
%         
%         LevelWrtFlg = zeros(size(RatesMatrix,1)+1,1);
%         for i = 1:size(RatesMatrix,1)+1
%           MaxN     = 250;
%           b        =  0;
%           ExitFlg  =  true;
%           while ExitFlg 
%             if sum(sum(WAllNormal(i,i+1:end) > b)) > MaxN
%               b = b + 20;
%             else
%               MinValW = b;
%               ExitFlg = false;
%             end
%           end
%           for j = i+1:size(RatesMatrix,2)+1
%             if WAllNormal(i,j) >= MinValW
%               fprintf(fileID,'%i,%i,%e,"direct"\n', i, j, WAllNormal(i,j) );
%               LevelWrtFlg(i) = LevelWrtFlg(i)+1;
%               LevelWrtFlg(j) = LevelWrtFlg(j)+1;
%             end
%           end
%         end
%         fclose(fileID);
%         clear WAllNormal
%         %       fclose(fileIDInverse);
%         timee = toc;
%         fprintf('All Rates Matrixes at time %e s written in %e s\n', t(StpInstants(iSteps)), timee)
% 
%         
%         FileName = strcat('/Users/sventuri/Dropbox/TempRes/LevelsAll_',num2str(StpInstants(iSteps)),'.csv');
%         fileID = fopen(FileName,'w');
%         fprintf(fileID,'id,v,Longitude,Latitude,rIn,Pop\n');
%         for i = 1:size(RatesMatrix,1)
%           if LevelWrtFlg(i) ~= 0
%             fprintf(fileID,'%i,%i,%i,%e,%e,%e\n', i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), log10(n(iSteps,i)));
%           end
%         end
%         fprintf(fileID,'%i,%i,%i,%e,%e,%e\n', size(RatesMatrix,1)+1, -1, -150, 30.d0, -1.d0, log10(n(iSteps,size(RatesMatrix,1)+1) .* n(iSteps,size(RatesMatrix,1)+2)));
%         fclose(fileID);   
%         
% %         FileName = strcat('/Users/sventuri/Dropbox/TempRes/LevelsAll_',num2str(StpInstants(iSteps)),'.csv');
% %         fileID = fopen(FileName,'w');
% %         fprintf(fileID,'id,v,Longitude,Latitude,rIn,Pop\n');
% %         for i = 1:size(RatesMatrix,1)
% %           if LevelWrtFlg(i) ~= 0
% %             fprintf(fileID,'%i,%i,%i,%e,%e,%e\n', i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), log10(n(iSteps,i)));
% %           end
% %         end
% %         fprintf(fileID,'%i,%i,%i,%e,%e,%e\n', size(RatesMatrix,1)+1, -1, -150, 30.d0, -1.d0, log10(n(iSteps,size(RatesMatrix,1)+1)));
% %         fprintf(fileID,'%i,%i,%i,%e,%e,%e\n', size(RatesMatrix,1)+2, -2, -170, 30.d0, -2.d0, log10(n(iSteps,size(RatesMatrix,1)+2)));
% %         fclose(fileID);   
      
      end

    end
    
    

  end
  
  
  
end