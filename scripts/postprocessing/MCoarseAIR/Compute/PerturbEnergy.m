function [iFigure] = PerturbEnergy(iFigure, Levelvqn, Leveljqn, LevelEeV)

  global System AtomMass iLevelToPerturb
  
  DeltaE   = 1.d-2;
  NAngMom  = 5;
  NEint    = 5;
  
  
  mu = AtomMass(1) * AtomMass(2) / ( AtomMass(1) + AtomMass(2) );
  Two    = 2.d0;
  One    = 1.d0
  Half   = 5.d-1;
  Zero   = 0.d0;
  
  filename = '/Users/sventuri/WORKSPACE/CoarseAIR/run_N3/Test/N3/N2/levels_cut.inp';
  startRow = 16;
  formatSpec = '%*116s%15f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  ri = dataArray{:, 1};
  clearvars filename startRow formatSpec fileID dataArray ans;
  
  filename = '/Users/sventuri/Downloads/elevels_g (1).dat';
  formatSpec = '%*39s%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN, 'ReturnOnError', false);
  fclose(fileID);
  LevToBin = dataArray{:, 1};
  clearvars filename formatSpec fileID dataArray ans;
  LevToBin = LevToBin + 1;
  for i = 1:length(LevToBin)
    qnToBin(Levelvqn(i,1)+1,Leveljqn(i,1)+1) = LevToBin(i);
  end
  for i = length(LevToBin)+1:size(LevelEeV,1)
    qnToBin(Levelvqn(i,1)+1,Leveljqn(i,1)+1) = 0;
  end

  
  for iLevel = iLevelToPerturb
    
    jqnIn    = Leveljqn(iLevel,1);
    vqnIn    = Levelvqn(iLevel,1);
    EintIn   = LevelEeV(iLevel,1);
    rBond    = ri(iLevel,1);
    AngMomIn = jqnIn + 0.5d0;
    MomIneIn = mu * rBond^2;
    ErotIn   = 0.5d0 * AngMomIn^2 / MomIneIn;
    
    figure
    cmap = colormap(lines(max(LevToBin)+1));
    
    for iEint = [-NEint:NEint]
      jEint = iEint + NEint + 1;

      for iAngMom = [-NAngMom:NAngMom]
        jAngMom = iAngMom + NAngMom +1;
        
        Erot                   = ErotIn + DeltaE * iAngMom;
        AngMom                 = sqrt(2.d0 * Erot * MomIneIn);
        jqn0                   = AngMom - Half;
        Eint                   = 0.d0 + DeltaE * iEint;
        EintMat(jAngMom,jEint) = Eint;
        jqnMat(jAngMom,jEint)  = jqn0;
        
        if Erot > 0.d0
          
          [riFin, roFin, jqnFin, vqnFin, iTypeFin] = AnalyzeFinalState(jqnIn, AngMom, Eint, rBond);
          iTypeFin

          if iTypeFin == 0
            
            vFin(jEint,jAngMom) = vqnFin;
            jFin(jEint,jAngMom) = jqnFin;

            BinMat(jAngMom,jEint)     = qnToBin(vqnFin+1,jqnFin+1);
            
          else
          
            BinMat(jAngMom,jEint)     = 0;

          end
          
        else
          
          BinMat(jAngMom,jEint)     = 0;
        
        end
        
        Colorr(jAngMom,jEint,1:3) = cmap(BinMat(jAngMom,jEint)+1,1:3);

      end

    end

    figure(iFigure)
    surf(EintMat,jqnMat,BinMat,Colorr)
    cmap = colormap(lines(max(LevToBin)+1));
    colorbar
    ylabel('AngMom');
    zlabel('FinalBin');
    iFigure = iFigure + 1;

    figure(iFigure)
    surf(EintMat,jqnMat,vFin)
    xlabel('Eint');
    ylabel('AngMom');
    iFigure = iFigure + 1;

    figure(iFigure)
    surf(EintMat,jqnMat,jFin)
    xlabel('Eint');
    ylabel('AngMom');
    iFigure = iFigure + 1;
    
  end

end