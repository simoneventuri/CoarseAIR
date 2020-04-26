function [DissRates, ProcessesRates] = CorrectRatesCG(iT, TypeCG, NLevels, Levelvqn, Levelg, LevelEeV0, DeltaEintDiss, DissRates, ProcessesRates)

  global KeV T0_Vec
  
  iBinnedMol = 1;
  
  if TypeCG == 0
    
    ProcessesRates(:,1,iT) = ProcessesRates(:,1,iT) .* 16.0/3.0;
    DissRates(:)           = ProcessesRates(:,1,iT) .* 16.0/3.0;
  elseif TypeCG == 1
    KDiss        = ProcessesRates(:,1,iT);
    NVqn         = max(Levelvqn(:,iBinnedMol))+1;
    QBinsVqn     = zeros(NVqn,iBinnedMol);
    NJqn         = zeros(NVqn,iBinnedMol);
    KDissVqnTemp = zeros(NVqn,iBinnedMol);
    EBin         = zeros(NVqn,iBinnedMol);
    
    for iLevels = 1:NLevels(iBinnedMol)
      vqn           = Levelvqn(iLevels,iBinnedMol) + 1;
      q(iLevels)    = Levelg(iLevels,iBinnedMol) .* exp( - (LevelEeV0(iLevels,iBinnedMol)) ./ (KeV*T0_Vec(iT)) );
      QBinsVqn(vqn) = QBinsVqn(vqn) + q(iLevels);
      NJqn(vqn)     = NJqn(vqn)     + 1;
    end

    for iLevels = 1:NLevels(iBinnedMol)
      vqn               = Levelvqn(iLevels,iBinnedMol) + 1;
      KDissVqnTemp(vqn) = KDissVqnTemp(vqn) + KDiss(iLevels) * q(iLevels) / QBinsVqn(vqn);
    end

    for iLevels = 1:NLevels
      vqn               = Levelvqn(iLevels,iBinnedMol) + 1;
      KDissNew(iLevels) = KDissVqnTemp(vqn);
    end
    
    ProcessesRates(:,1,iT) = KDissNew(:) .* 16.0/3.0;
    DissRates(:)           = KDissNew(:) .* 16.0/3.0;
  elseif TypeCG == 2
    NBins         = 45
    iBinnedMol    = 1
    LSpace        = linspace(0.0,1.0,NBins+1);
    TSpace        = abs(min(DeltaEintDiss(:,iBinnedMol))).*LSpace + min(DeltaEintDiss(:,iBinnedMol));
    ModSpace      = LSpace.^1;%.^(1/3);
    Extr          = abs(min(DeltaEintDiss(:,iBinnedMol))).*ModSpace + min(DeltaEintDiss(:,iBinnedMol));
    KDiss        = ProcessesRates(:,1,1); %max(ProcessesRates(:,1),1.e-20);
    Q            = zeros(NBins,iBinnedMol);
    QEn          = zeros(NBins,iBinnedMol);
    KDissBin     = zeros(NBins,iBinnedMol);
    
    for iLevels = 1:NLevels(iBinnedMol)
      iBin = 1;
      while (DeltaEintDiss(iLevels,iBinnedMol) >= Extr(iBin))
        iBin = iBin + 1;
      end
      iBin = iBin - 1;
      LevelToBin(iLevels) = iBin;
      q(iLevels)          = Levelg(iLevels,iBinnedMol) .* exp( - LevelEeV0(iLevels,iBinnedMol) ./ (KeV*T0_Vec(iT)) );
      Q(iBin)             = Q(iBin)       + q(iLevels);
    end
    q       = q';

    for iLevels = 1:NLevels(iBinnedMol)
      iBin           = LevelToBin(iLevels);
      KDissBin(iBin) = KDissBin(iBin) + KDiss(iLevels) * q(iLevels) / Q(iBin);
    end

    for iLevels = 1:NLevels
      iBin              = LevelToBin(iLevels);
      KDissNew(iLevels) = KDissBin(iBin);
    end
    
    ProcessesRates(:,1,iT) = KDissNew(:) .* 16.0/3.0;
    DissRates(:)           = KDissNew(:) .* 16.0/3.0;
  end
  

end