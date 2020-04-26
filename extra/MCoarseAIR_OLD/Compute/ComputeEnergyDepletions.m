function [CDInt, CDVib, CDRot] = ComputeEnergyDepletions(NLevels, Levelg, LevelEeV0, LevelEeV, LevelEeVVib0, LevelEeVRot, Pop, ProcessesRates, MolFracs, nd, LevelToBin, QBins, Steps)
    
  global NBins AvN MoleculeMu KinMthd Ue T0_Vec UKb

  iBinnedMol = 1;
  
  KDisss   = max(ProcessesRates(1:NBins(iBinnedMol),1,1),1.e-30);
  TempMat = zeros(NLevels, size(Pop,1));

  jStep = 1
  for iStep = 1:size(Pop,1)
    iStep
    NLevels
    
    if (sum(KinMthd(1,:) == 'CGM') == 3) || (sum(KinMthd(1,:) == 'VIB') == 3) 
      Popp(1:NLevels(iBinnedMol),iBinnedMol) = Pop(iStep,LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol)' ./ QBins(LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol) .* Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* exp( - LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol) .* Ue ./ (T0_Vec(1) .* UKb) );
    else
      Popp(1:NLevels(iBinnedMol),iBinnedMol) = Pop(iStep,1:NLevels(iBinnedMol),iBinnedMol)';
    end 
        
    rhoA(iStep)                        = nd(iStep) .* MolFracs(iStep,1);% .* MoleculeMu/2.0 ./ AvN;
    rhoM(iStep)                        = nd(iStep) .* MolFracs(iStep,2);% .* MoleculeMu     ./ AvN;
    PopM(iStep)                        = sum( Popp );
    PopA(iStep)                        = PopM(iStep) ./ MolFracs(iStep,2) .* MolFracs(iStep,1);
    rhoI(1:NLevels(iBinnedMol),1)      = rhoM(iStep) .* Popp' ./ PopM(iStep);
    
    TempVec(1:NLevels(iBinnedMol),1)   = KDisss(LevelToBin(1:NLevels(iBinnedMol),iBinnedMol)) .* PopA(iStep) .* Popp(1:NLevels(iBinnedMol),1) .* (rhoA(iStep).^2 - rhoI(1:NLevels(iBinnedMol)));
    %TempVec(1:NBins(iBinnedMol),1)     = KDisss(1:NBins(iBinnedMol)) .* PopA(iStep) .* Popp(1:NBins(iBinnedMol)) .* (rhoA(iStep).^2 - rhoM(iStep));
    %TempVec(1:NBins(iBinnedMol),1)     = KDisss(1:NBins(iBinnedMol)) .* Popp(1:NBins(iBinnedMol)) ./ PopM(iStep) .* (rhoI(1:NBins(iBinnedMol)) - rhoA(iStep).^2);

    Determ        = sum( TempVec(1:NLevels(iBinnedMol),1) );
    
    CDInt(jStep)  = sum( TempVec(1:NLevels(iBinnedMol),1) .* LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol) )    ./ Determ;

    CDVib(jStep)  = sum( TempVec(1:NLevels(iBinnedMol),1) .* LevelEeVVib0(1:NLevels(iBinnedMol),iBinnedMol) ) ./ Determ;

    CDRot(jStep)  = sum( TempVec(1:NLevels(iBinnedMol),1) .* LevelEeVRot(1:NLevels(iBinnedMol),iBinnedMol) )  ./ Determ;
    
    jStep = jStep + 1;
  end
  
  CDInt = CDInt ./ abs(LevelEeV(1));
  CDVib = CDVib ./ abs(LevelEeV(1));
  CDRot = CDRot ./ abs(LevelEeV(1));

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
end