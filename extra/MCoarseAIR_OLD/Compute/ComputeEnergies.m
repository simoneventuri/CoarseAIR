function [eInt, eRot, eVib] = ComputeEnergies(NLevels, LevelQ, Levelg, Pop, QBins, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelvqn, LevelToBin, Steps)    
    
  global NBins BinnedMolName NBinnedMol KeV KinMthd Ue T0_Vec UKb

  vPop = 1.d-100 * ones(size(Pop,1), max(max(Levelvqn))+1, NBinnedMol);
  eInt = zeros(size(Pop,1),1);
  eRot = zeros(size(Pop,1),1);
  eVib = zeros(size(Pop,1),1);
  
  for iBinnedMol = 1:NBinnedMol

    for iSteps = Steps
      iSteps
      
      if sum(KinMthd(1,:) == 'CGM') == 3 || sum(KinMthd(1,:) == 'VIB') == 3
        LevelPop(iSteps,1:NLevels(iBinnedMol),iBinnedMol) = Pop(iSteps,LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol)' ./ QBins(LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol) .* Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* exp( - LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol) .* Ue ./ (T0_Vec(1) .* UKb) );
      else
        LevelPop(iSteps,1:NLevels(iBinnedMol),iBinnedMol) = Pop(iSteps,LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol)';
      end 

      %LevelPop(iSteps,1:NLevels(iBinnedMol),iBinnedMol) = Pop(iSteps,LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol)';

      for iLevels = 1:NLevels(iBinnedMol)
        vPop(iSteps,Levelvqn(iLevels,iBinnedMol)+1,iBinnedMol) = vPop(iSteps,Levelvqn(iLevels,iBinnedMol)+1,iBinnedMol) + LevelPop(iSteps,iLevels,iBinnedMol);
      end      

      PopTot(iSteps,iBinnedMol) = sum(LevelPop(iSteps,1:NLevels(iBinnedMol),iBinnedMol));
      
      
      eInt(iSteps,iBinnedMol) = sum( LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol) .* LevelPop(iSteps,1:NLevels(iBinnedMol),iBinnedMol)' ) ./ PopTot(iSteps,iBinnedMol);


      eRot(iSteps,iBinnedMol) = sum( LevelEeVRot(1:NLevels(iBinnedMol),iBinnedMol) .* LevelPop(iSteps,1:NLevels(iBinnedMol),iBinnedMol)' ) ./ PopTot(iSteps,iBinnedMol);


      eVib(iSteps,iBinnedMol) = sum( vEeVVib0(1:max(Levelvqn(1:NLevels(iBinnedMol),iBinnedMol))+1,iBinnedMol) .* vPop(iSteps,1:max(Levelvqn(:,iBinnedMol))+1,iBinnedMol)' ) ./ PopTot(iSteps,iBinnedMol);
      
    end

  end
  
  clear vPop PopMol LevelPop
end