function [qssInt, qssRot, qssVib] = ComputeQssDistr(NLevels, LevelQ, Levelg, Pop, QBins, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelvqn, LevelToBin, Steps)    
    
  global NBins BinnedMolName NBinnedMol KeV KinMthd

  vPop = 1.d-100 * ones(size(Pop,1), max(max(Levelvqn))+1, NBinnedMol);
%  jPop = 1.d-100 * ones(size(Pop,1), max(max(Levelvqn))+1, NBinnedMol);
      
  for iBinnedMol = 1:NBinnedMol

    for iSteps = Steps
      iSteps
      
     LevelPop(iSteps,1:NLevels(iBinnedMol),iBinnedMol) = Pop(iSteps,LevelToBin(1:NLevels(iBinnedMol),iBinnedMol),iBinnedMol)';

      for iLevels = 1:NLevels(iBinnedMol)
        vPop(iSteps,Levelvqn(iLevels,iBinnedMol)+1,iBinnedMol) = vPop(iSteps,Levelvqn(iLevels,iBinnedMol)+1,iBinnedMol) + LevelPop(iSteps,iLevels,iBinnedMol);
      end      

      PopTot(iSteps,iBinnedMol) = sum(LevelPop(iSteps,1:NLevels(iBinnedMol),iBinnedMol));
      
      
      qssInt(iSteps,iBinnedMol) = LevelPop(iSteps,1:NLevels(iBinnedMol),iBinnedMol) ./ PopTot(iSteps,iBinnedMol);

      qssRot = 0.d0;
%      eRot(iSteps,iBinnedMol) = sum( LevelEeVRot(1:NLevels(iBinnedMol),iBinnedMol) .* LevelPop(iSteps,1:NLevels(iBinnedMol),iBinnedMol)' ) ./ PopTot(iSteps,iBinnedMol);


      qssVib(iSteps,iBinnedMol) = vPop(iSteps,1:max(Levelvqn(:,iBinnedMol))+1,iBinnedMol) ./ PopTot(iSteps,iBinnedMol);
      
    end

  end
  
  clear vPop PopMol LevelPop
end