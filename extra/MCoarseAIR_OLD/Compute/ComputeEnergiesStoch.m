function [eInt, eRot, eVib] = ComputeEnergiesStoch(NLevels, LevelQ, Levelg, PopStoch, QBins, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelvqn, LevelToBin, Steps)    
    
  global NBins BinnedMolName NBinnedMol KeV KinMthd Ue T0_Vec UKb iPESEnd iPESStart
  
  NPESTemp = iPESEnd - iPESStart + 1; 

  NSteps = length(Steps);
  
  vPop = 1.d-100 * ones(NSteps, max(max(Levelvqn))+1, NPESTemp);
  eInt = zeros(NSteps,NPESTemp);
  eRot = zeros(NSteps,NPESTemp);
  eVib = zeros(NSteps,NPESTemp);

  jSteps=1;
  for iSteps=Steps
    iSteps

    for iPES=1:NPESTemp
      
      LevelPop(jSteps,1:NLevels(1),iPES) = PopStoch(iSteps,LevelToBin(1:NLevels(1),1),iPES+iPESStart-1)' ./ QBins(LevelToBin(1:NLevels(1),1),1) .* Levelg(1:NLevels(1),1) .* exp( - LevelEeV0(1:NLevels(1),1) .* Ue ./ (T0_Vec(1) .* UKb) );

      for iLevels = 1:NLevels(1)
        vPop(jSteps,Levelvqn(iLevels,1)+1,iPES) = vPop(jSteps,Levelvqn(iLevels,1)+1,iPES) + LevelPop(jSteps,iLevels,iPES);
      end      

      PopTot(jSteps,iPES) = sum(LevelPop(jSteps,1:NLevels(1),iPES));


      eInt(jSteps,iPES) = sum( LevelEeV0(1:NLevels(1),1) .* LevelPop(jSteps,1:NLevels(1),iPES)' ) ./ PopTot(jSteps,iPES);


      eRot(jSteps,iPES) = sum( LevelEeVRot(1:NLevels(1),1) .* LevelPop(jSteps,1:NLevels(1),iPES)' ) ./ PopTot(jSteps,iPES);


      eVib(jSteps,iPES) = sum( vEeVVib0(1:max(Levelvqn(1:NLevels(1),1))+1,1) .* vPop(jSteps,1:max(Levelvqn(:,1))+1,iPES)' ) ./ PopTot(jSteps,iPES);

    end
      
    jSteps=jSteps+1;
  end
  
  clear vPop PopMol LevelPop
end