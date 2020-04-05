function [dAdE, P, Gamma, Life] = ComputeTunneling(NLevels, Levelvqn, Leveljqn, LevelEh, VMax, rIn, rOut)

  global AtomMass MoleculesName Pair_to_Atoms EhToeV NMolecules
  
  EpsE  = 1.d-4;
  rMin  = 1.d0;
  rMax  = 30.d0;
  rHalf = (rMax + rMin) ./ 2.d0;
  
  for iMol = 1:NMolecules
  
    for iLevels = 1:NLevels(iMol)
      if LevelEh(iLevels,iMol) > 0.d0
        
        iLevels
        
        ErrorFlg = 1;
        rga  = rOut(iLevels,iMol)+2.d0
        rgb  = rMax * 1.d0
        pause
        Iter = 1;
        while ErrorFlg == 1 && Iter < 10
          rga = rga + 1.d0;
          rgb = rga * 2.d0;
          [rOutTemp, r3(iLevels,iMol), ErrorFlg] = TurningPoints( [rOut(iLevels,iMol)-1.d0, rMax], [rOut(iLevels,iMol)+1.d0, rMax.*2], LevelEh(iLevels,iMol), Leveljqn(iLevels,iMol), iMol);
          Iter = Iter + 1;
        end
        if ErrorFlg == 1
          stop
        end
        
        dE = min( (VMax(iLevels,iMol) ./ EhToeV - LevelEh(iLevels,iMol)), EpsE )
        Em = LevelEh(iLevels,iMol) - dE;
        Ep = LevelEh(iLevels,iMol) + dE;

        [rOutm, r3m, ErrorFlg] = TurningPoints( [rMin, rHalf], [rHalf, rMax], Em, Leveljqn(iLevels,iMol), iMol);
        [Action3m]             = ActionIntegral(rOutm, r3m, Em, Leveljqn(iLevels,iMol), iMol);

        [rOutp, r3p, ErrorFlg] = TurningPoints( [rMin, rHalf], [rHalf, rMax], Ep, Leveljqn(iLevels,iMol), iMol);
        [Action3p]             = ActionIntegral(rOutp, r3p, Ep, Leveljqn(iLevels,iMol), iMol);


        dAdE = (Action3p - Action3m) ./ (2.d0 * dE );

        [Action3]  = ActionIntegral(rOut, r3, LevelEh(iLevels,iMol), Leveljqn(iLevels,iMol), iMol);

        P(iLevel,iMol) = exp( - Action3);

        Gamma(iLevel,iMol) = P(iLevel,iMol) / dAdE(iLevel,iMol);

        Life(iLevel,iMol) = 1.d0 ./ Gamma(iLevel,iMol);

      end
    end 
  end
  
end