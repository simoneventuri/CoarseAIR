function [LevToBin, BinEnergy, BinEnergyMax, qnToBin] = ComputeLevToBin(NLevels, Levelvqn, Leveljqn, LevelEeV, DeltaEintDiss, DeltaEintDepth)

  global SystemPath BinnedMolName NBinnedMol Ue T0_Vec UKb NBinsNew NBinsNewQB NBinsNewBLow MinEnQB NBinsNewBBack jqnBack
  
  BinEnergy = 0.d0 * Levelvqn;
  LevToBin  = 1 + 0.d0 * Levelvqn;
  
  for iMol = 1:1%NBinnedMol
    
    NBinsNew(iMol)
    if NBinsNew(iMol) ==  max(Levelvqn(:,iMol)) + 1
      
      for iv = 0:max(Levelvqn(:,iMol))
        for iLevel = 1:NLevels(iMol)
          if Levelvqn(iLevel,iMol) == iv && Leveljqn(iLevel,iMol) == 0
            EnTemp(iv+1) = LevelEeV(iLevel,iMol);
          end 
        end
      end
      
      for iBin = 1:NBinsNew(iMol)-1
        DeltaVDissBin(iBin) = EnTemp(iBin+1) - 1.d-10;
      end
      DeltaVDissBin(NBinsNew(iMol)) = 0.d0
      
      NBinsNewBBack(iMol) = 0;
      jqnBack(iMol)       = 300;
      
    else
      
      NBinsNewBHigh(iMol) = NBinsNew(iMol) - (NBinsNewQB(iMol) + NBinsNewBLow(iMol) + NBinsNewBBack(iMol))
      
      for iv = 1:NBinsNewBLow(iMol)
        for iLevel = 1:NLevels(iMol)
          if Levelvqn(iLevel,iMol) == iv && Leveljqn(iLevel,iMol) == 0
            EnTemp(iv) = LevelEeV(iLevel,iMol)-1.d-10;
          end 
        end
      end

      hEnQB = MinEnQB(iMol) / NBinsNewQB(iMol);
      EnTemp(NBinsNew(iMol)-NBinsNewBBack(iMol)) = 0.d0;
      for iEeV = NBinsNew(iMol)-1-NBinsNewBBack(iMol):-1:(NBinsNew(iMol)-NBinsNewBBack(iMol)-NBinsNewQB(iMol))
        EnTemp(iEeV) = EnTemp(iEeV+1) + hEnQB;
      end
      
      hEnBHigh = (MinEnQB(iMol) - EnTemp(NBinsNewBLow(iMol))) / NBinsNewBHigh(iMol);
      for iEeV = NBinsNew(iMol)-NBinsNewBBack(iMol)-NBinsNewQB(iMol)-1:-1:NBinsNewBLow(iMol)+1
        EnTemp(iEeV) = EnTemp(iEeV+1) - hEnBHigh;
      end
      
      for iBin = 1:NBinsNew(iMol)-NBinsNewBBack(iMol)
        DeltaVDissBin(iBin) = EnTemp(iBin);
      end
      DeltaVDissBin
      
    end 
    
    BinEnergyMax(1:NBinsNew(iMol)-NBinsNewBBack(iMol),iMol) = EnTemp(:);
      
    [Temp, BinsOrder] = sort(rand(NBinsNew(iMol),1));
    %BinsOrder = [1:NBinsNew(iMol)];
    
    for iLevels = 1:NLevels(iMol)
      if Leveljqn(iLevels,iMol) < jqnBack(iMol) || DeltaEintDiss(iLevels,iMol) >= DeltaVDissBin(NBinsNew(iMol)-NBinsNewBBack(iMol)-NBinsNewQB(iMol))
        iBin = 1;
        while DeltaEintDiss(iLevels,iMol) >= DeltaVDissBin(iBin)
          iBin = iBin + 1;
        end
      else
        iBin = NBinsNew(iMol);
      end 
      LevToBin(iLevels,iMol)                 = BinsOrder(iBin); 
      BinEnergy(LevToBin(iLevels,iMol),iMol) = BinEnergy(LevToBin(iLevels,iMol),iMol) + exp( LevelEeV(iLevels,iMol) .* Ue ./ (T0_Vec(1) .* UKb) );
    end

%     for iLevels = 1:NLevels(iMol)
%       if Leveljqn(iLevels,iMol) < jqnBack(iMol) || DeltaEintDiss(iLevels,iMol) < abs(DeltaVDissBin(NBinsNew(iMol)-NBinsNewBBack(iMol)-NBinsNewQB(iMol)))
%         iBin = 1;
%         while DeltaEintDepth(iLevels,iMol) < abs(DeltaVDissBin(iBin))
%           iBin = iBin + 1;
%         end
%       else
%         iBin = NBinsNew(iMol);
%       end 
%       LevToBin(iLevels,iMol)                 = BinsOrder(iBin); 
%       BinEnergy(LevToBin(iLevels,iMol),iMol) = BinEnergy(LevToBin(iLevels,iMol),iMol) + exp( LevelEeV(iLevels,iMol) .* Ue ./ (T0_Vec(1) .* UKb) );
%     end
    
    clear EnTemp BinsOrder DeltaVDissBin
    
% 
%     [Temp, BinsOrder] = sort(rand(NBinsNew(iMol),1));
% %     %BinsOrder = [1:NBinsNew(iMol)];
%     for iLevels = 1:NLevels(iMol)
%       iBin = 1;
%       if Levelvqn(iLevels,iMol) <= 0 %&& DeltaEintDiss(iLevels,iMol) < -8.d0
%         iBin = 2;
%       end
%       LevToBin(iLevels,iMol)                 = BinsOrder(iBin); 
%       BinEnergy(LevToBin(iLevels,iMol),iMol) = BinEnergy(LevToBin(iLevels,iMol),iMol) + exp( LevelEeV(iLevels,iMol) .* Ue ./ (T0_Vec(1) .* UKb) );
%     end
%     BinEnergyMax = BinEnergy .* 0.d0;
%      

    for i = 1:NLevels(iMol)
     qnToBin(Levelvqn(i,iMol)+1,Leveljqn(i,iMol)+1) = LevToBin(i,iMol);
    end

    fid=fopen(['./' strcat(BinnedMolName(iMol,:), '_LevToBin') strcat('.dat')],'wt');
    fprintf(fid,'# Levels To Bins\n');

    for i = 1:NLevels(iMol)
     fprintf(fid,'%g\n',LevToBin(i,iMol));
    end
    
  end
  
end