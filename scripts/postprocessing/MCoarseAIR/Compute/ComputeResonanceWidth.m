function [Egam] = ComputeResonanceWidth( rMinProp, rMaxProp, EProp, jqn, Tau, iMol)
  
  global MoleculesName AtomMass Pair_to_Atoms mu
  
  rEpsResWidthGlobal = 0.01;
  
  Range1 = [rMinProp, rMaxProp];

  rTemp     = rMaxProp;
  [Ve, dVe] = DiatPot(rTemp, jqn, iMol);
  VTemp     = Ve./27.2113839712790;
  while (VTemp >= EProp)
    rTemp     = rTemp + rEpsResWidthGlobal;
    [Ve, dVe] = DiatPot(rTemp, jqn, iMol);
    VTemp     = Ve./27.2113839712790;
  end 
  
  Range2 = [rMaxProp, rTemp];
  %Range2 = [rMaxProp, 1.e2];
  
  [rIn, rOut, ErrorFlg] = TurningPoints( Range1, Range2, EProp*27.2113839712790, jqn, iMol);

%   AProp = ActionIntegral( rIn, rOut, EProp, jqn, iMol);
%   Egam  = exp( -AProp ) / Tau;

  AProp = ActionIntegralGamma( rIn, rOut, EProp, jqn, iMol);
  Egam  = AProp / Tau;

end
