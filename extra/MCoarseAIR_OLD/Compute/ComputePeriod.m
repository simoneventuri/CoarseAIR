function [Tau] = ComputePeriod(EProp, hE, rIn, rMin, rMax, VMaxEh, jqn, iMol)
  
  global MoleculesName AtomMass Pair_to_Atoms mu

  EP = EProp + hE;
  if ( EP > VMaxEh )
    hE = 0.5 * (VMaxEh - EProp);
    EP = EProposed + hE;
  end
  EM = EProp - hE;

  Range1 = [0.5*rIn, rMin];
  Range2 = [rMin,    rMax];

  [rInP, rOutP, ErrorFlg] = TurningPoints( Range1, Range2, EP*27.2113839712790, jqn, iMol);
  AP = ActionIntegral( rInP, rOutP, EP, jqn, iMol);

  [rInM, rOutM, ErrorFlg] = TurningPoints( Range1, Range2, EM*27.2113839712790, jqn, iMol);
  AM = ActionIntegral( rInM, rOutM, EM, jqn, iMol);

  Tau = 0.5 * ( AP - AM ) / hE;

end
