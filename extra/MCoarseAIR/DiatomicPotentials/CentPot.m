function [Vc, dVc] = CentPot(r, jqn, iMol)
  
  global AtomMass Pair_to_Atoms
  
  mass  = [AtomMass(Pair_to_Atoms(iMol,1)), AtomMass(Pair_to_Atoms(iMol,2))];
  mu    = mass(1) * mass(2) / ( mass(1) + mass(2) );
  Vc_R2 = 1.d0 / mu * 0.5d0 * (jqn+0.5d0)^2;

  Vc    = Vc_R2 ./ r.^2;  
  dVc   = - 2.d0 .* Vc ./ r;  
  
end