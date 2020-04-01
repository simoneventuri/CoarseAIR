function [dV, ddV, dVdJ] = FindDerivativeDiatPot(r, jqn0, iMol, Leveljqn)

  global MoleculesName AtomMass Pair_to_Atoms
  
  [Ve, dV] = DiatPot(r, jqn0, iMol);
 
  coeff   = [1/12	-2/3	0	2/3	-1/12];
  Epsilon = 1.d-6;
  Mask    = [-2,-1,0,1,2];
  for iR = 1:5
    rVec(iR)          = r + Epsilon * Mask(iR);
    [Ve(iR), dVe(iR)] = DiatPot(rVec(iR), jqn0, iMol);
  end
  ddV      = ( coeff(1) * dVe(1) + coeff(2) * dVe(2) + coeff(4) * dVe(4) + coeff(5) * dVe(5) ) / Epsilon;
  
  
  Epsilon = 1;
  if jqn0 > 1 && jqn0 < max(Leveljqn(:,iMol)-1)
    coeff   = [1/12	-2/3	0	2/3	-1/12];
    Mask    = [jqn0-2, jqn0-1, jqn0, jqn0+1, jqn0+2];
  elseif jqn0 > 3
    coeff   = [3/12	-16/12	+36/12	-48/12	+25/12];
    Mask    = [jqn0-4, jqn0-3, jqn0-2, jqn0-1, jqn0];
  else
    coeff   = [-25/12	48/12	-36/12	16/12	-3/12];
    Mask    = [jqn0, jqn0+1, jqn0+2, jqn0+3, jqn0+4];
  end
  
  for iJ = 1:5
    [Ve(iJ), dVe(iJ)] = DiatPot(r, Mask(iJ), iMol);
  end
  dVdJ     = ( coeff(1) * Ve(1) + coeff(2) * Ve(2) +  coeff(3) * Ve(3) + coeff(4) * Ve(4) + coeff(5) * Ve(5) ) / Epsilon;
  
end