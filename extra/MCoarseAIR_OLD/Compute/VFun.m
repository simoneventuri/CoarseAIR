function [y, yJ] = VFun(x, jqn, Eint, iMol)

  global AtomMass MoleculesName Pair_to_Atoms
  
  [y, dVe]  = DiatPot(x, jqn, iMol);
  
  [yJ, dVc] = CentPot(x, jqn, iMol);
  
  y = Eint - y;
  
end