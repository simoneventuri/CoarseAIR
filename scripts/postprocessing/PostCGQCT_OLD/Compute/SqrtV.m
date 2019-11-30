function y = SqrtV(x, Eint, jqn, iMol)

  global AtomMass MoleculesName Pair_to_Atoms
  
  r = x;
    
  mass = [AtomMass(1), AtomMass(2)];
  mu = mass(1) * mass(2) / ( mass(1) + mass(2) );
  
  [Ve, dVe] = DiatPot(x, jqn, iMol);

  y  =  sqrt(abs(Eint - Ve));

end