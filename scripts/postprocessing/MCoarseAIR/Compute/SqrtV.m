function y = SqrtV(x, Eint, jqn, iMol)

  global mu MoleculesName Pair_to_Atoms
  
  r = x;
      
  [Ve, dVe] = DiatPot(x, jqn, iMol);

  y  =  sqrt(abs(Eint - Ve./27.2113839712790));

end