function [V, dV] = MorsePot(r)

  De   = 9.385;
  Beta = 2.799 / 1.88973;
  re   = 1.105 * 1.88973;
  
  V  = De * ( (1.0 - exp( -Beta*(r-re)) ).^2 - 1.0);
  dV = 2.0 * Beta * De * exp( -Beta*(r-re)) .* (1.0 - exp( -Beta*(r-re)));

end