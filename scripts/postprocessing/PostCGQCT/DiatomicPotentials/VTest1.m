function [V, dV] = VTest2(r)
  
  V  = - exp( - abs(r) );
  dV = r .* exp( - abs(r) ) ./ abs(r);
 
end