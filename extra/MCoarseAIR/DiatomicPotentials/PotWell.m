function [V, dV] = VTest1(r)

  % Potential Well with Infinitely Repulsive Walls
  
  if abs(r) > 1
    V = 1e200;
  else
    V = 0.0;
  end
  dV = 0;

end