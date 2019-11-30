function [p] = ComputeBondOrder(R, Lambda, re)

  global BondOrderFun
  
  if strcmp(BondOrderFun,'MorseFun')
    p = exp( - Lambda(1) .* (R - re(1)) )';
  elseif strcmp(BondOrderFun(:),'GaussFun')
    p = exp( - Lambda(1) .* (R - re(1)).^2 )';
  elseif strcmp(BondOrderFun,'MEGFun')
    p = exp( - Lambda(1) .* (R - re(1)) - Lambda(2) .* (R - re(2)).^2 )';
  end
  
end
