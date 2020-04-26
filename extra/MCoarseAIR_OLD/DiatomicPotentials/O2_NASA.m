function [V, dV] = O2_NASA(r)
  
  Zero = 0.d0;
  One  = 1.d0;
  Six  = 6.d0;
  
  nfitco  = 6;
  coefo2  = [-187.888473490235d0, -47.6881525811959d0, 13.8451858338530d0, -42.2149670400162d0,  15.9948962333531d0, -35.3271903135257d0, 32.0805760535816d0]; 
  c6o2    = 97.3185530097047d0;
  damp6co = 418.892d0;
  ao2     = 5.68d0;
  R0      = 2.15d0;
  alphaco = -2.234d0;
  npowco  = 9;

  vsr   = 64.0d0 .* exp(alphaco .* r) ./ r;
  dvsr  = vsr .* (alphaco - (One ./ r));
  vlr   = -c6o2  ./ (r.^6  + damp6co);
  dvlr  = -6.0d0 .*  r.^5 .* vlr ./ (r.^6 + damp6co);
  fact  = r.^npowco .* exp(-ao2 .* r);
  dfact = fact .* ((npowco ./ r) - ao2);
  
  cfcn  = coefo2(nfitco+1);
  dcfcn = Zero;
  for k = nfitco:-1:2
   dcfcn = cfcn      + dcfcn .* (r-R0);
   cfcn  = coefo2(k) + cfcn  .* (r-R0);
  end
  dcfcn = dfact     .* cfcn  + fact  .* dcfcn;
  cfcn  = coefo2(1)  + fact .* cfcn;
  
  V  = cfcn  + vsr  + vlr  + 187.888473490235d0;
  dV = dcfcn + dvlr + dvsr;

end