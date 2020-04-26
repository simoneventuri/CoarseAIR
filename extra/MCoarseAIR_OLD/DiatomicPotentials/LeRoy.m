function [V, dV] = LeRoy(r)
  
  nrep  = 12; 
  co    = 49.d0;
  ao    = 2.4d0;
  aOne  = 5.38484562702061d0;
  c     = 23.9717901220746d0;
  cT    = c * 20.d0;
  cF    = c * 500.d0;
  d     = 3.d0;
  re    = 2.1d0;
  csr   = [ -3.89327515161896d2, -6.19410598346194d2, -5.51461129947346d2, -3.54896837660797d2, -1.08347448451266d2, -6.59348244094835d1, ...
             1.30457802135760d1,  7.20557758084161d1, -1.81062803146583d1, -2.84137057950277d1,  1.40509401686544d1, -1.84651681798865d0];

%   suma    =   csr(nrep);
%   dr      =   r - re;
%   r       =   ri;
%   
%   for i = nrep-1:-1:1
%     suma  =   csr(i) + suma * dr;
%   end
%   
%   Vrep = 0.d0;
%   Vdlr = 0.d0;
%   Vsr  = 0.d0;
%   
%   Vrep    =   co * exp(-ao * r) / r;
%   Vdlr    =   - c / ( r^6 + d^6 ) - cT  / ( r^4 + d^4 )^2 - cF / ( r^2 + d^2 )^5;
%   Vsr     =   suma * exp( - aOne * r ) * r^6;
%   V       =   Vrep + Vdlr + Vsr;
  

  dV      =   r - re;
  ri      =   1.d0 ./ r;
  
  Vrep    =   co .* exp(-ao .* r) ./ r;
  drep    = - Vrep .* ( ao + ri );
  
  for iR = 1:length(r)
    suma(iR)    =   csr(nrep);
    dsuma(iR)   =   0.d0;
    for n = nrep-1:-1:1
      dsuma(iR) =   suma(iR) + dV(iR) .* dsuma(iR);
      suma(iR)  =   csr(n) + suma(iR) .* dV(iR);
    end
  end

  ex      =   exp( - aOne .* r ) .* r.^6;
  Vsr     =   suma .* ex;

  drep    =   drep + ( (-aOne + (6.0d0 .* ri) ) .* suma + dsuma ) .* ex;

  Vdlr    =        - c  ./ ( r.^6 + d.^6 );
  Vdlr    =   Vdlr - cT ./ ( r.^4 + d.^4 ).^2;
  Vdlr    =   Vdlr - cF ./ ( r.^2 + d.^2 ).^5;

  dV      =   drep + 6.0  .* c  .* r.^5 ./ ( r.^6 + d.^6 ).^2;
  dV      =   dV   + 8.0  .* cT .* r.^3 ./ ( r.^4 + d.^4 ).^3;
  dV      =   (dV   + 10.0 .* cF .* r    ./ ( r.^2 + d.^2 ).^6) * 27.2113839712790;
  
  V       =   (Vrep + Vdlr + Vsr) .* 27.2113839712790;

end