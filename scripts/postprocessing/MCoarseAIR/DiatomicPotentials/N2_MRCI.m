function [V, dV] = N2_MRCI(R)

% N2 Min @ 2.088828 (V=-9.3437497)

  csr  = [-7.64388144433896e+02, -2.11569028013955e+03, -2.67675739355454e+03, -2.46303302042839e+03, -1.53331085357664e+03, -6.03848288404665e+02, 1.27445322988364e+02, 1.49427946215489e+02];
  co   = 49.0;
  ao   = 2.4;
  a1   = 6.0;
  clr  = 21.8526654564325;
  d    = 2.5;
  nrep = 8;
  nlr  = 10;
  re   = 2.1;
  
  c6   = clr;
  c8   = clr*20.0;
  c10  = clr*500.0;

  vrep =    co .* exp(-ao.*R) ./ R;
  dvsr = -vrep .* (ao + (1.0./R));

  suma  = csr(nrep);
  dsuma = 0.0;
  dr    = R - re;
  for n=nrep-1:-1:1
      dsuma = suma   + dsuma .* dr;
      suma  = csr(n) + suma  .* dr;
  end
  
  ex   = exp(-a1.*R) .* (R.^6);
  vsr  = suma .* ex;
  drep = ((-a1 + (6.0 ./R)) .*suma + dsuma) .* ex;

  vdlr =      - c6  ./ (R.^6 + d.^6);
  vdlr = vdlr - c8  ./ (R.^4 + d.^4).^2;
  vdlr = vdlr - c10 ./ (R.^2 + d.^2).^5;
  
  dvdlr =       +  6.0 .*  c6 .* R.^5 ./ (R.^6 + d.^6).^2;
  dvdlr = dvdlr +  8.0 .*  c8 .* R.^3 ./ (R.^4 + d.^4).^3;
  dvdlr = dvdlr + 10.0 .* c10 .* R    ./ (R.^2 + d.^2).^6;
  
  V    = (vrep + vdlr + vsr)'    .* 27.2113839712790;
  dV   = (drep + dvdlr +  dvsr)' .* 27.2113839712790;
  
end