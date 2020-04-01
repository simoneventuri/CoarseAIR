function [V, dV] = N2_UMN(R)    

% N2 UMN Min @ 2.074067 (V=-9.90436)

  B_To_Ang       = 0.52917721092;
  Kcm_To_Hartree = 0.159360144e-2;
  Hartree_To_eV  = 27.2113839712790;
  
  %% FOR N2+O2
  cs   = [2.71405774451e0, 1.32757649829e-1, 2.66756890408e-1, 1.95350725241e-1, -4.08663480982e-1, 3.92451705557e-1, 1.13006674877e0];
  de   = 228.4;
  
  %% FOR N4
%   cs   = [2.70963254293, 1.32620177271e-1, 2.96757048793e-1, 1.97112432229e-1, -5.02002309588e-1, 3.80734244606e-1, 1.21001628750];
%   de   = 228.7;
  red  = 1.098;
  VRef = 0.0;%0.191619504727;
  
  red4   = red   .^ 4;
  RAng   = R     .* B_To_Ang;
  RAng3  = RAng  .^3;
  RAng4  = RAng3 .* RAng;
  
  TempSum = (RAng4 + red.^4);
  y       = (RAng4 - red.^4) ./ TempSum;
  y2      = y.^2;
  y3      = y2.*y;
  y4      = y2.*y2;
  y5      = y3.*y2;
  y6      = y3.*y3;
  
  fy      =   cs(1) + cs(2)*y + cs(3) .* y2 + cs(4) .* y3 + cs(5) .* y4 + cs(6) .* y5 + cs(7) .* y6;
  
  u       =   exp(-fy .* (RAng-red));
  minu    =   (1.0-u);

  dfdy    =   cs(2) + 2.0 .* cs(3) .* y + 3.0 .* cs(4) .* y2 + 4.0 .* cs(5) .* y3 + 5.0 .* cs(6) .* y4 + 6.0 .* cs(7) .* y5;
  dydr    =    8.0 .* RAng3 .* red4 ./ (TempSum.^2);
  dfdr    =   dfdy .* dydr;

  V       =   de .* minu.^2 - de;
  V       =   V .* Kcm_To_Hartree .* Hartree_To_eV + VRef; 
  
  dV      =   2.0 .* de .* minu .* u .* (dfdr .* (RAng - red) + fy);
  dV      =   dV .* Kcm_To_Hartree .* Hartree_To_eV .* B_To_Ang;
  
end