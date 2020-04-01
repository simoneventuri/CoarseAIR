function [V, dV] = NO_UMN(R)    

% NO UMN Min @ 2.17661 (V=-6.617426)

  B_To_Ang       = 0.52917721067;
  Kcm_To_Hartree = 0.159360144e-2;
  Hartree_To_eV  = 27.2113839712790;
 
  cs   = [ 0.322338e0, 5.878590e0, -12.790761e0, 13.320811e0, -7.516309e0, 1.875839e0, -0.052723e0, -0.037783e0, 0.48294e0, 1.98697e0]; 
  red  = 1.1508;
  de   = 152.6;
  VRef = 0.0;%0.191619504727d0;
  
  RAng = R .* B_To_Ang;
 
  u    = exp( -(RAng-red) ./ cs(9) - (RAng-red) .^2.0 ./ (cs(10)) );

  dfdr =  ( -2.0 .* (RAng-red)./cs(10) - 1.0./cs(9) );

  V    =  - de*(cs(1) .* u + cs(2) .* u.^2 + cs(3) .* u.^3 + cs(4) .* u.^4 + cs(5) .* u.^5.0 + cs(6) .* u.^6 + cs(7) .* u.^7 + cs(8) .* u.^8);
  V    = (V' * 0.159360144e-2 + VRef).*Hartree_To_eV;
  
  dV   =  - de .* (cs(1) .*          dfdr .* u        + ...
                   cs(2) .* 2.0d0 .* dfdr .* u.^2     + ...
                   cs(3) .* 3.0d0 .* dfdr .* u.^3     + ...
                   cs(4) .* 4.0d0 .* dfdr .* u.^4     + ...
                   cs(5) .* 5.0d0 .* dfdr .* u.^5     + ...
                   cs(6) .* 6.0d0 .* dfdr .* u.^6     + ...
                   cs(7) .* 7.0d0 .* dfdr .* u.^7     + ...
                   cs(8) .* 8.0d0 .* dfdr .* u.^8);
  
  dV = dV' .* Kcm_To_Hartree.*Hartree_To_eV .* B_To_Ang;

end