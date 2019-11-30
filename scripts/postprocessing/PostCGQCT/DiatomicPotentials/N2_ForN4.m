function [V, dV] = N2_ForN4(r)

% N2 Min 

  idp  = 10;
  idq  = 5;
  idd  = 5;
  npow = 6;
  ao   = 2.0;
  a1   = 6.0;
  a0   = 2.0;
  einf = -109.036558873442;
  c6   = -58.893858362275;
  d    = 2.5;
  co   = 49;
  pes  = [146.736533087545,101.160591312956,-251.323812509419,-499.020252425371,-210.074594675137, 135.174682961986,214.7693593595,425.317131479472,570.993394092895,266.199319038817];
  a2   = 2.0;
  qm   = [17.3398942423333,10.2017829193229,-4.7194296899574,-9.61643152091987,-2.12839249419926];
  a3   = 2.0;
  dip  = [-42.6175278814367,873.328188249918,984.492872053711,518.057930455541,124.185667297002];
  dp1   = (873.238*exp(-2.2)*2.1 - 42.6175);
  dp2   = dp1;
  preqq = 1.5*sqrt(2.0);
  re     =  2.1;
  re2    =  2.0*re; 

  iwant = [1, 0, 0]
  V1  = 0.0
  V2  = 0.0
  V3  = 0.0
  dV1 = 0.0
  dV2 = 0.0
  dV3 = 0.0
  
  
  if (iwant(1) == 1)

    % csd diatomic potential
    r2 = r  .* r;
    r4 = r2 .* r2;
    r6 = r2 .* r4;
    
    d2 = d  .* d;
    d4 = d2 .* d2;
    d6 = d2 .* d4;
    
    c8  = c6 .* 20.0;
    c10 = c6 .* 500.0;
    vlr = einf + c6 ./ (r6+d6) + c8 ./ (r4+d4).^2 + c10 ./ (r2+d2).^5;
    
    vrep    = co .* exp( -ao .* r ) ./ r;
    dvlrdr  = -6.0  .* c6 * (r4.*r) ./ (r6+d6).^2 - 8.0 * c8 * (r2.*r) ./ (r4+d4).^3 - 10.0 .* c10 .* r ./ (r2+d2).^6;
    dvrepdr = -vrep .* ao - vrep ./ r;
    dr    = r-re;
    ex    = exp( -a1 .* r ) .* r.^npow;
    dexdr = -a1 .* ex + ex .* npow ./ r;
    poly  = pes(idp)
    for i=idp-1:-1:1
      poly = pes(i) + poly .* dr;
    end

    dpoly  = (((((((( 9.0  .* pes(10) .* dr + ...
                      8.0  .* pes(9)) .* dr + ...
                      7.0  .* pes(8)) .* dr + ...
                      6.0  .* pes(7)) .* dr + ...
                      5.0  .* pes(6)) .* dr + ...
                      4.0  .* pes(5)) .* dr + ...
                      3.0  .* pes(4)) .* dr + ...
                      2.0  .* pes(3)) .* dr + ...
                              pes(2))
    V1  = vrep + poly .* ex + vlr;
    dV1 = dvrepdr + dvlrdr + poly .* dexdr + ex .* dpoly;
  end

  
  if (iwant(2) == 1)
    % quadrupole
    ex  = exp( -a2 .* r );
    dr  = r-re;
    sum = qm(idq);
    for i=idq-1:-1:1
      sum = qm(i) + sum .* dr;
    end
    V2    = sum .* ex .* r;
    dV2 = sum .* ex - a2 .* V2 + r .* ex .* ( ( (4.0 .* qm(5) .* dr + 3.0 .* qm(4) ) .* dr + 2.0 .* qm(3) ) .* dr + qm(2) );
  end

  
  if (iwant(3) == 1)
    % dipole polarizability
    ex   = exp( -a3 .* r );
    dr   = r-re;
    sum  = dip(idd);
    dsum = 0.0;                                                 
    for i=idd-1:-1:2
      dsum =    sum + dsum .* dr;                                              
      sum  = dip(i)  + sum .* dr;
    end
    V3    = sum .* ex .* r + dip(1);
    dV3   = ex .* r .* dsum + sum .* ex - a3 .* r .* ex .* sum;
  end
  
  V       = (V1+V2+V3)'    .*27.2113839712790;
  dV      = (dV1+dV2+dV3)' .*27.2113839712790;
  
 
end