function [EPred] = ComputeBNNPES(R, PreLogShift, Params)

  global System AbscissaConverter
  
  NData = size(R,1);
  
  p1 = ComputeBondOrder(R(:,1), Params.Lambda, Params.re);
  p2 = ComputeBondOrder(R(:,2), Params.Lambda, Params.re);
  p3 = ComputeBondOrder(R(:,3), Params.Lambda, Params.re);
     
  G  = ComputePIP(p1, p2, p3, Params.G_MEAN, Params.G_SD);

  b1Mat = repmat(Params.b1',[NData,1]);
  b2Mat = repmat(Params.b2',[NData,1]);
  b3Mat = repmat(Params.b3',[NData,1]);

  z1     = G * Params.W1 + b1Mat;
  y1     = tanh(z1);
  z2     = y1 * Params.W2 + b2Mat;
  y2     = tanh(z2);
  EPred  = y2 * Params.W3 + b3Mat;

  EPred  = exp(EPred + Params.Noise);
    
  EPred  = EPred - PreLogShift;
  
  if strcmp(System,'N3')
    [E1, dE1] = N2_LeRoy(R(:,1)'./AbscissaConverter);
    [E2, dE2] = N2_LeRoy(R(:,2)'./AbscissaConverter);
    [E3, dE3] = N2_LeRoy(R(:,3)'./AbscissaConverter);
  elseif strcmp(System,'O3')
    [E1, dE1] = O2_UMN(R(:,1)'./AbscissaConverter);
    [E2, dE2] = O2_UMN(R(:,2)'./AbscissaConverter);
    [E3, dE3] = O2_UMN(R(:,3)'./AbscissaConverter);
  end
  EDiat = E1 + E2 + E3;
  
  EPred = EPred + EDiat;

end