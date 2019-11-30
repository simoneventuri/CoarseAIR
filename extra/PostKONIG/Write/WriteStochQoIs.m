function WriteStochQoIs(StochT, StochP, KDiss, PDF_KDiss, KDissQSS, PDF_KDissQSS, KExch, PDF_KExch, KExchQSS, PDF_KExchQSS, ...
                        tauIntPStoch_WExch,  PDF_tauIntP_WExch,  tauVibPStoch_WExch,  PDF_tauVibP_WExch,  tauRotPStoch_WExch,  PDF_tauRotP_WExch)
  
  NSamples = size(KDiss,2);
                      
  FileName = strcat('./StochasticQoIs_', num2str(StochT), '.csv');
  fileID1  = fopen(FileName,'w');
  
  fprintf(fileID1,'# T=%e K; P=%e atm\n', StochT, StochP)
  
  fprintf(fileID1,'#   MEAN: KDissQSS, KDissThermal, KExchQSS, KExchThermal, TauVib, TauRot, TauInt\n');
  fprintf(fileID1,'%e, %e, %e, %e, %e, %e, %e\n', PDF_KDissQSS.mu, PDF_KDiss.mu, PDF_KExchQSS.mu, PDF_KExch.mu, PDF_tauVibP_WExch.mu, PDF_tauRotP_WExch.mu, PDF_tauIntP_WExch.mu );
  
  fprintf(fileID1,'# ST DEV: KDissQSS, KDissThermal, KExchQSS, KExchThermal, TauVib, TauRot, TauInt\n');
  fprintf(fileID1,'%e, %e, %e, %e, %e, %e, %e\n', PDF_KDissQSS.sigma, PDF_KDiss.sigma, PDF_KExchQSS.sigma, PDF_KExch.sigma, PDF_tauVibP_WExch.sigma, PDF_tauRotP_WExch.sigma, PDF_tauIntP_WExch.sigma );
  
  fprintf(fileID1,'# iPES, KDissQSS, KDissThermal, KExchQSS, KExchThermal, TauVib, TauRot, TauInt\n');
  for i=1:NSamples
    fprintf(fileID1,'%i, %e, %e, %e, %e, %e, %e, %e\n', i, KDissQSS(i), KDiss(i), KExchQSS(i), KExch(i), tauVibPStoch_WExch(i), tauRotPStoch_WExch(i), tauIntPStoch_WExch(i) );
  end
  
  fclose(fileID1);
  
end