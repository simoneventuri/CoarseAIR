function [ErriE, NChngsSign] = ComputeNumerov(PsiRefi, PsiRefii, PsiReff, R, dR, EGuess, jqn, iMol, PlotFlg)

  global AtomMass MoleculesName Pair_to_Atoms mu

  Psi(1) = PsiRefi;
  Psi(2) = PsiRefii;

  %% Computing Diatomic Energy and Shifting it
  %VDiat = PotWell(R(1:2));
  [VDiat, dV] = DiatPot(R, jqn, iMol);
  VDiat = min(VDiat./27.2113839712790, 1e300);
  ff    = 2.0 * mu * (VDiat - EGuess);

  %% Initializing Wave Function at Left Bound
  Phi(1:2) = Psi(1:2) .* (1.0 - dR^2 .* ff(1:2) ./ 12.0);
  
  NChngsSign = 0;
  %% Starting Exploring Spatial Grid
  for iR=2:length(R)-1

    %% Computing Numerov Function @ k+1
    Phi(iR+1) = 2.0 .* Phi(iR) - Phi(iR-1) + dR^2 .* ff(iR) .* Psi(iR);

    %% Computing Wave Function @ k+1
    Psi(iR+1) = Phi(iR+1) ./ (1.0 - dR^2 .* ff(iR+1) ./ 12.0);
    
    %% Finding out how Many Times Wave Function Changes Sign
    if (Psi(iR+1)*Psi(iR) <= 0.0) && (Psi(iR) ~= 0)
      NChngsSign = NChngsSign + 1;
    end
  end
  
  %% Computing Error at the Right Bound
  ErriE   = Psi(end) - PsiReff;

  PsiTot  = trapz(R,Psi);
  Psi     = Psi ./ PsiTot;
  if (PlotFlg)
    figure(1)
    plot(R, Psi, 'linewidth', 1)
  end
  
end