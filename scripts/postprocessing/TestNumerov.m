close all
clear all
clc

global System T0_Vec
global NBins NTint NAtoms AtomsName MoleculesName NMolecules DegeneracyFactor BinnedMolName ColPartToComp NBinnedMol BinnedMolToComp NComp ...
       CompNames CompColor AtomColor AtomSize AllMoleculesName PairColor AtomMass ComponentMass ColorVec ComponentDeg MoleculeMu RxLxIdx MoleculedDissEn
global BCVec Pair_to_Atoms mu

%% Choose one of the Systems Implemented in UpdateAllInput.m; Select H3 for a System of Morse Potentials 
System        = 'N2O'
T0_Vec        = [];
Pair_to_Atoms = [1,2; 1,3; 2,3];

UpdateAllInput()
AtomsName


dR     =  5.e-3;          %% Spatial Grid Spacing
Ri     =    1.0;          %% Where to Locate Left Boundary Condition

PsiRefi  = 0.0;           %% Value of Wave Function at (the Left Boundary Condition)
PsiRefii = 1.e-100;       %% Value of Wave Function at (the Left Boundary Condition + 1*dR)
PsiReff  = 1.e-200;        %% Value of Wave Function at (the Right Boundary Condition)

StartEnergy  =  [-1.e10];    %% Minimum Energy Required to the Levels
DeltaESearch =  1.e-5;    %% Step for Increasing Energy in the Level Finding 

EpsE         = 1.e-13;    %% Stopping Criterium based on Energy; Maximum Error in Energy Required to the Levels
%%EpsPsi       = 1.e-15;  %% Stopping Criterium based on Wave Function; Maximum Error in Wave Functions Required to the Levels
 
jqnMin       = [103];       %% Minimum Rotational Q.N. Required to the Levels
jqnMax       = [103];       %% Maximum Rotational Q.N. Required to the Levels
vqnMax       = [100];    %% Maximum Vibrational Q.N. Required to the Levels

PlotLevelsFlg = true;     %% Flag for 3D Plot of Levels on Top of Diatomic Potential
PlotPsiFlg    = true;     %% Flag for Plot of Wave Functions


for iMol=1:1
  mu = (AtomMass(Pair_to_Atoms(iMol,1)) * AtomMass(Pair_to_Atoms(iMol,2))) / (AtomMass(Pair_to_Atoms(iMol,1)) + AtomMass(Pair_to_Atoms(iMol,2)));
  MoleculesName(iMol,:)
  
  rMinTemp = [];
  VMinTemp = [];
  rMaxTemp = [];
  VMaxTemp = [];
  rMinOld  = 1.5;   %% Right Extreme for Potential Minima Finding
  rMaxOld  = 10.0;  %% Left Extreme  for Potential Maxima Finding

  iLevels = 0;
  for jqn = jqnMin(iMol):jqnMax(iMol)
    

    %% Computing Maxima and Minima of Diatomic Potential
    [rMinTemp(jqn+1), VMinTemp(jqn+1), rMaxTemp(jqn+1), VMaxTemp(jqn+1), VInf] = FindMinMaxDiatPot(jqn, iMol, rMinOld, rMaxOld);
    rMinOld = rMinTemp(jqn+1);
    rMaxOld = rMaxTemp(jqn+1);
    VMinTemp(jqn+1) = VMinTemp(jqn+1) ./27.2113839712790;
    Ri     = fzero(@(x) DiatPotShift(x, jqn, iMol, VMaxTemp(jqn+1)), 0.5);
    VMaxTemp(jqn+1) = VMaxTemp(jqn+1) ./27.2113839712790;
    Rf      = rMaxTemp(jqn+1);   %% Where we Locate Right Boundary Condition
    R       = [Ri-0.5:dR:Rf];
    vqn     = -1;

    %% Plotting Diatomic Potential
    if (PlotLevelsFlg)
      rMinPlot = 0.8;
      rMaxPlot = 30.0;
      RPlot    = linspace(rMinPlot, rMaxPlot, 5000);
      [VDiat, dV] = DiatPot(RPlot, jqn, iMol);
      figure(100)
      plot(RPlot, VDiat, '-k', 'LineWidth', 3)
      ylim([0.0 5.0])
      figure(1+iMol)
      hold on
      plot3(RPlot, jqn+RPlot.*0.0, VDiat, '-k', 'LineWidth', 3)
    end 
  
    
    if StartEnergy(iMol) > VMinTemp(jqn+1)
      ECurrent = StartEnergy(iMol);
    else
      ECurrent = VMinTemp(jqn+1);       
    end
    while (ECurrent <= VMaxTemp(jqn+1)) && (vqn <= vqnMax(iMol))
      
      if (PlotPsiFlg)
        fig = figure(1)
        hold on  
      end
      
      %% Finding Nb of Times the Wave Functions Changes Signs for ECurrent
      iE = 1;
      [ErriE, NChngsSigniECurrent] = ComputeNumerov(PsiRefi, PsiRefii, PsiReff, R, dR, ECurrent, jqn, iMol, false);
      Err(iE) = ErriE;
      NChngsSign(iE) = NChngsSigniECurrent;
 
      %% Making sure that the  Nb of Times the Wave Functions Changes Signs for NChngsSigniETry is Larger than NChngsSigniECurrent
      iE = 2;
      ETry            = ECurrent;
      NChngsSigniETry = NChngsSigniECurrent;
      %fprintf('%e,%e,%e\n', ECurrent, NChngsSigniECurrent, ErriE)
      while (NChngsSigniETry == NChngsSigniECurrent)
        if (VMaxTemp(jqn+1) - ECurrent)     > 1.0
          ETry                     = ETry + DeltaESearch*1.e1;
        elseif (VMaxTemp(jqn+1) - ECurrent) > 1.e-1
          ETry                     = ETry + DeltaESearch;
        elseif (VMaxTemp(jqn+1) - ECurrent) > 1.e-2
          ETry                     = ETry + DeltaESearch/1.e1;
        elseif (VMaxTemp(jqn+1) - ECurrent) > 1.e-3
          ETry                     = ETry + DeltaESearch/1.e2;
        elseif (VMaxTemp(jqn+1) - ECurrent) > 1.e-4
          ETry                     = ETry + DeltaESearch/1.e3;
        end
        [ErriE, NChngsSigniETry] = ComputeNumerov(PsiRefi, PsiRefii, PsiReff, R, dR, ETry, jqn, iMol, false);
        Err(iE)        = ErriE;
        NChngsSign(iE) = NChngsSigniETry;
        NChngsSign(iE);
        %fprintf('%e,%e,%e\n', ETry, NChngsSigniETry, ErriE)
      end
     
                
      %% %% BISECTION METHOD for ROOT-FINDING
      %% Check if 2 Initial Guesses for Energy are Correct
      EStart  = [ECurrent, ETry];
      if prod(Err(1:2)) > 0.0
        fprintf('ERRROR! Initial Guesses for Energies Bring to Psi(end) with Same Signs!')
        stop
      else
        if Err(1) > 0.0
          Ep = EStart(1);
          Em = EStart(2);
        else
          Ep = EStart(2);
          Em = EStart(1);
        end
      end
      %% Iterating
      while abs(Ep - Em) > EpsE   %% OR abs(Err(iE) > EpsPsi), if Stopping Criterium based on Wave Function
        iE      = iE+1;
        E       = (Ep + Em) ./ 2.0;
        [ErriE, NChngsSigniE] = ComputeNumerov(PsiRefi, PsiRefii, PsiReff, R, dR, E, jqn, iMol, false);
        Err(iE)        = ErriE;
        NChngsSign(iE) = NChngsSigniE;
        %% Updating Energy
        if Err(iE) > 0.0
          Ep = E;
        else
          Em = E;
        end
      end 
      EProposed = (Ep + Em) / 2.0;
      %% %%
      
      
      if EProposed < VMaxTemp(jqn+1)
        iLevels = iLevels + 1;
        %% Plotting Wave Function
        if (PlotPsiFlg)
          [ErriETemp, NChngsSigniETemp] = ComputeNumerov(PsiRefi, PsiRefii, PsiReff, R, dR, E, jqn, iMol, true);
        end

        %% Storing Final Value for Energy
        LevelEh(iLevels,iMol)  = EProposed;
        LevelEeV(iLevels,iMol) = EProposed .* 27.2113839712790;
        
        %% Storing Level Properties
        Levelvqn(iLevels,iMol) = vqn+1;
        Leveljqn(iLevels,iMol) = jqn;
        rMin(iLevels,iMol)     = rMinTemp(jqn+1);
        rMax(iLevels,iMol)     = rMaxTemp(jqn+1);
        VMineV(iLevels,iMol)   = VMinTemp(jqn+1) .* 27.2113839712790;
        VMaxeV(iLevels,iMol)   = VMaxTemp(jqn+1) .* 27.2113839712790;
        VMin(iLevels,iMol)     = VMinTemp(jqn+1);
        VMax(iLevels,iMol)     = VMaxTemp(jqn+1);

        
        %% Computing Inner and Outer Turning Points
        [rIn(iLevels,iMol), rOut(iLevels,iMol), ErrorFlg] = TurningPoints( [0.5d0, rMin(iLevels,iMol)], [rMin(iLevels,iMol), rMax(iLevels,iMol)], LevelEeV(iLevels,iMol), Leveljqn(iLevels,iMol), iMol);
        
        %% Computing Vibrational Period
        [Tau(iLevels,iMol)] = ComputePeriod(LevelEh(iLevels,iMol), EpsE, rIn(iLevels,iMol), rMin(iLevels,iMol), rMax(iLevels,iMol), VMax(iLevels,iMol), Leveljqn(iLevels,iMol), iMol);
        
        %% Computing Resonance Width
        Egam(iLevels,iMol) = 0.0;
        if LevelEh(iLevels,iMol) > 0.0
          [Egam(iLevels,iMol)] = ComputeResonanceWidth(rMin(iLevels,iMol), rMax(iLevels,iMol), LevelEh(iLevels,iMol), Leveljqn(iLevels,iMol), Tau(iLevels,iMol), iMol);
        end
        
        %% Plotting Levels on Top of Diatomic Potential
        if (PlotLevelsFlg)
          figure(1+iMol)
          hold on
          plot3([rIn(iLevels,iMol), rOut(iLevels,iMol)], [jqn, jqn], [LevelEeV(iLevels,iMol), LevelEeV(iLevels,iMol)], '-', 'LineWidth', 2)
        end
        
        %% Printing Result @ Screen
        Strr = strcat('Numerov for %i-th Level (v=%i,j=%i) converged in %i Bisection Iterations; Energy = %.',num2str(floor(abs(log10(DeltaESearch)))-1),'fEh \n');
        fprintf(Strr, iLevels, Levelvqn(iLevels,iMol), Leveljqn(iLevels,iMol), iE, LevelEh(iLevels,iMol) )
        fprintf('rIn  for %i-th Level = %fa0 \n',   iLevels, rIn(iLevels,iMol) )
        fprintf('rOut for %i-th Level = %fa0 \n\n', iLevels, rOut(iLevels,iMol) )
        
        
        %% Uncomment for Testing the Calculation of Levels through the Morse Potential's Analytical Solution                
        % De   = 5.0/27.2113839712790;
        % Beta = 1.2;
        % re   = 1.8;
        % EMorse = -De + sqrt(0.5 * Beta^2 / mu * 4 * De) * (Levelvqn(iLevels,iMol) + 0.5) - 0.5 * Beta^2 / mu * (Levelvqn(iLevels,iMol) + 0.5)^2;
        % ErrorMorse = abs(LevelEh(iLevels,iMol) - EMorse);
        % fprintf('Comparison to Morse Levels: Error = %e Eh \n\n', ErrorMorse )
        %% 
        
        vqn = vqn + 1;
      end
        
      
      %% Restarting ~ from here
      ECurrent = EProposed + EpsE*1.e2;
    end % Level
    
  end % Jqn
  
end % iMol