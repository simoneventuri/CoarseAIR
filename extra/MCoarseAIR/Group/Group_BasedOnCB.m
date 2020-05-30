%% Group the Levels Based on Centrifugal Barrier
%
function [LevelToGroup] = Group_BasedOnCB(Controls, iMol)    
    
    %%==============================================================================================================
    % 
    % Coarse-Grained method for Quasi-Classical Trajectories (CG-QCT) 
    % 
    % Copyright (C) 2018 Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
    %
    % Based on "VVTC" (Vectorized Variable stepsize Trajectory Code) by David Schwenke (NASA Ames Research Center). 
    % 
    % This program is free software; you can redistribute it and/or modify it under the terms of the 
    % Version 2.1 GNU Lesser General Public License as published by the Free Software Foundation. 
    % 
    % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
    % without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
    % See the GNU Lesser General Public License for more details. 
    % 
    % You should have received a copy of the GNU Lesser General Public License along with this library; 
    % if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
    % 
    %---------------------------------------------------------------------------------------------------------------
    %%==============================================================================================================

    global Input Kin Param Syst Temp Rates

    fprintf('  = Group_BasedOnCB ==================================\n')
    fprintf('  ====================================================\n')
  
    NBins = Controls.NBins(iMol);
    fprintf('  Nb Bins         = %i \n', NBins );

    alpha = Controls.alpha(iMol);
    fprintf('  Alpha Parameter = %e \n', alpha );


    Extr(1) = abs(max(Syst.Molecule(iMol).LevelECB));
    for i=1:NBins
       Extr(i+1) =  (1.0 - (i/NBins)^alpha) * Extr(1);
    end
    
    LevelToGroup = zeros(Syst.Molecule(iMol).NLevels,1);
    for iLevels = 1:Syst.Molecule(iMol).NLevels
      iBin = 1;
      while (-Syst.Molecule(iMol).LevelECB(iLevels) >= -Extr(iBin)+1e-20)
        iBin = iBin + 1;
      end
      iBin = iBin - 1;
      LevelToGroup(iLevels) = iBin;
    end    
    
    fprintf('  ====================================================\n\n')
    
end