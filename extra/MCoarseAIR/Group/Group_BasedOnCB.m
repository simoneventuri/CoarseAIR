%% Group the Levels Based on Centrifugal Barrier
%
function [LevelToGroup] = Group_BasedOnCB(Syst, Controls, iMol)    
    
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


    fprintf('    = Group_BasedOnCB ==================================\n')
    fprintf('    ====================================================\n')
  
    NGroups = Controls.NGroups_CB(iMol);
    fprintf('    Nb Groups       = %i \n', NGroups );

    alpha = Controls.alpha(iMol);
    fprintf('    Alpha Parameter = %e \n', alpha );
    
    MinEeV = min(Controls.MinEeV(iMol), abs(max(Syst.Molecule(iMol).LevelECB)));
    fprintf('    Maximum Distance from Centrifugal Barrier = %e eV\n', MinEeV)

    Extr(1) = MinEeV;
    for i=1:NGroups
       Extr(i+1) = (1.0 - i/(NGroups))^(alpha) * MinEeV;
    end
%     figure
%     plot(Extr,'o-')
    
    LevelToGroup = zeros(Syst.Molecule(iMol).NLevels,1);
    for iLevels = 1:Syst.Molecule(iMol).NLevels
        if (Syst.Molecule(iMol).LevelECB(iLevels) <= MinEeV)
            iGroup = 1;
            while (Syst.Molecule(iMol).LevelECB(iLevels) < Extr(iGroup+1))
                iGroup = iGroup + 1;
            end
            LevelToGroup(iLevels) = iGroup;
        end
    end    
    
    fprintf('    ====================================================\n\n')
%     
%     min(LevelToGroup)
%     max(LevelToGroup)
end