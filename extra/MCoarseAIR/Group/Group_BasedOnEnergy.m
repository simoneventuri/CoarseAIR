%% Group the Levels Based on RoVibrational Energy
%
function [LevelToGroup] = Group_BasedOnEnergy(Syst, Controls, iMol)    
    
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


    fprintf('    = Group_BasedOnEnergy ==============================\n')
    fprintf('    ====================================================\n')
  
    NGroups         = Controls.NGroups_E(iMol);
    NGroups_Bound   = Controls.NGroups_Bound(iMol);
    NGroups_QBound  = NGroups - NGroups_Bound; 
    fprintf('    Nb Groups = %i, of which Bound = %i \n', NGroups, NGroups_Bound );
    
    DissEn          = min(Syst.Molecule(iMol).LevelVMax);
    fprintf('    Dissociation Energy            = %e eV\n', DissEn );
    
    DeltaEeV_Bound  = DissEn - min(Syst.Molecule(iMol).LevelEeV);
    hEeV_Bound      = DeltaEeV_Bound/NGroups_Bound;
    EeVVec_Bound    = [[min(Syst.Molecule(iMol).LevelEeV):hEeV_Bound:DissEn], 1e10];
    fprintf('    Bound Energy Range = %e eV; Bound Bins Energy Range = %e\n', DeltaEeV_Bound, hEeV_Bound)
    
    DeltaEeV_QBound = max(Syst.Molecule(iMol).LevelEeV) - DissEn;
    hEeV_QBound     = DeltaEeV_QBound/NGroups_QBound;
    EeVVec_QBound   = [[DissEn:hEeV_QBound:max(Syst.Molecule(iMol).LevelEeV)], 1e10];
    fprintf('    Q. Bound Energy Range = %e eV; Q. Bound Bins Energy Range = %e\n', DeltaEeV_QBound, hEeV_QBound)

    LevelToGroup = zeros(Syst.Molecule(iMol).NLevels,1);
    for iLevels = 1:Syst.Molecule(iMol).NLevels    
        if (Syst.Molecule(iMol).LevelEeV(iLevels) <= DissEn)
            iGroup = 1;
            while (Syst.Molecule(iMol).LevelEeV(iLevels) > EeVVec_Bound(iGroup+1))
                iGroup = iGroup + 1;
            end
        else
            iGroup = + 1;
            while (Syst.Molecule(iMol).LevelEeV(iLevels) > EeVVec_QBound(iGroup+1))
                iGroup = iGroup + 1;
            end
            iGroup = iGroup +  NGroups_Bound;
        end
        LevelToGroup(iLevels) = iGroup;
    end    
    
    fprintf('    ====================================================\n\n')
    
end