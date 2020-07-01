%% Group the Levels in Input
%
function [Syst] = Group_In(Syst)    
    
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

    global Input Kin Param  Temp Rates

    fprintf('  = Group_In =========================================\n')
    fprintf('  ====================================================\n')
    
    
    for iMol = 1:Syst.NMolecules
        fprintf(['  Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );        
        
        if (strcmp(Input.Kin.MolResolutionIn(iMol), 'StS'))
            fprintf('  Not Grouping the Molecule\n')

            LevelToGroup = [1:Syst.Molecule(iMol).NLevels]';
            
        elseif (strcmp(Input.Kin.MolResolutionIn(iMol), 'VSM'))
            fprintf('  Grouping based on Vibrational Quantum Number\n\n')
            
            LevelToGroup = Group_BasedOnVib(Syst, iMol);
            
        elseif (strcmp(Input.Kin.MolResolutionIn(iMol), 'CGM'))
            fprintf('  Grouping for Coarse-Grained Model\n')
                            
            if (strcmp(Input.Kin.CGM_Strategy(iMol), 'CBM'))
                fprintf('  Gropuing based on Energy-Distance from Centrifugal Barrier\n\n')

                Controls.NGroups_CB(iMol) = Input.Kin.NGroupsIn(iMol);
                Controls.alpha(iMol)      = Input.Kin.ParamsGroupsIn(iMol);
                LevelToGroup = Group_BasedOnCB(Syst, Controls, iMol);  
                
            elseif (strcmp(Input.Kin.CGM_Strategy(iMol), 'DPM'))
                fprintf('Gropuing based on Diatomic Potential\n\n')
                
                Controls.NGroups(iMol)       = Input.Kin.NGroupsIn(iMol);
                Controls.NGroups_Excit(iMol) = ceil(Input.Kin.ParamsGroupsIn(iMol) * Controls.NGroups(iMol));
                LevelToGroup = Group_BasedOnDP(Syst, Controls, iMol);  
                
            elseif (strcmp(Input.Kin.CGM_Strategy(iMol), 'File'))
                fprintf('  Reading Levels-Groups Mapping from File\n\n')
                
                Controls.FilePath(iMol) = Input.Kin.PathToMappingIn(iMol);
                LevelToGroup = Group_FromFile(Syst, Controls, iMol);
                
            end
        
        end
        
        Syst.Molecule(iMol).LevelToGroupIn = LevelToGroup;
        if min(Syst.Molecule(iMol).LevelToGroupIn) < 1
            fprintf('  ERROR! Some Levels are not assigned to Groups! \n\n')
            break;
        end
        Syst.Molecule(iMol).NGroupsIn      = max(Syst.Molecule(iMol).LevelToGroupIn);
        
        
        if not(strcmp(Input.Kin.PathToWriteMappingIn(iMol), ''))
            
            Controls.WriteFldr                   = Input.Kin.PathToWriteMappingIn(iMol);
            Controls.Molecule(iMol).LevelToGroup = Syst.Molecule(iMol).LevelToGroupIn;
            Controls.Molecule(iMol).Strategy     = Input.Kin.MolResolutionIn(iMol);
            Write_GroupsMapping(Controls, iMol)
            
        end
        
    end
    
    
    fprintf('====================================================\n\n')
    
end