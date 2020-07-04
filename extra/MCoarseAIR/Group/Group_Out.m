%% Group the Levels in Output
%
function Group_Out()    
    
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

    fprintf('= Group_Out ========================================\n')
    fprintf('====================================================\n')
    
    
    for iMol = 1:Syst.NMolecules
        fprintf(['Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );
        clear LevelToGroup
        
        if (strcmp(Input.Kin.MolResolutionOut(iMol), 'StS'))
            fprintf('Not Grouping the Molecule\n')

            LevelToGroup = [1:Syst.Molecule(iMol).NLevels]';
            
        elseif (strcmp(Input.Kin.MolResolutionOut(iMol), 'VSM'))
            fprintf('Grouping based on Vibrational Quantum Number\n\n')
            
            LevelToGroup = Group_BasedOnVib(Syst, iMol);
            
        elseif (strcmp(Input.Kin.MolResolutionOut(iMol), 'CGM'))
            fprintf('Grouping for Coarse-Grained Model\n')
            
            if (strcmp(Input.Kin.CGM_Strategy(iMol), 'RVE'))
                fprintf('  Gropuing based on RoVibrational Energy\n\n')

                Controls.NGroups_E(iMol) = Input.Kin.NGroupsOut(iMol);
                Controls.DissEn(iMol)        = 0.0;
                Controls.NGroups_Bound(iMol) = Input.Kin.ParamsGroupsOut(iMol);
                LevelToGroup = Group_BasedOnEnergy(Syst, Controls, iMol);  
            
            elseif (strcmp(Input.Kin.CGM_Strategy(iMol), 'CBM'))
                fprintf('Gropuing based on Energy-Distance from Centrifugal Barrier\n\n')

                Controls.NGroups_CB(iMol) = Input.Kin.NGroupsOut(iMol);
                Controls.alpha(iMol)      = Input.Kin.ParamsGroupsOut(iMol);
                LevelToGroup = Group_BasedOnCB(Syst, Controls, iMol);  
                
            elseif (strcmp(Input.Kin.CGM_Strategy(iMol), 'DPM'))
                fprintf('Gropuing based on Diatomic Potential\n\n')

                Controls.NGroups(iMol)       = Input.Kin.NGroupsOut(iMol);
                Controls.NGroups_Excit(iMol) = ceil(Input.Kin.ParamsGroupsOut(iMol) * Controls.NGroups(iMol));
                LevelToGroup = Group_BasedOnDP(Syst, Controls, iMol);  
                
            elseif (strcmp(Input.Kin.CGM_Strategy(iMol), 'File'))
                fprintf('Reading Levels-Groups Mapping from File\n\n')
                
                Controls.FilePath(iMol) = Input.Kin.PathToMappingOut(iMol);
                LevelToGroup = Group_FromFile(Syst, Controls, iMol);
                
            end
        
        end
        
        Syst.Molecule(iMol).LevelToGroupOut = LevelToGroup;
        if min(Syst.Molecule(iMol).LevelToGroupOut) < 1
            fprintf('ERROR! Some Levels are not assigned to Groups! \n\n')
            break;
        end
        Syst.Molecule(iMol).NGroupsOut      = max(Syst.Molecule(iMol).LevelToGroupOut);
        
        
        if not(strcmp(Input.Kin.PathToWriteMappingOut, '')) && not(strcmp(char(Input.Kin.MolResolutionOut(iMol)), 'StS'))
            
            clear Controls.Strategy
            Controls.WriteFldr    = Input.Kin.PathToWriteMappingOut;
            Controls.LevelToGroup = Syst.Molecule(iMol).LevelToGroupOut;
            Controls.Strategy     = Input.Kin.MolResolutionOut(iMol);
            Write_GroupsMapping(Controls, iMol)
            
        end
        
    end
    
    
    fprintf('====================================================\n\n')
    
end