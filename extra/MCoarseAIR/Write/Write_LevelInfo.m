%% The Function Writes the Molecules' Level Info from the list used for QCT
%        
function Write_LevelInfo(Flg)
      
    % -- MATLAB --
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
        
    global Syst
    
    fprintf('= Write_LevelInfo ==================================\n')
    fprintf('====================================================\n')
    
    for iMol = 1:Syst.NMolecules       
        fprintf('Writing Level Quantities for Molecule Nb %i \n',  iMol )

        fileID  = fopen(char(Flg.Path(iMol)),'w');
        fprintf(fileID,'######################################################################################################################################################\n');
        fprintf(fileID,'# jqn   : the rotational q.n. of the i-th quantum state                                                                                               \n');
        fprintf(fileID,'# vqn   : the vibrational q.n. of the i-th quantum state                                                                                              \n');
        fprintf(fileID,'# eint  : internal energy of i-th quantum state                                                                                                       \n');
        fprintf(fileID,'# egam  : Half width of i-th quantum state                                                                                                            \n');
        fprintf(fileID,'# rmin  : the position of the potential minimum (included centrifugal potential) for i-th quantum state                                               \n');
        fprintf(fileID,'# vmin  : the value of the potential minimun (inc. cent. pot.)                                                                                        \n');
        fprintf(fileID,'# vmax  : the value of the local potential maximum (inc. cent. pot.)                                                                                  \n');
        fprintf(fileID,'# tau   : the vibrational period of the i-th quantum state                                                                                            \n');
        fprintf(fileID,'# ri    : inner turning point                                                                                                                         \n');
        fprintf(fileID,'# ro    : outter turning point                                                                                                                        \n');
        fprintf(fileID,'# rmax  : location of maximum in centrifugal barrier                                                                                                  \n');
        fprintf(fileID,'######################################################################################################################################################\n');
        fprintf(fileID,'#   vqn  jqn      eint           egam           rmin           rmax           vmin           vmax            tau             ri             ro        \n');
        fprintf(fileID,'######################################################################################################################################################\n');
        for iLevels = 1:Syst.Molecule(iMol).NLevels
            fprintf(fileID,'%6d %4d %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E\n',  Syst.Molecule(iMol).Levelvqn(iLevels),  ...
                                                                                                        Syst.Molecule(iMol).Leveljqn(iLevels),  ...
                                                                                                        Syst.Molecule(iMol).LevelEEh(iLevels),  ...
                                                                                                        Syst.Molecule(iMol).LevelEgam(iLevels), ...
                                                                                                        Syst.Molecule(iMol).LevelrMin(iLevels), ...
                                                                                                        Syst.Molecule(iMol).LevelrMax(iLevels), ...
                                                                                                        Syst.Molecule(iMol).LevelVMin(iLevels), ...
                                                                                                        Syst.Molecule(iMol).LevelVMax(iLevels), ...
                                                                                                        Syst.Molecule(iMol).LevelTau(iLevels),  ...
                                                                                                        Syst.Molecule(iMol).LevelrIn(iLevels),  ...
                                                                                                        Syst.Molecule(iMol).LevelrOut(iLevels) );
        end
        fclose(fileID);
        
    end
    
    fprintf('====================================================\n')
    
end