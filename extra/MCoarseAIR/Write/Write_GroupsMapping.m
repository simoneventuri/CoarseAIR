%% Write the Mapping Level-to-Group
%
function Write_GroupsMapping(Controls, iMol)    
    
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

    fprintf('    = Write_GroupsMapping ==============================\n')
    fprintf('    ====================================================\n')
    
    
    
    WriteFldr = strcat(char(Controls.WriteFldr));
    [status,msg,msgID] = mkdir(WriteFldr);
    WriteFldr = strcat(WriteFldr, '/', Syst.NameLong, '/');
    [status,msg,msgID] = mkdir(WriteFldr);
    WriteFldr = strcat(WriteFldr, '/', Syst.Molecule(iMol).Name, '/');
    [status,msg,msgID] = mkdir(WriteFldr);
        
    LevelToGroup = Controls.LevelToGroup;
    if strcmp(char(Controls.Strategy), 'CGM')
        Strategy = char(Input.Kin.CGM_Strategy(iMol));
    else
        Strategy = char(Controls.Strategy);
    end
    NbGroups     = max(LevelToGroup);
    fprintf(['    Grouping Strategy: ' Strategy '\n'] );
    fprintf(['    Nb Groups:         ' num2str(NbGroups) '\n'] );
    fprintf(['    Writing in Folder: ' WriteFldr '\n'] );


    fprintf('    Writing Mapping\n')
    %FileName1 = strcat(WriteFldr, '/LevelsMap_', Strategy ,'.csv');
    FileName1 = strcat(WriteFldr, '/LevelsMap_', Strategy , num2str(NbGroups), '.csv');
    %FileName1 = strcat(WriteFldr, '/LevelsMap_DP20.csv');
    fileID1   = fopen(FileName1,'w');
    fprintf(fileID1,'#Idx,Group\n');

    for iLevel = 1:Syst.Molecule(iMol).NLevels
        fprintf(fileID1,'%i,%i\n', iLevel, LevelToGroup(iLevel));            
    end
    fclose(fileID1);



    fprintf('    Writing Level Properties\n')
    %FileName2 = strcat(WriteFldr, '/LevelsInfo_', Strategy ,'.csv');
    FileName2 = strcat(WriteFldr, '/LevelsInfo_', Strategy , num2str(NbGroups), '.csv');
    %FileName2 = strcat(WriteFldr, '/LevelsInfo_DP20.csv');
    fileID2   = fopen(FileName2,'w');
    fprintf(fileID2,'#Idx,EeV,g,rIn,v,J,ECB,Group\n');

    for iLevel = 1:Syst.Molecule(iMol).NLevels
        fprintf(fileID2,'%i,%e,%e,%e,%i,%i,%e,%i\n',   iLevel,                                      ...
                                                       Syst.Molecule(iMol).LevelEeV(iLevel),        ...
                                                       Syst.Molecule(iMol).Levelg(iLevel),          ...
                                                       Syst.Molecule(iMol).LevelrIn(iLevel),        ...
                                                       Syst.Molecule(iMol).Levelvqn(iLevel),        ...
                                                       Syst.Molecule(iMol).Leveljqn(iLevel),        ...
                                                       Syst.Molecule(iMol).LevelECB(iLevel),        ...
                                                       LevelToGroup(iLevel));            
    end
    fclose(fileID2);

    
    fprintf('    ====================================================\n\n')
     
end