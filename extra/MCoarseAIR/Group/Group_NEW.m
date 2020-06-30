%% The Function writes the Rates in the Format for being read by Amal's Clustering Algorithm
%
function Group_NEW(Controls)

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

    
    fprintf('= Group_NEW ============================ T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    
    NbGroups      = 20; 
    NbGroups_Inel = 10;
    NbGroups_CB   = NbGroups - NbGroups_Inel;
    
    


    WriteFldr = strcat(Controls.WriteFldr);
    [status,msg,msgID] = mkdir(WriteFldr);
    WriteFldr = strcat(WriteFldr, '/', Syst.NameLong, '/');
    [status,msg,msgID] = mkdir(WriteFldr);
    WriteFldr = strcat(WriteFldr, '/T', Temp.TNowChar, 'K/');
    [status,msg,msgID] = mkdir(WriteFldr);

    
    for iMol = 1:1
        fprintf(['Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );
        
        
        Nvqn     = Syst.Molecule(iMol).Nvqn; 
        NLevels  = Syst.Molecule(iMol).NLevels; 
        Levelvqn = Syst.Molecule(iMol).Levelvqn;
        Leveljqn = Syst.Molecule(iMol).Leveljqn;
        LevelEeV = Syst.Molecule(iMol).LevelEeV;
        
        MinState = 1; %max(Controls.MinState, 1);
        MaxState = NLevels; %min(Controls.MaxState, Syst.Molecule(iMol).NLevels);
        

        Controls.MinEeV(iMol) = abs(Syst.Molecule(iMol).DissEn) / 2.0;
        Controls.NBins        = NbGroups_CB;
        Controls.alpha(iMol)  = 1.0/2.0; 
        LevelToGroupCB        = Group_BasedOnCB(Syst, Controls, iMol);
    
        
        vPlus_1 = 3;
        vPlus_2 = 10;
        
        EeV_J0 = [];
        MaxVqn = 0;
        for iLevel = 1:NLevels
            if (Leveljqn(iLevel) == 0) 
                EeV_J0(Levelvqn(iLevel)+1) = LevelEeV(iLevel);
                if (LevelToGroupCB(iLevel) == 0)
                    MaxVqn = max(MaxVqn,Levelvqn(iLevel)+1);
                end
            end
        end
        
        
        MaxSplit1 = 0;
        MaxSplit2 = 0;
        MaxSplit3 = 0;
        for iv=0:MaxVqn
            for iLevel = 1:NLevels
                if (Levelvqn(iLevel) == iv) && (LevelToGroupCB(iLevel) == 0)
                   if LevelEeV(iLevel)     < EeV_J0(iv+1+vPlus_1)
                       Split(iLevel,:) = [iv+1, 1];
                       MaxSplit1 = max(MaxSplit1,iv+1);
                   elseif LevelEeV(iLevel) < EeV_J0(iv+1+vPlus_2)
                       Split(iLevel,:) = [iv+1, 2];
                       MaxSplit2 = max(MaxSplit2,iv+1);
                   else
                       Split(iLevel,:) = [iv+1, 3];
                       MaxSplit3 = max(MaxSplit3,iv+1);
                   end
                end
            end
        end
        
%         PreGroups(1) = 0;
%         PreGroups(2) = MaxSplit1;
%         PreGroups(3) = MaxSplit1 + MaxSplit2;
%         PreGroups(4) = MaxSplit1 + MaxSplit2 + MaxSplit3;
%         LevelToGroup = zeros(NLevels,1);
%         for iLevel = 1:NLevels
%             if (LevelToGroupCB(iLevel) == 0)
%                 LevelToGroup(iLevel) = Split(iLevel,1) + PreGroups(Split(iLevel,2));
%             else
%                 LevelToGroup(iLevel) = PreGroups(4)    + LevelToGroupCB(iLevel);
%             end
%         end
        
        
        
        NbGroups_Inel_1 = ceil(NbGroups_Inel / 2    );
        NbGroups_Inel_2 = ceil(NbGroups_Inel * 2 / 3);
        NbGroups_Inel_3 = NbGroups_Inel - (NbGroups_Inel_1+NbGroups_Inel_2);
        
        MaxVqn_1        = MaxSplit1;
        NbVqnPerGroup_1 = ceil( ((MaxVqn_1) + 1 - NbGroups_Inel_1) / (NbGroups_Inel_1 - 3) );
        VqnExtr_1       = [0, 1, 2, 3];
        iAdd            = 0;
        while (VqnExtr_1(end)+NbVqnPerGroup_1+iAdd < MaxVqn_1)
            VqnExtr_1 = [VqnExtr_1, (VqnExtr_1(end)+NbVqnPerGroup_1+iAdd)];
            iAdd      = iAdd+1;
        end    
        %VqnExtr_1

        
        MaxVqn_2        = MaxSplit2;
        NbVqnPerGroup_2 = ceil( ((MaxVqn_2) + 1 - NbGroups_Inel_2) / (NbGroups_Inel_2 - 1) );
        VqnExtr_2       = [0, 1];
        iAdd            = 0;
        while (VqnExtr_2(end)+NbVqnPerGroup_2+iAdd < MaxVqn_2)
            VqnExtr_2 = [VqnExtr_2, (VqnExtr_2(end)+NbVqnPerGroup_2+iAdd)];
            iAdd      = iAdd+1;
        end    
        %VqnExtr_2
        

        MaxVqn_3        = MaxSplit3;
        NbVqnPerGroup_3 = ceil( ((MaxVqn_3) + 1 - NbGroups_Inel_3) / (NbGroups_Inel_3) );
        VqnExtr_3       = [0];
        iAdd            = 0;
        while (VqnExtr_3(end)+NbVqnPerGroup_3+iAdd < MaxVqn_3)
            VqnExtr_3 = [VqnExtr_3, (VqnExtr_3(end)+NbVqnPerGroup_3+iAdd)];
            iAdd    = iAdd+1;
        end    
        %VqnExtr_3
        
        


        PreGroups(1) = 0;
        PreGroups(2) = NbGroups_Inel_1;
        PreGroups(3) = NbGroups_Inel_1 + NbGroups_Inel_2;
        PreGroups(4) = NbGroups_Inel_1 + NbGroups_Inel_2 + NbGroups_Inel_3;
        LevelToGroup = zeros(NLevels,1);
        for iLevel = 1:NLevels
            if (LevelToGroupCB(iLevel) == 0)
                if Split(iLevel,2)     == 1
                    iGroup = sum(Split(iLevel,1) > VqnExtr_1);
                elseif Split(iLevel,2) == 2
                    iGroup = sum(Split(iLevel,1) > VqnExtr_2);
                elseif Split(iLevel,2) == 3
                    iGroup = sum(Split(iLevel,1) > VqnExtr_3);
                end
                LevelToGroup(iLevel) = iGroup + PreGroups(Split(iLevel,2));
            else
                LevelToGroup(iLevel) = PreGroups(4)    + LevelToGroupCB(iLevel);
            end
        end
        
        
        
        fprintf('Writing Level Properties\n')
        
        FileName2 = strcat(WriteFldr, '/LevelsInfo_NEW.csv');
        fileID2   = fopen(FileName2,'w');
        fprintf(fileID2,'#Idx,EeV,g,rIn,v,J,ECB,Group\n');
         
        for iLevel = MinState:MaxState
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


     end
     
     
end