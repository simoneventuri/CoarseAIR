%% The Function writes the Rates in the Format for being read by Amal's Clustering Algorithm
%
function [LevelToGroup] = Group_BasedOnDP(Syst, Controls, iMol)

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

    global Temp

    
    fprintf('    = Group_BasedOnDP ==================================\n')
    fprintf('    ====================================================\n')


    NGroups       = Controls.NGroups(iMol); 
    NGroups_Excit = Controls.NGroups_Excit(iMol);
    NGroups_CB    = NGroups - NGroups_Excit;
    fprintf(['    Nb Groups: ' num2str(NGroups) ', of which: ' num2str(NGroups_Excit) ' for Excitation and ' num2str(NGroups_CB) ' for Dissociation \n\n'] );

    Nvqn     = Syst.Molecule(iMol).Nvqn; 
    NLevels  = Syst.Molecule(iMol).NLevels; 
    Levelvqn = Syst.Molecule(iMol).Levelvqn;
    Leveljqn = Syst.Molecule(iMol).Leveljqn;
    LevelEeV = Syst.Molecule(iMol).LevelEeV;

    MinState = 1; %max(Controls.MinState, 1);
    MaxState = NLevels; %min(Controls.MaxState, Syst.Molecule(iMol).NLevels);


    Controls.MinEeV(iMol)     = abs(Syst.Molecule(iMol).DissEn) / 2.0;
    Controls.NGroups_CB(iMol) = NGroups_CB;
    Controls.alpha(iMol)      = 1.0/2.0; 
    LevelToGroupCB            = Group_BasedOnCB(Syst, Controls, iMol);


    vPlus_1 = 3;
    vPlus_2 = 12;

    EeV_J0 = ones(200,1) .* 100.0;
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



    NGroups_Excit_1 = floor(NGroups_Excit / 2    );
    NGroups_Excit_2 = floor(NGroups_Excit * 2 / 5);
    NGroups_Excit_3 = NGroups_Excit - (NGroups_Excit_1+NGroups_Excit_2);

    
    VqnExtr_1_ = ones(NGroups_Excit_1-2,1);
    iTry       = 0;
    TotVqn     = 0;
    while (TotVqn < MaxSplit1-NGroups_Excit_1)
        iTry   = iTry + 1;
        TotVqn = iTry^2;
    end
    iTry = iTry-1;
    VqnExtr_1_ = [1; 1; VqnExtr_1_];
    ii=0;
    for jTry=iTry:-1:1
        VqnExtr_1_(end-ii) = VqnExtr_1_(end-ii) + (2*jTry-1); 
        ii=ii+1;
    end
    VqnExtr_1_(end)   = VqnExtr_1_(end)   + floor(2/3*(MaxSplit1 - sum(VqnExtr_1_))); 
    VqnExtr_1_(end-1) = VqnExtr_1_(end-1) +           (MaxSplit1 - sum(VqnExtr_1_)); 
    VqnExtr_1  = [0];
    for iGroup = 1:NGroups_Excit_1
        VqnExtr_1 = [VqnExtr_1, VqnExtr_1(iGroup)+VqnExtr_1_(iGroup)];
    end
    
    
    
    VqnExtr_2_ = ones(NGroups_Excit_2-1,1);
    iTry       = 0;
    TotVqn     = 0;
    while (TotVqn < MaxSplit2-NGroups_Excit_2)
        iTry   = iTry + 1;
        TotVqn = iTry^2;
    end
    iTry = iTry-1;
    VqnExtr_2_ = [1; VqnExtr_2_];
    ii=0;
    for jTry=iTry:-1:1
        VqnExtr_2_(end-ii) = VqnExtr_2_(end-ii) + (2*jTry-1); 
        ii=ii+1;
    end
    VqnExtr_2_(end)   = VqnExtr_2_(end)   + floor(2/3*(MaxSplit2 - sum(VqnExtr_2_))); 
    VqnExtr_2_(end-1) = VqnExtr_2_(end-1) +           (MaxSplit2 - sum(VqnExtr_2_)); 
    VqnExtr_2  = [0];
    for iGroup = 1:NGroups_Excit_2
        VqnExtr_2 = [VqnExtr_2, VqnExtr_2(iGroup)+VqnExtr_2_(iGroup)];
    end
    
    
    
    if (NGroups_Excit_3 > 1)    
        VqnExtr_3_ = ones(NGroups_Excit_3,1);
        iTry       = 0;
        TotVqn     = 0;
        while (TotVqn < MaxSplit3-NGroups_Excit_3)
            iTry   = iTry + 1;
            TotVqn = iTry^2;
        end
        iTry = iTry-1;
        ii=0;
        for jTry=iTry:-1:1
            VqnExtr_3_(end-ii) = VqnExtr_3_(end-ii) + (2*jTry-1); 
            ii=ii+1;
        end
        VqnExtr_3_(end)   = VqnExtr_3_(end)   + floor(2/3*(MaxSplit3 - sum(VqnExtr_3_))); 
        VqnExtr_3_(end-1) = VqnExtr_3_(end-1) +           (MaxSplit3 - sum(VqnExtr_3_)); 
        VqnExtr_3  = [0];
        for iGroup = 1:NGroups_Excit_3
            VqnExtr_3 = [VqnExtr_3, VqnExtr_3(iGroup)+VqnExtr_3_(iGroup)];
        end
    else
       VqnExtr_3 = [0, MaxSplit3] 
    end
 
    
    
    PreGroups(1) = 0;
    PreGroups(2) = NGroups_Excit_1;
    PreGroups(3) = NGroups_Excit_1 + NGroups_Excit_2;
    PreGroups(4) = NGroups_Excit_1 + NGroups_Excit_2 + NGroups_Excit_3;
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
    
    
    fprintf('    ====================================================\n\n')
          
end