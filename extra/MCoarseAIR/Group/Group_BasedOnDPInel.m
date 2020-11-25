%% The Function writes the Rates in the Format for being read by Amal's Clustering Algorithm
%
function [LevelToGroup] = Group_BasedOnDPInel(Syst, Controls, iMol)

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

    DbgFlg = false;
    
    fprintf('    = Group_BasedOnDPInel ==============================\n')
    fprintf('    ====================================================\n')


    NGroups       = Controls.NGroups(iMol); 
    NGroups_Excit = NGroups;
    fprintf(['    Nb Groups: ' num2str(NGroups) ', of which: ' num2str(NGroups_Excit) ' for Excitation \n\n'] );

    Nvqn     = Syst.Molecule(iMol).Nvqn; 
    NLevels  = Syst.Molecule(iMol).NLevels; 
    Levelvqn = Syst.Molecule(iMol).Levelvqn;
    Leveljqn = Syst.Molecule(iMol).Leveljqn;
    LevelEeV = Syst.Molecule(iMol).LevelEeV;

    MinState = 1; %max(Controls.MinState, 1);
    MaxState = NLevels; %min(Controls.MaxState, Syst.Molecule(iMol).NLevels);

    EeV_J0 = ones(200,1) .* 100.0;
    MaxVqn = 0;
    for iLevel = 1:NLevels
        if (Leveljqn(iLevel) == 0) 
            EeV_J0(Levelvqn(iLevel)+1) = LevelEeV(iLevel);
            MaxVqn = max(MaxVqn,Levelvqn(iLevel)+1);
        end
    end

    vPlus_1 = 3;
    vPlus_2 = 8;
    
    fprintf(['    For Excitation: \n'] );

    MaxSplit1 = 0;
    MaxSplit2 = 0;
    MaxSplit3 = 0;
    for iv=0:MaxVqn
        for iLevel = 1:NLevels
            if (Levelvqn(iLevel) == iv)
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


    NGroups_Excit_1 = floor(NGroups_Excit / 2    );
    NGroups_Excit_2 = floor(NGroups_Excit * 2 / 5);
    NGroups_Excit_3 = NGroups_Excit - (NGroups_Excit_1+NGroups_Excit_2);
    fprintf(['    Nb Groups in Region 1: ' num2str(NGroups_Excit_1) '; Nb Groups in Region 2: ' num2str(NGroups_Excit_2)  '; Nb Groups in Region 3: ' num2str(NGroups_Excit_3)  '\n'] );

    
    NBinss    = NGroups_Excit_1;
    maxvqn    = double(MaxSplit1);
    x         = [1:1:NBinss];
    Deltap    = 1.e-3;
    p         = 0.7-Deltap;
    SumValues = 0;
    while SumValues ~= maxvqn
        p         = p + Deltap;
        GeoDist   = 1 + round(p.*(1-p).^(max(x)-x) .* (maxvqn+0.1 - NBinss));
        SumValues = sum(GeoDist);
    end
    fprintf(['      Region 1, p                    = ' num2str(p) '\n'] );
    fprintf(['      Region 1, Nb of vqns per Group = ' num2str(GeoDist) '\n'] );
    if (DbgFlg)
        figure(4321)
        hold on
        plot(x,GeoDist)
    end
    VqnExtr_1  = [0];
    for iGroup = 1:NGroups_Excit_1
        VqnExtr_1 = [VqnExtr_1, GeoDist(iGroup)+VqnExtr_1(iGroup)];
    end
    

    clear NBins maxvqn x Deltap p SumValues GeoDist
    NBins     = NGroups_Excit_2;
    maxvqn    = double(MaxSplit2);
    x         = [1:1:NBins];
    Deltap    = 1.e-3;
    p         = 0.7-Deltap;
    SumValues = 0;
    while SumValues ~= maxvqn
        p         = p + Deltap;
        GeoDist   = 1 + round(p*(1-p).^(max(x)-x) .* (maxvqn - NBins));
        SumValues = sum(GeoDist);
    end
    fprintf(['      Region 2, p                    = ' num2str(p) '\n'] );
    fprintf(['      Region 2, Nb of vqns per Group = ' num2str(GeoDist) '\n'] );
    if (DbgFlg)
        figure(4321)
        hold on
        plot(x,GeoDist)
    end
    VqnExtr_2  = [0];
    for iGroup = 1:NGroups_Excit_2
        VqnExtr_2 = [VqnExtr_2, GeoDist(iGroup)+VqnExtr_2(iGroup)];
    end
    
    
    clear NBins maxvqn x Deltap p SumValues GeoDist
    NBins     = NGroups_Excit_3;
    maxvqn    = double(MaxSplit3);
    x         = [1:1:NBins];
    Deltap    = 1.e-3;
    p         = 0.7-Deltap;
    SumValues = 0;
    while SumValues ~= maxvqn
        p         = p + Deltap;
        GeoDist   = 1 + round(p*(1-p).^(max(x)-x) .* (maxvqn - NBins));
        SumValues = sum(GeoDist);
    end
    fprintf(['      Region 3, p                    = ' num2str(p) '\n'] );
    fprintf(['      Region 3, Nb of vqns per Group = ' num2str(GeoDist) '\n'] );
    if (DbgFlg)
        figure(4321)
        hold on
        plot(x,GeoDist)
        hold off
    end
    VqnExtr_3  = [0];
    for iGroup = 1:NGroups_Excit_3
        VqnExtr_3 = [VqnExtr_3, GeoDist(iGroup)+VqnExtr_3(iGroup)];
    end

    
    PreGroups(1) = 0;
    PreGroups(2) = NGroups_Excit_1;
    PreGroups(3) = NGroups_Excit_1 + NGroups_Excit_2;
    PreGroups(4) = NGroups_Excit_1 + NGroups_Excit_2 + NGroups_Excit_3;
    LevelToGroup = zeros(NLevels,1);
    for iLevel = 1:NLevels
        if Split(iLevel,2)     == 1
            iGroup = sum(Split(iLevel,1) > VqnExtr_1);
        elseif Split(iLevel,2) == 2
            iGroup = sum(Split(iLevel,1) > VqnExtr_2);
        elseif Split(iLevel,2) == 3
            iGroup = sum(Split(iLevel,1) > VqnExtr_3);
        end
        LevelToGroup(iLevel) = iGroup + PreGroups(Split(iLevel,2));
    end
    
    
    fprintf('    ====================================================\n\n')
          
end