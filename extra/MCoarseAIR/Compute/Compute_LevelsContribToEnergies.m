%% Computing the Energy Stored in each of the Internal Modes
%
function Compute_LevelsContribToEnergies(Controls)    

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
    
    global Input Kin Param Syst Temp

    
    fprintf('= Compute_LevelsContribToEnergies ====== T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    

    for iMol = Controls.MoleculesOI
        fprintf(['Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );

        clear LevelToBin Levelvqn LevelEeV LevelEeV0 LevelEeVRot Nvqn NLevels LevelPop PotTot eInt eRot eVib vPop
        LevelToBin  = Syst.Molecule(iMol).LevelToGroupIn;
        Levelvqn    = Syst.Molecule(iMol).Levelvqn;
        LevelEeV    = Syst.Molecule(iMol).LevelEeV;
        LevelEeV0   = Syst.Molecule(iMol).LevelEeV0;
        LevelEeVRot = Syst.Molecule(iMol).LevelEeVRot;
        LevelECB    = Syst.Molecule(iMol).LevelECB;
        vEeVVib0    = Syst.Molecule(iMol).vEeVVib0;
        Nvqn        = Syst.Molecule(iMol).Nvqn;
        NLevels     = Syst.Molecule(iMol).NLevels;
        
        clear NGroups hEeV GroupsEeV LevelToGroup
%         NGroups      = Controls.NGroups;
%         hEeV         = (max(LevelEeV0) - min(LevelEeV0))/NGroups;
%         GroupsEeV    = linspace(min(LevelEeV0), max(LevelEeV0), NGroups+1);
%         LevelToGroup = zeros(NLevels,1);
%         for iLevel=1:NLevels
%            iGroup=2;
%            while (LevelEeV0(iLevel)> GroupsEeV(iGroup))
%                iGroup = iGroup+1;
%            end
%            LevelToGroup(iLevel) = iGroup-1;
%         end
        
        NGroups        = Controls.NGroups;
        hEeV           = (max(LevelECB) - min(LevelECB))/NGroups;
        GroupsEeV      = linspace(-max(LevelECB), -min(LevelECB), NGroups+1);
        LevelToGroup   = zeros(NLevels,1);
        for iLevel=1:NLevels
           iGroup=2;
           while (-LevelECB(iLevel)> GroupsEeV(iGroup))
               iGroup = iGroup+1;
           end
           LevelToGroup(iLevel) = iGroup-1;
        end
        
        
        
        eInt = zeros(Kin.T(Temp.iT).NSteps,NGroups);
        %eRot = zeros(Kin.T(Temp.iT).NSteps,1);
        %eVib = zeros(Kin.T(Temp.iT).NSteps,1);        
        %vPop = zeros(Nvqn,1);
        for iStep = 1:Kin.T(Temp.iT).NSteps

            
            if strcmp(Syst.Molecule(iMol).KinMthdIn, 'StS')
                LevelPop(:,1) = Kin.T(Temp.iT).Molecule(iMol).Pop(iStep,:);            
            else
                LevelPop(:,1) = Kin.T(Temp.iT).Molecule(iMol).PopOverg(iStep,LevelToBin(:))' .* Syst.Molecule(iMol).T(Temp.iT).Levelq(:);
            end
            
            %vPop = vPop.*0.0;
            %for iLevels = 1:NLevels
            %    vPop(Levelvqn(iLevels)+1) = vPop(Levelvqn(iLevels)+1) + LevelPop(iLevels);
            %end
            
            PopTot = zeros(NGroups,1);
            for iLevel=1:NLevels
                iGroup              = LevelToGroup(iLevel);
                eInt(iStep, iGroup) = eInt(iStep, iGroup) + (LevelEeV0(iLevel) .* LevelPop(iLevel) );
                PopTot(iGroup)      = PopTot(iGroup)      + LevelPop(iLevel);
            end

            for iGroup=2:NGroups
                eInt(iStep, iGroup) = eInt(iStep, iGroup) + eInt(iStep, iGroup-1);
                PopTot(iGroup)      = PopTot(iGroup)      + PopTot(iGroup-1);
            end
            eInt(iStep, :) = eInt(iStep, :) ./ PopTot(:,1)';
       
            %eRot(iStep) = sum( LevelEeVRot' .* LevelPop ) ./ PopTot;
            %eVib(iStep) = sum( vEeVVib0'    .* vPop )     ./ PopTot;
        end
        
        
        figure(Input.iFig)
        fig = gcf;
        screensize   = get( groot, 'Screensize' );
        %fig.Position = screensize;
        %fig.Color='None';
        for iGroup = 1:NGroups

            h1=semilogx(Kin.T(Temp.iT).t, eInt(:, iGroup)./eInt(:, NGroups), '-', 'Color', Param.KCVec, 'LineWidth', Param.LineWidth);
            hold on

            xt = get(gca, 'XTick');
            set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
            yt = get(gca, 'YTick');
            set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

            str_x = ['t [s]'];
            xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
            xlab.Interpreter = 'latex';
            %xlim([max(min(LevelEeV)), MinEvPlot, min(max(LevelEeV)), MaxEvPlot]);

            str_y = ['Internal Energy [eV]'];
            ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
            ylab.Interpreter = 'latex';
            %ylim([1.d5, 1.d23]);
            %set(gca, 'YScale', 'log')


            pbaspect([1 1 1])

            if Input.SaveFigsFlgInt > 0
                [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
                FolderPath = strcat(Input.Paths.SaveFigsFldr, '/T_', Temp.TNowChar, 'K_', Input.Kin.Proc.OverallFlg, '/');
                [status,msg,msgID] = mkdir(FolderPath);
                FileName = strcat(Syst.Molecule(iMol).Name,'_LevelsContribToEnergies');
                if Input.SaveFigsFlgInt == 1
                    FileName   = strcat(FolderPath, FileName);
                    export_fig(FileName, '-pdf')
                elseif Input.SaveFigsFlgInt == 2
                    FileName   = strcat(FolderPath, strcat(FileName,'.fig'));
                    savefig(FileName)
                end
                %close
            end
            Input.iFig = Input.iFig + 1;
            
        end
        
       
    end

    fprintf('====================================================\n\n')        

    
end