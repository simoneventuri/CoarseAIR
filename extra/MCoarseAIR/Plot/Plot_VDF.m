%% The Function plots the Vibrational Distribution Function
%
function Plot_VDF(Controls)    
    
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

    fprintf('= Plot_VDF ============================= T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    
    
    for iMol = Controls.MoleculesOI
        fprintf(['Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );

        LevelToBin = Syst.Molecule(iMol).LevelToGroupIn;
        Levelvqn   = Syst.Molecule(iMol).Levelvqn;
        LevelEeV   = Syst.Molecule(iMol).LevelEeV;
        Levelq     = Syst.Molecule(iMol).T(Temp.iT).Levelq;
        Levelg     = Syst.Molecule(iMol).Levelg;
        
        
        for tStep = Controls.tSteps
            iStep = 1;
            while Kin.T(Temp.iT).t(iStep) < tStep
                iStep = iStep + 1;
            end     
            fprintf(['Plotting Time Step Nb ' num2str(iStep) ', t = ' num2str(Kin.T(Temp.iT).t(iStep)) ' s (' num2str(tStep) ' s)\n'] );


            if strcmp(Syst.Molecule(iMol).KinMthdIn, 'StS')
                LevelPop(:) = Kin.T(Temp.iT).Molecule(iMol).PopOverg(iStep,:);            
            else
                LevelPop(:) = Kin.T(Temp.iT).Molecule(iMol).PopOverg(iStep,LevelToBin(:)) .* Levelq(:) ./ Levelg(:);
            end


            figure(Input.iFig)
            fig = gcf;
            screensize   = get( groot, 'Screensize' );
            %fig.Position = screensize;
            %fig.Color='None';

            
            qq = zeros(Syst.Molecule(iMol).Nvqn + 1,1);
            e0 = zeros(Syst.Molecule(iMol).Nvqn + 1,1); 
            for iLevel = 1:Syst.Molecule(iMol).NLevels
              qq(Levelvqn(iLevel)+1) = qq(Levelvqn(iLevel)+1) + Levelq(iLevel);
              if Syst.Molecule(iMol).Leveljqn(iLevel) == 0
                e0(Levelvqn(iLevel)+1) = Syst.Molecule(iMol).LevelEeV(iLevel);
              end
            end
            ff = zeros(Syst.Molecule(iMol).Nvqn + 1,1);
            ee = zeros(Syst.Molecule(iMol).Nvqn + 1,1);
            for iLevel = 1:Syst.Molecule(iMol).NLevels
              ff(Levelvqn(iLevel)+1) = ff(Levelvqn(iLevel)+1) + LevelPop(iLevel);
              ee(Levelvqn(iLevel)+1) = ee(Levelvqn(iLevel)+1) + LevelEeV(iLevel) * Levelq(iLevel) / qq(Levelvqn(iLevel)+1);
            end
            ff = ff / sum(LevelPop);

            plot(e0,ff,'o')
            hold on
            
            
            xt = get(gca, 'XTick');
            set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
            yt = get(gca, 'YTick');
            set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

            str_x = ['$\epsilon_i$ [eV]'];
            xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
            xlab.Interpreter = 'latex';
            %xlim([max(min(LevelEeV)), MinEvPlot, min(max(LevelEeV)), MaxEvPlot]);

            str_y = ['VDF'];
            ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
            ylab.Interpreter = 'latex';
            %ylim([1.d5, 1.d23]);
            set(gca, 'YScale', 'log')

            
            pbaspect([1 1 1])

            if Input.SaveFigsFlgInt > 0
                [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
                FolderPath = strcat(Input.Paths.SaveFigsFldr, '/T_', Temp.TNowChar, 'K_', Input.Kin.Proc.OverallFlg, '/');
                [status,msg,msgID] = mkdir(FolderPath);
                FileName = strcat('VDF_t', num2str(tStep), 's');
                if Input.SaveFigsFlgInt == 1
                    FileName   = strcat(FolderPath, FileName);
                    export_fig(FileName, '-pdf')
                elseif Input.SaveFigsFlgInt == 2
                    FileName   = strcat(FolderPath, strcat(FileName,'.fig'));
                    savefig(FileName)
                end
                close
            end
            Input.iFig = Input.iFig + 1;


        end
        
        
    end


    fprintf('====================================================\n\n')

end