%% The Function plots the Ro-Vibrational Populations at Given Time Steps
%
function Plot_PopulationsVqnSpecific(Controls)    
    
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

    fprintf('= Plot_PopulationsVqnSpecific ========== T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    
    
    for iMol = Controls.MoleculesOI
        fprintf(['Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );
        
        clear LevelToBin Levelvqn LevelEeV LevelPop
        if strcmp(Syst.Molecule(iMol).KinMthdIn, 'StS')
            LevelToBin = Syst.Molecule(iMol).LevelToGroupIn;
        else
            LevelToBin = Syst.Molecule(iMol).LevelToGroupOut;
        end
        Levelvqn   = Syst.Molecule(iMol).Levelvqn;
        LevelEeV   = Syst.Molecule(iMol).LevelEeV;
        iComp      = Syst.MolToCFDComp(iMol);
        
        
        ivVec = [0]
        for iv = ivVec
        
            
            figure(Input.iFig)
            fig = gcf;
            screensize   = get( groot, 'Screensize' );
            %fig.Position = screensize;
            %fig.Color='None';

            jStep = 1;
            for tStep = Controls.tSteps
                iStep = 1;
                while Kin.T(Temp.iT).t(iStep) < tStep
                    iStep = iStep + 1;
                end     
                fprintf(['Plotting Time Step Nb ' num2str(iStep) ', t = ' num2str(Kin.T(Temp.iT).t(iStep)) ' s (' num2str(tStep) ' s)\n'] );


                if strcmp(Syst.Molecule(iMol).KinMthdIn, 'StS')
                    LevelPop(:) = Kin.T(Temp.iT).Molecule(iMol).PopOverg(iStep,:);            
                else
                    LevelPop(:) = Kin.T(Temp.iT).Molecule(iMol).PopOverg(iStep,LevelToBin(:))' .* Syst.Molecule(iMol).T(Temp.iT).Levelq(:) ./ Syst.Molecule(iMol).Levelg(:);
                end


                iivM = 0;
                iivP = Syst.Molecule(iMol).Nvqn-1;
                ColorMat = distinguishable_colors(Syst.Molecule(iMol).Nvqn);

               
                jj = 0;
                for iLevel = 1:Syst.Molecule(iMol).NLevels
                    if Levelvqn(iLevel) == iv
                        jj = jj + 1;
                        LevelEeVTemp(jj) = LevelEeV(iLevel) - LevelEeV(1);% ./ Syst.Molecule(iMol).DissEn;
                        LevelPopTemp(jj) = LevelPop(iLevel) / LevelPop(1);% .* Syst.Molecule(iMol).Levelg(iLevel) ./ Syst.Molecule(iMol).T(Temp.iT).Levelq(iLevel);
                    end
                end
                scatter(LevelEeVTemp, LevelPopTemp, 50, '.', 'MarkerEdgeColor', ColorMat(jStep,:), 'MarkerFaceColor', ColorMat(jStep,:), 'LineWidth', 1.5)
                hold on
                plot(LevelEeVTemp', LevelPopTemp', 'Color', ColorMat(jStep,:), 'LineWidth', 1.5)
                set(gca, 'YScale', 'log')
                clear LevelEeVTemp LevelPopTemp
                
%                 jj = 0;
%                 for iLevel = 1:Syst.Molecule(iMol).NLevels
%                     if Levelvqn(iLevel) == iv+1
%                         jj = jj + 1;
%                         LevelEeVTemp(jj) = LevelEeV(iLevel);
%                         LevelPopTemp(jj) = LevelPop(iLevel);
%                     end
%                 end
%                 scatter(LevelEeVTemp, LevelPopTemp, 10, '^', 'MarkerEdgeColor', ColorMat(jStep,:), 'MarkerFaceColor', ColorMat(jStep,:), 'LineWidth', 1.5)
%                 hold on
%                 plot(LevelEeVTemp', LevelPopTemp', ':', 'Color', ColorMat(jStep,:), 'LineWidth', 1.5)
%                 set(gca, 'YScale', 'log')
%                 clear LevelEeVTemp LevelPopTemp
                
                
                jStep = jStep + 1;
            end
            
            
%             for jv=1:Syst.Molecule(iMol).Nvqn
%               semilogy([Syst.Molecule(iMol).vEeVVib(jv),Syst.Molecule(iMol).vEeVVib(jv)],[1.e0,1e25],'-')  
%             end
                
            xt = get(gca, 'XTick');
            set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
            yt = get(gca, 'YTick');
            set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

            str_x = ['$\epsilon_i$ [eV]'];
            xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
            xlab.Interpreter = 'latex';
            %xlim([max(min(LevelEeV)), MinEvPlot, min(max(LevelEeV)), MaxEvPlot]);

            str_y = ['$N_{i} / g_{i}$ $[m^{-3}]$'];
            ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
            ylab.Interpreter = 'latex';
            ylim([1.d-20, 1.d0]);
            set(gca, 'YScale', 'log')

            pbaspect([1 1 1])

            if Input.SaveFigsFlgInt > 0
                [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr);
                FolderPath = strcat(Input.Paths.SaveFigsFldr, '/T_', Temp.TNowChar, 'K_', Input.Kin.Proc.OverallFlg, '/');
                [status,msg,msgID] = mkdir(FolderPath);
                FileName = strcat(Syst.Molecule(iMol).Name, 'v', num2str(iv));
                if Input.SaveFigsFlgInt == 1
                    FileName   = strcat(FolderPath, FileName);
                    export_fig(FileName, '-pdf');
                elseif Input.SaveFigsFlgInt == 2
                    FileName   = strcat(FolderPath, strcat(FileName,'.fig'));
                    savefig(FileName);
                end
                close
            else
                str_title = [Syst.Molecule(iMol).Name, ', v = ',  num2str(iv)];
                title(str_title, 'interpreter', 'latex');               
            end
            Input.iFig = Input.iFig + 1;

        end

    end


    fprintf('====================================================\n\n')

end