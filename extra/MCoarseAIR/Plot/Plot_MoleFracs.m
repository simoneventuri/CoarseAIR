%% The Function plots the Mole Fractions of the Chemical System's Components 
%
function Plot_MoleFracs(Controls)    
    
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

    fprintf('= Plot_MoleFracs ======================= T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    
    
    figure(Input.iFig)
    fig = gcf;
    screensize   = get( groot, 'Screensize' );
    %fig.Position = screensize;
    %fig.Color='None';
    
    if Controls.Normalize(1) == 0
        for iComp = Controls.CompStart:Controls.CompEnd
          semilogx(Kin.T(Temp.iT).t(:), Kin.T(Temp.iT).MolFracs(:,iComp), 'Color', Syst.CFDComp(iComp).Color, 'linestyle', Syst.CFDComp(iComp).LineStyle, 'LineWidth', Param.LineWidth)
          hold on
        end
        hold on

        if (Input.Tasks.Plot_Populations.Flg)
            jStep = 1;
            for tStep = Input.Tasks.Plot_Populations.tSteps
                iStep = 1;
                while Kin.T(Temp.iT).t(iStep) < tStep
                    iStep = iStep + 1;
                end  
                Controls.iSteps(jStep) = iStep;
                jStep = jStep + 1;
            end
            Controls.iSteps(jStep) = Kin.T(Temp.iT).QSS.i;

            for iMol = Input.Tasks.Plot_Populations.MoleculesOI
                iComp = Syst.MolToCFDComp(iMol);
                semilogx(Kin.T(Temp.iT).t(Controls.iSteps), Kin.T(Temp.iT).MolFracs(Controls.iSteps,iComp), 'o', 'MarkerFaceColor', Syst.CFDComp(iComp).Color, 'MarkerEdgeColor', Syst.CFDComp(iComp).Color)
            end
        end

        xt = get(gca, 'XTick');
        set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
        yt = get(gca, 'YTick');
        set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

        clab             = legend(Syst.CFDComp(Controls.CompStart:Controls.CompEnd).Name, 'Location', 'Best');
        clab.Interpreter = 'latex';
        set(clab,'FontSize', Param.LegendFontSz, 'FontName', Param.LegendFontNm, 'Interpreter', 'latex');

        str_x = ['t [s]'];
        xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
        xlab.Interpreter = 'latex';
        %xlim(XLimPlot);

        str_y = ['Mole Fraction'];
        ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
        ylab.Interpreter = 'latex';
        %ylim(YLimPlot);

        pbaspect([1 1 1])
        
    else
        
        
        for iComp = Controls.Normalize
          semilogx(Kin.T(Temp.iT).t(:), Kin.T(Temp.iT).MolFracs(:,iComp)./Kin.T(Temp.iT).MolFracs(1,iComp), 'Color', Syst.CFDComp(iComp).Color, 'linestyle', Syst.CFDComp(iComp).LineStyle, 'LineWidth', Param.LineWidth)
          hold on
        end
        hold on

        xt = get(gca, 'XTick');
        set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
        yt = get(gca, 'YTick');
        set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

        clab             = legend(Syst.CFDComp(Controls.CompStart:Controls.CompEnd).Name, 'Location', 'Best');
        clab.Interpreter = 'latex';
        set(clab,'FontSize', Param.LegendFontSz, 'FontName', Param.LegendFontNm, 'Interpreter', 'latex');

        str_x = ['t [s]'];
        xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
        xlab.Interpreter = 'latex';
        %xlim(XLimPlot);

        str_y = ['[X]/[X]${}_0$'];
        ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
        ylab.Interpreter = 'latex';
        %ylim(YLimPlot);

        pbaspect([1 1 1])
        
    end


    if Input.SaveFigsFlgInt > 0
        [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr);
        FolderPath = strcat(Input.Paths.SaveFigsFldr, '/T_', Temp.TNowChar, 'K_', Input.Kin.Proc.OverallFlg, '/');
        [status,msg,msgID] = mkdir(FolderPath);
        if Input.SaveFigsFlgInt == 1
            FileName   = strcat(FolderPath, 'MoleFractions');
            export_fig(FileName, '-pdf');
        elseif Input.SaveFigsFlgInt == 2
            FileName   = strcat(FolderPath, 'MoleFractions.fig');
            savefig(FileName);
        end
        close
    end
    Input.iFig = Input.iFig + 1;

    
    fprintf('====================================================\n\n')

end