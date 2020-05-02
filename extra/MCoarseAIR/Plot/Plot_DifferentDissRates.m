%% The Function plots the Contributions to the Dissociation Rate Coefficients from the Different Pairs
%
function Plot_DifferentDissRates()    
    
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

    fprintf('= Plot_DifferentDissRates ============== T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')

    iMol       = 1;
    LevelToBin = Syst.Molecule(iMol).LevelToBin;
    Levelvqn   = Syst.Molecule(iMol).Levelvqn;
    LevelEeV   = Syst.Molecule(iMol).LevelEeV;
    
    
    figure(Input.iFig)
    fig = gcf;
    screensize   = get( groot, 'Screensize' );
    %fig.Position = screensize;
    %fig.Color='None';
    
    ProcNames = {};
    for iP = 1:3
        scatter(LevelEeV, Rates.T(Temp.iT).Diss(:,iP+1), 20, 'Filled', 'MarkerFaceColor', Param.CMat(iP,:));
        hold on
        ProcNames = [ProcNames, strcat('Pair 1,{ }', Syst.Molecule(Syst.Pair(iP).ToMol).Name)];
    end
    

    clab             = legend(ProcNames, 'Location', 'Best');
    clab.Interpreter = 'latex';
    set(clab,'FontSize', Param.LegendFontSz, 'FontName', Param.LegendFontNm, 'Interpreter', 'latex');    

    xt = get(gca, 'XTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

    str_x = ['$\epsilon_i$ [eV]'];
    xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    xlab.Interpreter = 'latex';
    %xlim([max(min(LevelEeV)), MinEvPlot, min(max(LevelEeV)), MaxEvPlot]);

    str_y = ['$k_i^D$ $[cm^{3}/s]$'];
    ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    ylab.Interpreter = 'latex';
    %ylim([1.d5, 1.d23]);
    set(gca, 'YScale', 'log')


    pbaspect([1 1 1])

    if Input.SaveFigsFlgInt > 0
        [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
        FolderPath = strcat(Input.Paths.SaveFigsFldr, '/T_', Temp.TNowChar, 'K_', Input.Kin.Proc.OverallFlg, '/');
        [status,msg,msgID] = mkdir(FolderPath);
        FileName = strcat('Pops_t', num2str(tStep), 's');
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

    
    fprintf('====================================================\n\n')

end