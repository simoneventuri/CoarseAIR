%% The Function plots the Ro-Vibrational Populations at Given Time Steps
%
function Plot_DiatPot(Controls)    
    
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

    fprintf('= Plot_DiatPot =====================================\n')
    fprintf('====================================================\n')
    
    
    for iMol = 1:Controls.MoleculesOI
        fprintf(['Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );

               
%         figure(Input.iFig)
%         fig = gcf;
%         screensize   = get( groot, 'Screensize' );
%         %fig.Position = screensize;
%         %fig.Color='None';
%         
%         
%         rVec = linspace(Controls.Extremes(iMol,1), Controls.Extremes(iMol,2), 10000);
%         
%         ij = 0;
%         for jqn = Controls.jqnVec
%             jqn
%             ij = ij + 1;
%             
%             [Ve, dVe] = DiatPot(rVec, jqn, iMol);
%             
%             plot(rVec, Ve, '-', 'Color', Param.CMat(ij,:), 'LineWidth', Param.LineWidth);
%             hold on
%         end
%             
% 
%         xt = get(gca, 'XTick');
%         set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
%         yt = get(gca, 'YTick');
%         set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
% 
%         str_x = ['r [$a_0$]'];
%         xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
%         xlab.Interpreter = 'latex';
%         %xlim([max(min(LevelEeV)), MinEvPlot, min(max(LevelEeV)), MaxEvPlot]);
% 
%         str_y = ['Energy [eV]'];
%         ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
%         ylab.Interpreter = 'latex';
%         %ylim([1.d5, 1.d23]);
%         %set(gca, 'YScale', 'log')
% 
% 
%         pbaspect([1 1 1])
% 
%         if Input.SaveFigsFlgInt > 0
%             [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
%             FolderPath = strcat(Input.Paths.SaveFigsFldr, '/T_', Temp.TNowChar, 'K_', Input.Kin.Proc.OverallFlg, '/');
%             [status,msg,msgID] = mkdir(FolderPath);
%             FileName = strcat(Syst.Molecule(iMol).Name,'_Energies');
%             if Input.SaveFigsFlgInt == 1
%                 FileName   = strcat(FolderPath, FileName);
%                 export_fig(FileName, '-pdf')
%             elseif Input.SaveFigsFlgInt == 2
%                 FileName   = strcat(FolderPath, strcat(FileName,'.fig'));
%                 savefig(FileName)
%             end
%             close
%         end
%         Input.iFig = Input.iFig + 1;
        
        
        
        [status,msg,msgID] = mkdir(Input.Paths.SaveDataFldr);
        FileName           = strcat(Input.Paths.SaveDataFldr, '/DiatPot_', Syst.Molecule(iMol).Name, '.csv');
        fileID = fopen(FileName,'w');
        fprintf(fileID,'id,v,J,EeV,rIn,rOut,EeVVib,EeVRot,dCentBarr,Bin\n');
        for iLevel = 1:Syst.Molecule(iMol).NLevels
            
            
            fprintf(fileID,'%i,%i,%i,%e,%e,%e,%e,%e,%e,%i\n',   iLevel, 	                              ...
                                                                Syst.Molecule(iMol).Levelvqn(iLevel),     ...
                                                                Syst.Molecule(iMol).Leveljqn(iLevel),     ...
                                                                Syst.Molecule(iMol).LevelEeV(iLevel),     ...
                                                                Syst.Molecule(iMol).LevelrIn(iLevel),     ...
                                                                Syst.Molecule(iMol).LevelrOut(iLevel),    ...
                                                                Syst.Molecule(iMol).LevelEeVVib0(iLevel), ...
                                                                Syst.Molecule(iMol).LevelEeVRot(iLevel),  ...
                                                                Syst.Molecule(iMol).LevelECB(iLevel),     ...
                                                                Syst.Molecule(iMol).LevelToGroupOut(iLevel)  ...
                                                                );
                                                                   
        end
        fclose(fileID);
        
        
        [status,msg,msgID] = mkdir(Input.Paths.SaveDataFldr);
        FileName           = strcat(Input.Paths.SaveDataFldr, '/DiatPot_', Syst.Molecule(iMol).Name, '_CentBarr.csv');
        fileID = fopen(FileName,'w');
        fprintf(fileID,'id,v,J,EeV,rIn,rOut,EeVVib,EeVRot,dCentBarr,VMax,rMax\n');
        for iLevel = 1:Syst.Molecule(iMol).NLevels
            
            if Syst.Molecule(iMol).Levelvqn(iLevel) == 0

                fprintf(fileID,'%i,%i,%i,%e,%e,%e,%e,%e,%e,%e,%e\n',  iLevel, 	                          ...
                                                                Syst.Molecule(iMol).Levelvqn(iLevel),     ...
                                                                Syst.Molecule(iMol).Leveljqn(iLevel),     ...
                                                                Syst.Molecule(iMol).LevelEeV(iLevel),     ...
                                                                Syst.Molecule(iMol).LevelrIn(iLevel),     ...
                                                                Syst.Molecule(iMol).LevelrOut(iLevel),    ...
                                                                Syst.Molecule(iMol).LevelEeVVib0(iLevel), ...
                                                                Syst.Molecule(iMol).LevelEeVRot(iLevel),  ...
                                                                Syst.Molecule(iMol).LevelECB(iLevel),     ...
                                                                Syst.Molecule(iMol).LevelVMax(iLevel),     ...
                                                                Syst.Molecule(iMol).LevelrMax(iLevel)     ...
                                                                );

            end
            
        end
        fclose(fileID);
        
    end
    
    
    fprintf('====================================================\n\n')
    
end