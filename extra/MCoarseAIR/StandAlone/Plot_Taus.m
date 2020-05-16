%% The Function plots the Relaxation Time Constants as Function of Translational Temperature
%
  
    
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


clear all
%close all
clc

global Input Syst Temp Param Kin Rates



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SPECIFYING INPUT 

%% System Inputs
Input.Paths.ToQCTFldr       = '/home/venturi/WORKSPACE/CoarseAIR/CO2_ALL/Test/';
Input.Paths.ToKinMainFldr   = '/home/venturi/WORKSPACE/Mars_Database/Run_0D/';
Input.Paths.ToHDF5Fldr      = '/home/venturi/WORKSPACE/Mars_Database/HDF5_Database/';
Input.TranVec               = [12500 15000 20000];
Input.SystNameLong          = 'CO2_NASA';
Input.iPES                  = 0;
Input.Kin.MolResolutionIn   = ['StS'; 'StS'];
Input.Kin.MinStateIn        = [    1,     1];
Input.Kin.MaxStateIn        = [13521,  6078];
Input.Kin.NGroupsIn         = [    0,     0];
Input.Kin.Proc.DissFlg      = 0;
Input.Kin.DissCorrFactor    = 1.0;
Input.Kin.Proc.DissInelFlg  = 0;
Input.Kin.Proc.InelFlg      = 1;
Input.Kin.Proc.ExchFlg1     = 0;
Input.Kin.Proc.ExchFlg2     = 0;
Input.Kin.RateSource        = 'HDF5'; % CoarseAIR / CG-QCT / HDF5 / PLATO
Input.FigureFormat          = 'PrePrint';
Input.ReLoad                = 1;


%% Inputs for Plotting
Input.iFig               = 201;
Input.SaveFigsFlgInt     = 0;
Input.Paths.SaveFigsFldr = '/home/venturi/WORKSPACE/CO2_Paper/Figures/';


%% Inputs for Saving Data
Input.Paths.SaveDataFldr = '/home/venturi/WORKSPACE/CO2_Paper/Data/';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initializing
Initialize_ChemicalSyst()
Initialize_Input()
Initialize_Parameters()



fprintf('= Plot_Taus ========================================\n')
fprintf('====================================================\n')


for iMol = 1:1%Syst.NMolecules
    fprintf(['Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );

    
    
    opts = delimitedTextImportOptions("NumVariables", 5);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["T", "P", "tauIntP", "tauRotP", "tauVibP"];
    opts.VariableTypes = ["double", "double", "double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    FileName = strcat(Input.Paths.SaveDataFldr, '/Taus_', Syst.Molecule(iMol).Name, '_', Input.Kin.Proc.OverallFlg, '.csv');
    tbl = readtable(FileName, opts);
    Vec.T       = tbl.T;
    Vec.P       = tbl.P;
    Vec.tauIntP = tbl.tauIntP;
    Vec.tauRotP = tbl.tauRotP;
    Vec.tauVibP = tbl.tauVibP;
    clear opts tbl
    
    xx=Vec.T
    yy=Vec.tauVibP


    figure(Input.iFig)
    fig = gcf;
    screensize   = get( groot, 'Screensize' );
    %fig.Position = screensize;
    %fig.Color='None';


    semilogy(Vec.T.^(-1/3), Vec.tauVibP, '-', 'Color', Param.CMat(1,:), 'LineWidth', Param.LineWidth);
    hold on


    xt = get(gca, 'XTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca,'FontSize', Param.AxisFontSz, 'FontName', Param.AxisFontNm, 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

    str_x = ['T [K]'];
    xlab             = xlabel(str_x, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    xlab.Interpreter = 'latex';
    %xlim([max(min(LevelEeV)), MinEvPlot, min(max(LevelEeV)), MaxEvPlot]);

    str_y = ['$P*\tau_V$ [atm*s]'];
    ylab             = ylabel(str_y, 'Fontsize', Param.AxisLabelSz, 'FontName', Param.AxisLabelNm);
    ylab.Interpreter = 'latex';
    %ylim([1.d5, 1.d23]);
    %set(gca, 'YScale', 'log')


    pbaspect([1 1 1])

    if Input.SaveFigsFlgInt > 0
        [status,msg,msgID]  = mkdir(Input.Paths.SaveFigsFldr)
        FolderPath = strcat(Input.Paths.SaveFigsFldr, '/', Input.Kin.Proc.OverallFlg, '/');
        [status,msg,msgID] = mkdir(FolderPath);
        FileName = strcat(Syst.Molecule(iMol).Name,'_Taus');
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


fprintf('====================================================\n\n')