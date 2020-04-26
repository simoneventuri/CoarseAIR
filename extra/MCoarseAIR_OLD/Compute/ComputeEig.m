%% The Function computes the eigen-functions of the rates matrix and plots the evolutions of the distribution functions caused by ...
%      Inelasti Processes
%
%  Input Arguments:  - iFigure              Index for the first figure to plot
%                    - RatesMatrix          Matrix of Rates
%                    - NLevels              Number of Molecules Levels
%                    - Levelg               Levels Degeneracies
%                    - LevelEeV             Levels energy in eV, with the 0 corresponding to the dissociation energy
%   
%  Input Global Var: - ComputeEigFlg        Flag 1/2; if 1, the Eigen-functions are computed; if 2, they are red from .mat file
%                    - tInstants            Vector of Time Instants for computing and plotting levels' populations 
%                    - PairsToSum           Vector of Pair Indexes of Inelastic Exchange (e.g.: CO+O = [2]) 
%

function [iFigure, EigenVals] = ComputeEig(iFigure, iT, RatesMatrix, NLevels, Levelg, LevelEeV)

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
  
  global ComputeEigFlg NBins RatesPath T0_Vec Ue UKb tInstants PairsToSum
  global SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz LegendFontNm
  
  T0 = 300;
  
  
  iBinnedMol = 1;
  Pop0(1:NLevels(iBinnedMol),iBinnedMol) = Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* exp( - LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) .* Ue ./ (T0 .* UKb) );  
  Pop0(1:NLevels(iBinnedMol),iBinnedMol) = Pop0(1:NLevels(iBinnedMol),iBinnedMol) ./ Levelg(1:NLevels(iBinnedMol),iBinnedMol) ./ sum(Pop0(1:NLevels(iBinnedMol),iBinnedMol));
  sum(Pop0(1:NLevels(iBinnedMol),iBinnedMol))
  
  fig1=figure(iFigure);
  fig = gcf;
  if SaveFigs > 0
    screensize = get( groot, 'Screensize' );
    fig.Position=screensize;
    fig.Color='None';
  end 
  
  semilogy(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),Pop0(1:NLevels(iBinnedMol),iBinnedMol));
  
  xlab = xlabel('Energy [eV]');
  xlab.Interpreter = 'latex';
  ystr = strcat('$N_i / (N * g_i)$')
  ylab = ylabel(ystr);
  ylab.Interpreter = 'latex';
  
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  %daspect([1 1 1])

  if SaveFigs == 1
     FolderPath = strcat(FigDirPath, '/Eigen/');
     [status,msg,msgID] = mkdir(FolderPath);
     FileName   = strcat(FolderPath, 'Eigen');
     export_fig(FileName, '-pdf')
     close
  elseif SaveFigs == 2
     FolderPath = strcat(FigDirPath, '/Eigen/');
     [status,msg,msgID] = mkdir(FolderPath);
     FileName   = strcat(FolderPath, 'Eigen.fig');
     saveas(fig,FileName);
  end
  iFigure = iFigure + 1;
  
  
  if ComputeEigFlg == 1
    
    InelProcesses(1:NBins(1),1:NBins(1)) = RatesMatrix(1:NBins(1),1:NBins(1),1,iT);
    for iP = PairsToSum
      InelProcesses(1:NBins(1),1:NBins(1)) = InelProcesses(1:NBins(1),1:NBins(1)) + RatesMatrix(1:NBins(1),1:NBins(1),iP,iT);
    end 
    
    [REigenVect,EigenVals] = eig(InelProcesses);
  
    filenameEigenVals = strcat(RatesPath,'/T_',num2str(T0_Vec(iT)),'_',num2str(T0_Vec(iT)),'/Eigen')
    save(filenameEigenVals,'EigenVals','REigenVect','-v7.3');
    
  else
    
    filenameEigenVals = strcat(RatesPath,'/T_',num2str(T0_Vec(iT)),'_',num2str(T0_Vec(iT)),'/Eigen.mat')
    load(filenameEigenVals,'EigenVals','REigenVect');
    
  end
  
  fig1=figure(iFigure);
  fig = gcf;
  if SaveFigs > 0
    screensize = get( groot, 'Screensize' );
    fig.Position=screensize;
    fig.Color='None';
  end 
  
  plot(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),real(EigenVals),'ko')
  
  xlab = xlabel('Level Index');
  xlab.Interpreter = 'latex';
  ystr = strcat('Eigenvalue')
  ylab = ylabel(ystr);
  ylab.Interpreter = 'latex';
  
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  %daspect([1 1 1])

  if SaveFigs == 1
     FolderPath = strcat(FigDirPath, '/Eigen/');
     [status,msg,msgID] = mkdir(FolderPath);
     FileName   = strcat(FolderPath, 'EigenValues');
     export_fig(FileName, '-pdf')
     close
  elseif SaveFigs == 2
     FolderPath = strcat(FigDirPath, '/Eigen/');
     [status,msg,msgID] = mkdir(FolderPath);
     FileName   = strcat(FolderPath, 'EigenValues.fig');
     saveas(fig,FileName);
  end
  iFigure = iFigure + 1;
  
  
  
  LEigenVect = inv(REigenVect);
  
  
  
  for tInst = tInstants
  
    
    Popt = REigenVect * exp(EigenVals .* tInst) * LEigenVect * Pop0;

    
    fig1=figure(iFigure);
    fig = gcf;
    if SaveFigs > 0
      screensize = get( groot, 'Screensize' );
      fig.Position=screensize;
      fig.Color='None';
    end 
    semilogy(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol),Popt(1:NLevels(iBinnedMol),iBinnedMol),'o');
    
    xlab = xlabel('Energy [eV]');
    xlab.Interpreter = 'latex';
    ystr = strcat('$N_i / (N * g_i)$')
    ylab = ylabel(ystr);
    ylab.Interpreter = 'latex';

    set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
    %daspect([1 1 1])

    if SaveFigs == 1
       FolderPath = strcat(FigDirPath, '/Eigen/');
       [status,msg,msgID] = mkdir(FolderPath);
       FileName   = strcat(FolderPath, 'DistFunc-', num2str(tInst));
       export_fig(FileName, '-pdf')
       close
    elseif SaveFigs == 2
       FolderPath = strcat(FigDirPath, '/Eigen/');
       [status,msg,msgID] = mkdir(FolderPath);
       FileName   = strcat(FolderPath, 'DistFunc-', num2str(tInst), '.fig');
       saveas(fig,FileName);
    end
    iFigure = iFigure + 1;

    
  end
  
  
  
end 