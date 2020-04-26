%% The Function plots the Levels' Populations at different time instants  
%
%  Input Arguments:  - iT:                      Index for the current Translational Temperature
%                    - t:                       Vector of time instants
%                    - ProcessesRatesOverall:   Matrix of Processes Overall Rates (Levels X 4(i.e.: Diss,Pair1,Pair2,Pair3) X NTint)
%                    - MolFracs:                Matrix of Mole Fractions (Time-Instants X Components)
%                    ...
%                    - LevToBin:                (Optional) Mapping Levels-Bins
%
%  Input Global Var: - StepsOverlappingSteps:   Flag 0/1; if =1, different time steps are plotted in the same figure
%                    - MinEvPlot:               Min for x-axes plot
%                    - MaxEvPlot:               Max for x-axes plot
%                    - LevToBinFlg:             Flag 0/1; if =1, populations are colored based on the levels' bins (LevToBin is used)
%                    - vqnColor:                Flag 0/1; if =1, populations are colored based on the levels' vqn (Levelvqn is used)
%

function [iFigure] = PlotPartFunc(iT, iFigure, t, MolFracs, ProcessesRatesOverall, NLevels, LevelEeV, LevelEeV0, LevelQ, Levelg, LevelToBin, Pop, QBins, StpInstants, Levelvqn, Leveljqn, LevToBin, DeltaEintDiss, ProcessesRates)    

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
  
  global MinEvPlot MaxEvPlot StepsOverlappingSteps LevToBinFlg vqnColor KinMthd
  global MoleculesName
  global SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz AxisLabelSz AxisLabelNm LegendFontNm XLimPlot YLimPlot PlotPairUservqnColor 
  global Ue T0_Vec UKb
  
  iMol = 1

  if StepsOverlappingSteps == 1
    figure(iFigure)
    fig = gcf;
    screensize = get( groot, 'Screensize' );
    fig.Position=screensize;
    fig.Color='None';
    LegendText = [];
  end 
  
  ii = 0;
  for iSteps=StpInstants(1:end-1)
    ii = ii+1;
    
    %cmap     = colormap(lines(max(LevToBin(:,iMol,ii))));

    if StepsOverlappingSteps == 0
      figure(iFigure)
      fig = gcf;
      screensize = get( groot, 'Screensize' );
      fig.Position=screensize;
      fig.Color='None';
      title(['Time = ', num2str(t(iSteps)), ' s'])
      hold on
    end

    if (sum(KinMthd(1,:) == 'CGM') == 3) || (sum(KinMthd(1,:) == 'VIB') == 3) 
      LevelPop(1:NLevels(iMol),iMol) = Pop(iSteps,LevelToBin(1:NLevels(iMol),iMol),iMol)' ./ QBins(LevelToBin(1:NLevels(iMol),iMol),iMol) .* exp( - LevelEeV0(1:NLevels(iMol),iMol) .* Ue ./ (T0_Vec(1) .* UKb) );
      LegendPlot(iSteps) = LevelPop(1,iMol);
    else
      LevelPop(1:NLevels(iMol),iMol) = Pop(iSteps,LevelToBin(1:NLevels(iMol),iMol),iMol)' ./ Levelg(LevelToBin(1:NLevels(iMol),iMol),iMol);
      LegendPlot(iSteps) = LevelPop(1,iMol);
    end 
    
    qq = zeros(max(Levelvqn(:,iMol))+1,1);
    e0 = zeros(max(Levelvqn(:,iMol))+1,1); 
    for iLevel = 1:NLevels(iMol)
      qq(Levelvqn(iLevel,iMol)+1) = qq(Levelvqn(iLevel,iMol)+1) + Levelg(iLevel,iMol) * exp( - LevelEeV0(iLevel,iMol) .* Ue ./ (T0_Vec(1) .* UKb) );
      if Leveljqn(iLevel,iMol) == 0
        e0(Levelvqn(iLevel,iMol)+1) = LevelEeV(iLevel,iMol);
      end
    end
    ff = zeros(max(Levelvqn(:,iMol))+1,1);
    ee = zeros(max(Levelvqn(:,iMol))+1,1);
    for iLevel = 1:NLevels(iMol)
      ff(Levelvqn(iLevel,iMol)+1) = ff(Levelvqn(iLevel,iMol)+1) + Pop(iSteps,iLevel,iMol);
      ee(Levelvqn(iLevel,iMol)+1) = ee(Levelvqn(iLevel,iMol)+1) + LevelEeV(iLevel,iMol) * Levelg(iLevel,iMol) * exp( - LevelEeV0(iLevel,iMol) .* Ue ./ (T0_Vec(1) .* UKb) ) / qq(Levelvqn(iLevel,iMol)+1);
    end
    ff = ff / sum(Pop(iSteps,:,1));
    
    plot(e0,ff,'o')
    hold on
    
end