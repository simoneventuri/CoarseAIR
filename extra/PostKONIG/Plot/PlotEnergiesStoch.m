%% The Function plots the evolutions of the molecular energies  
%
%  Input Arguments:  - iT:                      Index for the current Translational Temperature
%                    - t:                       Vector of time instants
%                    - ProcessesRatesOverall:   Matrix of Processes Overall Rates (Levels X 4(i.e.: Diss,Pair1,Pair2,Pair3) X NTint)
%                    - MolFracs:                Matrix of Mole Fractions (Time-Instants X Components)
%                    - P:                       Vector of Pressures
%                    - eInt, eRot, eVib:        Vectors of Internal, Rotational and Vibrational Energies
%                    ...
%
%  Input Global Var: - XLimPlot:                Vector of [Min, Max] for x-axes plot
%                    - YLimPlot:                Vector of [Min, Max] for y-axes plot
%

function [iFigure, Pressure, tauIntStoch, tauIntPStoch, PDF_tauIntP, tauVibStoch, tauVibPStoch, PDF_tauVibP, tauRotStoch, tauRotPStoch, PDF_tauRotP] = PlotEnergiesStoch(iT, iFigure, tStoch, eIntStoch, eRotStoch, eVibStoch, MolFracsStoch, PStoch, Steps)    

  %%==============================================================================================================
  % 
  % Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
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
  
  global XLimPlot YLimPlot
  global NBinnedMol ColPartToComp T0_Vec linS linST MoleculesName
  global SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz AxisLabelSz AxisLabelNm LegendFontNm 
  global RCVec BCVec GCVec KCVec OCVec PCVec WCVec JCVec YCVec CCVec MCVec
  global iPESEnd iPESStart
  
  NPESTemp = iPESEnd - iPESStart + 1; 


  figure(iFigure)
  fig = gcf;
  screensize = get( groot, 'Screensize' );
  fig.Position=screensize;
  LegendText = [];
  
  for iPES=1:NPESTemp

    Pressure = MolFracsStoch(end,ColPartToComp,iPES+iPESStart-1) * PStoch(end,iPES+iPESStart-1) / 101325;
   

    h1=semilogx(tStoch(Steps(:)),eIntStoch(:,iPES),'k-','LineWidth',1);
    hold on
    h2=semilogx(tStoch(Steps(:)),eRotStoch(:,iPES),'r-','LineWidth',1);
    h3=semilogx(tStoch(Steps(:)),eVibStoch(:,iPES),'g-','LineWidth',1);

    
    eVibLim = (eVibStoch(end,iPES) - eVibStoch(1,iPES)) * 0.632 + eVibStoch(1,iPES);
    iVib=1;
    while eVibStoch(iVib, iPES) < eVibLim
      iVib = iVib+1;
    end
    Deltat = tStoch(Steps(iVib)) - tStoch(Steps(iVib-1));
    Deltae = eVibStoch(iVib,iPES) - eVibStoch(iVib-1,iPES);
    Semie  = eVibLim - eVibStoch(iVib-1,iPES);
    tauVibStoch(iPES)  = Deltat ./ Deltae * Semie + tStoch(Steps(iVib-1));
    tauVibPStoch(iPES) = Pressure * tauVibStoch(iPES);


    eRotLim = (eRotStoch(end,iPES) - eRotStoch(1,iPES)) * 0.632 + eRotStoch(1,iPES);
    iRot=1;
    while eRotStoch(iRot,iPES) < eRotLim
      iRot = iRot+1;
    end
    Deltat = tStoch(Steps(iRot)) - tStoch(Steps(iRot-1));
    Deltae = eRotStoch(iRot,iPES) - eRotStoch(iRot-1,iPES);
    Semie  = eRotLim - eRotStoch(iRot-1,iPES);
    tauRotStoch(iPES)  = Deltat ./ Deltae * Semie + tStoch(Steps(iRot-1));
    tauRotPStoch(iPES) = Pressure * tauRotStoch(iPES);
    
    
    eIntLim = (eIntStoch(end,iPES) - eIntStoch(1,iPES)) * 0.632 + eIntStoch(1,iPES);
    iInt=1;
    while eIntStoch(iInt,iPES) < eIntLim
      iInt = iInt+1;
    end
    Deltat = tStoch(Steps(iInt)) - tStoch(Steps(iInt-1));
    Deltae = eIntStoch(iInt,iPES) - eIntStoch(iInt-1,iPES);
    Semie  = eIntLim - eIntStoch(iInt-1,iPES);
    tauIntStoch(iPES)  = Deltat ./ Deltae * Semie + tStoch(Steps(iInt-1));
    tauIntPStoch(iPES) = Pressure * tauIntStoch(iPES)

  end

  xt = get(gca, 'XTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  yt = get(gca, 'YTick');
  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  str_x = ['Time [s]'];
  xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  xlab.Interpreter = 'latex';
  xlim(XLimPlot);
  str_y = ['Energy [eV]'];
  ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
  ylab.Interpreter = 'latex';
  iFigure=iFigure+1;   
  
  
  
  fig1=figure(iFigure);
  fig = gcf;
  if SaveFigs > 0
    screensize = get( groot, 'Screensize' );
    fig.Position=screensize;
    fig.Color='None';
  end 
  
  h            = histfit(tauVibPStoch,30,'lognormal');
  h(1).FaceAlpha=0.4;
  hold on
  PDF_tauVibP  = fitdist(tauVibPStoch','lognormal');
  
  h            = histfit(tauRotPStoch,30,'lognormal');
  PDF_tauRotP  = fitdist(tauRotPStoch','lognormal');
  
  set(gca, 'XScale', 'log')
  xlab = xlabel('$p \tau [atm*s]$');
  xlab.Interpreter = 'latex';
  ylab = ylabel('Posterior Samples');
  grid on

  set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
  %daspect([1 1 1])
  iFigure = iFigure+1;
  
  PDF_tauIntP = fitdist(tauIntPStoch','lognormal');
  
end