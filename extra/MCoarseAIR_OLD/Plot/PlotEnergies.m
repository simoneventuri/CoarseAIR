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

function [iFigure, tauInt, tauIntP, tauVib, tauVibP, tauRot, tauRotP] = PlotEnergies(iT, iFigure, t, eInt, eRot, eVib, MolFracs, P, Steps)    

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
  
  global XLimPlot YLimPlot
  global NBinnedMol ColPartToComp T0_Vec linS linST MoleculesName
  global SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz AxisLabelSz AxisLabelNm LegendFontNm 
  global RCVec BCVec GCVec KCVec OCVec PCVec WCVec JCVec YCVec CCVec MCVec
  global DissFlg InelFlg ExchFlg
  
  for iMol = 1:1 %NBinnedMol
    
    Pressure = MolFracs(end,ColPartToComp) * P(end) / 101325
    
    figure(iFigure)
    fig = gcf;
    screensize = get( groot, 'Screensize' );
    fig.Position=screensize;
    fig.Color='None';
    LegendText = [];

    h1=semilogx(t(Steps),eInt(Steps,iMol),'-','Color',KCVec,'LineWidth',3);
    hold on
    h2=semilogx(t(Steps),eRot(Steps,iMol),'-','Color',RCVec,'LineWidth',3);
    h3=semilogx(t(Steps),eVib(Steps,iMol),'-','Color',GCVec,'LineWidth',3);
    
    eVibLim = (eVib(Steps(end),iMol) - eVib(Steps(1),iMol)) * 0.632 + eVib(Steps(1),iMol);
    %semilogx([t(2),t(end)],[eVibLim, eVibLim],'-','Color','g','LineWidth',2,'linestyle',linST{iT});
    iVib=1;
    while eVib(Steps(iVib)) < eVibLim
      iVib = iVib+1;
    end
    tauVib(iMol)   = (t(Steps(iVib)) + t(Steps(iVib-1)))./2.d0;
%     eVibInfo       = stepinfo(eVib(Steps),t(Steps));
%     tauVib(iMol)   = eVibInfo.RiseTime;
    eVibLT(Steps)  = eVib(Steps(end),iMol) - (eVib(Steps(end),iMol)- eVib(Steps(1),iMol)) .* exp( -t(Steps) ./ tauVib(iMol));
    tauVibP(iMol)  = Pressure * tauVib(iMol)
    h6 = semilogx(t(Steps),eVibLT(Steps),':','Color',GCVec,'LineWidth',3);
    
    
    eIntLim = (eInt(Steps(end),iMol) - eInt(Steps(1),iMol)) * 0.632 + eInt(Steps(1),iMol);
    %semilogx([t(2),t(end)],[eIntLim, eIntLim],'-','Color','k','LineWidth',2,'linestyle',linST{iT});
    iInt=1;
    while eInt(Steps(iInt)) < eIntLim
      iInt = iInt+1;
    end
    tauInt(iMol)   = (t(Steps(iInt)) + t(Steps(iInt-1)))./2.d0;
%     eIntInfo       = stepinfo(eInt(Steps),t(Steps));
%     tauInt(iMol)   = eIntInfo.RiseTime;
    eIntLT(Steps)  = eInt(Steps(end),iMol) - (eInt(Steps(end),iMol)- eInt(Steps(1),iMol)) .* exp( -t(Steps) ./ tauInt(iMol));
    tauIntP(iMol)  = Pressure * tauInt(iMol)
    h4 = semilogx(t(Steps),eIntLT(Steps),':','Color',KCVec,'LineWidth',3);
    
    
    eRotLim = (eRot(Steps(end),iMol) - eRot(Steps(1),iMol)) * 0.632 + eRot(Steps(1),iMol);
    %semilogx([t(2),t(end)],[eRotLim, eRotLim],'-','Color','r','LineWidth',2,'linestyle',linST{iT});
    iRot=1;
    while eRot(Steps(iRot)) < eRotLim
      iRot = iRot+1;
    end
    if iRot~=1 
      tauRot(iMol)   = (t(Steps(iRot)) + t(Steps(iRot-1)))./2.d0;
  %     eRotInfo       = stepinfo(eRot(Steps),t(Steps));
  %     tauRot(iMol)   = eRotInfo.RiseTime;
    else
      tauRot(iMol)   = 0.0
    end
    eRotLT(Steps)  = eRot(Steps(end),iMol) - (eRot(Steps(end),iMol)- eRot(Steps(1),iMol)) .* exp( -t(Steps) ./ tauRot(iMol));
    tauRotP(iMol)  = Pressure * tauRot(iMol)
    h5 = semilogx(t(Steps),eRotLT(Steps),':','Color',RCVec,'LineWidth',3);
    
    title(['T^{-1/3} = ', num2str(T0_Vec(iT)^(-1/3)), 'K;  P \tau_{int} =', num2str(tauIntP(iMol)), ' atm*s;  P \tau_{vib} =', num2str(tauVibP(iMol)), ' atm*s; P \tau_{rot} =', num2str(tauRotP(iMol)), ' atm*s']);
    legend([h1,h4,h2,h5,h3,h6],'Internal Energy from Master Equation','Internal Energy from Landau Teller','Rotational Energy from Master Equation','Rotational Energy from Landau Teller','Vibrational Energy from Master Equation','Vibrational Energy from Landau Teller');

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
    %ylim(YLimPlot);

    if SaveFigs == 1
      FolderPath = strcat(FigDirPath, '/FlowQuantities/');
      [status,msg,msgID] = mkdir(FolderPath);
      FileName   = strcat(FolderPath, MoleculesName(1,:),'-Energies');
      export_fig(FileName, '-pdf')
      close
    elseif SaveFigs == 2
      FolderPath = strcat(FigDirPath, '/FlowQuantities/');
      [status,msg,msgID] = mkdir(FolderPath);
      FileName   = strcat(FolderPath, MoleculesName(1,:),'-Energies.fig');
      savefig(FileName)
      close
    end

    iFigure=iFigure+1;   
    
                        
    FileName = strcat('./Taus_',num2str(DissFlg),'_',num2str(InelFlg),'_',num2str(ExchFlg),'_',num2str(ExchFlg),'.csv');
    if exist(FileName, 'file')
      fileID1  = fopen(FileName,'a');
    else
      fileID1  = fopen(FileName,'w');
      fprintf(fileID1,'# T [K], P [atm], TauInt, TauVib, TauRot\n');
    end
    fprintf(fileID1,'%e,%e,%e,%e,%e\n', T0_Vec(1), Pressure, tauIntP, tauVibP, tauRotP)
    fclose(fileID1);

    
    

%     eIntOrig = eInt;
%     eRotOrig = eRot;
%     eVibOrig = eVib;
%     %save('./TempEnergies', 'eInt', 'eRot', 'eVib', '-v7.3')
%     %load('./TempEnergies', 'eInt', 'eRot', 'eVib')
%     %load('./InelExchDiss', 'eInt', 'eRot', 'eVib')
%     eIntTemp = eInt;
%     eRotTemp = eRot;
%     eVibTemp = eVib;
%     eIntVar = (eIntTemp - eIntOrig) ./ eIntOrig .* 1.d2;
%     eRotVar = (eRotTemp - eRotOrig) ./ eRotOrig .* 1.d2;
%     eVibVar = (eVibTemp - eVibOrig) ./ eVibOrig .* 1.d2;
%     
%     figure(iFigure)
%     fig = gcf;
%     screensize = get( groot, 'Screensize' );
%     fig.Position=screensize;
%     fig.Color='None';
%     LegendText = [];
% 
%     h1=semilogx(t(Steps),eIntVar(Steps,iMol),'-','Color',KCVec,'LineWidth',3);
%     hold on
%     h2=semilogx(t(Steps),eRotVar(Steps,iMol),'-','Color',RCVec,'LineWidth',3);
%     h3=semilogx(t(Steps),eVibVar(Steps,iMol),'-','Color',GCVec,'LineWidth',3);
%     
%     xt = get(gca, 'XTick');
%     set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
%     yt = get(gca, 'YTick');
%     set(gca,'FontSize',AxisFontSz, 'FontName',AxisFontNm,'TickDir','out','TickLabelInterpreter', 'latex');
% 
%     str_x = ['Time [s]'];
%     xlab = xlabel(str_x,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
%     xlab.Interpreter = 'latex';
%     xlim(XLimPlot);
% 
%     str_y = ['Energy [eV]'];
%     ylab = ylabel(str_y,'Fontsize',AxisLabelSz,'FontName',AxisLabelNm);
%     ylab.Interpreter = 'latex';
%     %ylim(YLimPlot);
% 
%     if SaveFigs == 1
%       FolderPath = strcat(FigDirPath, '/FlowQuantities/');
%       [status,msg,msgID] = mkdir(FolderPath);
%       FileName   = strcat(FolderPath, MoleculesName(1,:),'-Energies');
%       export_fig(FileName, '-pdf')
%       close
%     elseif SaveFigs == 2
%       FolderPath = strcat(FigDirPath, '/FlowQuantities/');
%       [status,msg,msgID] = mkdir(FolderPath);
%       FileName   = strcat(FolderPath, MoleculesName(1,:),'-Energies.fig');
%       savefig(FileName)
%       close
%     end
% 
%     iFigure=iFigure+1;   
    
  end
    
end