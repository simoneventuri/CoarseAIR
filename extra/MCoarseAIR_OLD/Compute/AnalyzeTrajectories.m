% -- MATLAB --
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

function [iFigure, xxFin, EintFin, RbondFin, AngMomFin, jqnFin, vqnFin, iTypeFin] = AnalyzeTrajectories(iFigure, iTint, iBins, Levelvqn, Leveljqn, LevelEeV)
  
  global SaveFigs AtomMass Pair_to_Atoms T0_Vec iNode iProc System
    
  filename = strcat('../Test/T_',num2str(T0_Vec(iTint)),'_',num2str(T0_Vec(iTint)),'/PaQSOl-Tot.out')
  %filename = strcat('../Test/T_',num2str(T0_Vec(iTint)),'_',num2str(T0_Vec(iTint)),'/Bins_',num2str(iBins),'_0/PaQSol.out')
  startRow = 2;
  formatSpec = '%*35s%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%18f%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  PaQSol    = [dataArray{1:end-1}];
  clearvars formatSpec fileID dataArray ans;
  formatSpec = '%17f%18f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  TrajIndx = dataArray{:, 1};
  tFin     = dataArray{:, 2};
  clearvars filename startRow formatSpec fileID dataArray ans;
  HIni(:)        = PaQSol(:,1);
  PaQIni(:,1:12) = PaQSol(:,2:13);
  HFin(:)        = PaQSol(:,14);
  PaQFin(:,1:12) = PaQSol(:,15:26);
  
  
  figure(iFigure)
  h1=plot(TrajIndx,HIni,'go');
  hold on
  h2=plot(TrajIndx,HFin,'ro');
  legend([h1,h2],'H_{Initial}','H_{Final}')
  set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
  set(gcf, 'PaperPositionMode', 'auto');
  if SaveFigs == 1
      FileName = strcat(MoleculesName(1,:),'-Dissociation');
      FilePathNamePng = strcat(FilePath, FileName, '.png');
      FilePathNameFig = strcat(FilePath, FileName, '.fig');
      saveas(gcf,strcat(FilePath, FileName),'pdf')
      saveas(gcf,FilePathNamePng)
      savefig(FilePathNameFig);
  end
  xlabel('Trajectory Indx')
  ylabel('H [Eh]');
  iFigure=iFigure+1;
  
  figure(iFigure)
  histogram(HFin,50,'Normalization','probability');
  set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
  set(gcf, 'PaperPositionMode', 'auto');
  if SaveFigs == 1
      FileName = strcat(MoleculesName(1,:),'-Dissociation');
      FilePathNamePng = strcat(FilePath, FileName, '.png');
      FilePathNameFig = strcat(FilePath, FileName, '.fig');
      saveas(gcf,strcat(FilePath, FileName),'pdf')
      saveas(gcf,FilePathNamePng)
      savefig(FilePathNameFig);
  end
  xlabel('H_{Final}')
  ylabel('PDF');
  iFigure=iFigure+1;
  
  
  filename = '/Users/sventuri/Downloads/elevels_g (1).dat';
  formatSpec = '%*39s%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN, 'ReturnOnError', false);
  fclose(fileID);
  LevToBin = dataArray{:, 1};
  clearvars filename formatSpec fileID dataArray ans;
  LevToBin = LevToBin + 1;
  NBins = max(LevToBin);
  for i = length(LevToBin)+1:size(Leveljqn,1)
    LevToBin(i) = NBins+1;
  end
  for i = 1:length(LevToBin)
    qnToBin(Levelvqn(i,1)+1,Leveljqn(i,1)+1) = LevToBin(i);
  end
 
  
  
  for iTraj = 1:length(TrajIndx)
    
    iTraj
    
    [xxIni(1:3,1:3,iTraj), xxdotIni(1:3,1:3,iTraj), RbsIni(1:3,iTraj)] = ComputeCoordAndVeloc(PaQIni(iTraj,1:12));
    [RbsIniSorted(:,iTraj), RbsIniIndx(:,iTraj)] = sort(RbsIni(:,iTraj),'ascend');
    
    [cmIni(1:3,1:2,iTraj), cmdotIni(1:3,1:2,iTraj), xreldistIni(1:3,iTraj), vreldistIni(1:3,iTraj), distIni(iTraj), vdistIni(iTraj), xkinIni(iTraj)] = ComputeTranslatedCoord(xxIni(1:3,1:3,iTraj), xxdotIni(1:3,1:3,iTraj));
    
    [EintIni(iTraj), RbondIni(iTraj), AngMomIni(iTraj)] = ComputeState(RbsIniIndx(1,iTraj), xxIni(1:3,1:3,iTraj), xxdotIni(1:3,1:3,iTraj));
    
    
    [xxFin(1:3,1:3,iTraj), xxdotFin(1:3,1:3,iTraj), RbsFin(1:3,iTraj)] = ComputeCoordAndVeloc(PaQFin(iTraj,1:12));
    [RbsFinSorted(:,iTraj), RbsFinIndx(:,iTraj)] = sort(RbsFin(:,iTraj),'ascend');
    
    [cmFin(1:3,1:2,iTraj), cmdotFin(1:3,1:2,iTraj), xreldistFin(1:3,iTraj), vreldistFin(1:3,iTraj), distFin(iTraj), vdistFin(iTraj), xkinFin(iTraj)] = ComputeTranslatedCoord(xxFin(1:3,1:3,iTraj), xxdotFin(1:3,1:3,iTraj));

    [EintFin(iTraj), RbondFin(iTraj), AngMomFin(iTraj)] = ComputeState(RbsFinIndx(1,iTraj), xxFin(1:3,1:3,iTraj), xxdotFin(1:3,1:3,iTraj));
    
    [riFin(iTraj), roFin(iTraj), jqnFin(iTraj), vqnFin(iTraj), iTypeFin(iTraj)] = AnalyzeFinalState(Leveljqn(iBins,1), AngMomFin(iTraj), EintFin(iTraj), RbondFin(iTraj));
    
  end
  

  figure(iFigure)
  for iTraj = 1:length(TrajIndx)
    plot3(xxIni(1,:,iTraj),xxIni(2,:,iTraj),xxIni(3,:,iTraj),'ro')
    hold on
    plot3([xxIni(1,Pair_to_Atoms(RbsIniIndx(1,iTraj),1),iTraj) xxIni(1,Pair_to_Atoms(RbsIniIndx(1,iTraj),2),iTraj)], [xxIni(2,Pair_to_Atoms(RbsIniIndx(1,iTraj),1),iTraj) xxIni(2,Pair_to_Atoms(RbsIniIndx(1,iTraj),2),iTraj)], [xxIni(3,Pair_to_Atoms(RbsIniIndx(1,iTraj),1),iTraj) xxIni(3,Pair_to_Atoms(RbsIniIndx(1,iTraj),2),iTraj)],'-ok')  
  end
  xlabel('x [a_{0}]')
  ylabel('y [a_{0}]')
  zlabel('z [a_{0}]')
  set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
  set(gcf, 'PaperPositionMode', 'auto');
  if SaveFigs == 1
      FileName = strcat(MoleculesName(1,:),'-Dissociation');
      FilePathNamePng = strcat(FilePath, FileName, '.png');
      FilePathNameFig = strcat(FilePath, FileName, '.fig');
      saveas(gcf,strcat(FilePath, FileName),'pdf')
      saveas(gcf,FilePathNamePng)
      savefig(FilePathNameFig);
  end
  iFigure=iFigure+1;
  
  figure(iFigure)
  cmap = colormap(lines(max(LevToBin)));
  for iTraj = 1:length(TrajIndx)
    plot3(xxFin(1,:,iTraj),xxFin(2,:,iTraj),xxFin(3,:,iTraj),'go')
    hold on
    if iTypeFin(iTraj) < 2
      plot3([xxFin(1,Pair_to_Atoms(RbsFinIndx(1,iTraj),1),iTraj) xxFin(1,Pair_to_Atoms(RbsFinIndx(1,iTraj),2),iTraj)], [xxFin(2,Pair_to_Atoms(RbsFinIndx(1,iTraj),1),iTraj) xxFin(2,Pair_to_Atoms(RbsFinIndx(1,iTraj),2),iTraj)], [xxFin(3,Pair_to_Atoms(RbsFinIndx(1,iTraj),1),iTraj) xxFin(3,Pair_to_Atoms(RbsFinIndx(1,iTraj),2),iTraj)],'o-','Color',cmap(LevToBin(qnToBin(vqnFin(iTraj)+1,jqnFin(iTraj))),1:3))
    else
      plot3([xxFin(1,Pair_to_Atoms(RbsFinIndx(1,iTraj),1),iTraj) xxFin(1,Pair_to_Atoms(RbsFinIndx(1,iTraj),2),iTraj)], [xxFin(2,Pair_to_Atoms(RbsFinIndx(1,iTraj),1),iTraj) xxFin(2,Pair_to_Atoms(RbsFinIndx(1,iTraj),2),iTraj)], [xxFin(3,Pair_to_Atoms(RbsFinIndx(1,iTraj),1),iTraj) xxFin(3,Pair_to_Atoms(RbsFinIndx(1,iTraj),2),iTraj)],'o-k')
    end
  end
  xlabel('x [a_{0}]')
  ylabel('y [a_{0}]')
  zlabel('z [a_{0}]')
  set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
  set(gcf, 'PaperPositionMode', 'auto');
  if SaveFigs == 1
      FileName = strcat(MoleculesName(1,:),'-Dissociation');
      FilePathNamePng = strcat(FilePath, FileName, '.png');
      FilePathNameFig = strcat(FilePath, FileName, '.fig');
      saveas(gcf,strcat(FilePath, FileName),'pdf')
      saveas(gcf,FilePathNamePng)
      savefig(FilePathNameFig);
  end
  iFigure=iFigure+1;
  
  
  figure(iFigure)
  Nr   = 1000;
  rMin = 1.5;
  rMax = 10.d0;
  figure
  cmap = colormap(lines(max(LevToBin)));
  for iTraj = 1:length(TrajIndx)
    if iTypeFin(iTraj) < 2
      rvec = linspace(rMin,rMax,Nr);
      for ir = 1:Nr
        [Ve(ir), dVe] = DiatPot(rvec(ir), jqnFin(iTraj), System);
      end
      plot3(rvec,jqnFin(iTraj)*ones(Nr,1),Ve,'k');
      hold on
      [vib, dVv] = LeRoy(RbondFin(iTraj));
      EintShiftFin(iTraj) = EintFin(iTraj) + vib;
      plot3([riFin(iTraj); roFin(iTraj)], [jqnFin(iTraj); jqnFin(iTraj)], [EintShiftFin(iTraj); EintShiftFin(iTraj)], '-','Color',cmap(LevToBin(qnToBin(vqnFin(iTraj)+1,jqnFin(iTraj))),1:3));
    end
  end
  h=colorbar;
  xlabel('R [a_{0}]')
  ylabel('J')
  zlabel('Ve [Eh]')
  set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
  set(gcf, 'PaperPositionMode', 'auto');
  if SaveFigs == 1
      FileName = strcat(MoleculesName(1,:),'-Dissociation');
      FilePathNamePng = strcat(FilePath, FileName, '.png');
      FilePathNameFig = strcat(FilePath, FileName, '.fig');
      saveas(gcf,strcat(FilePath, FileName),'pdf')
      saveas(gcf,FilePathNamePng)
      savefig(FilePathNameFig);
  end
  iFigure=iFigure+1;
 

  figure(iFigure)
  histogram(xkinIni,20,'Normalization','probability')
  hold on
  histogram(xkinFin,20,'Normalization','probability')
  set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
  set(gcf, 'PaperPositionMode', 'auto');
  xlabel('E_{Kin} [Eh]')
  ylabel('PDF')
  if SaveFigs == 1
      FileName = strcat(MoleculesName(1,:),'-Dissociation');
      FilePathNamePng = strcat(FilePath, FileName, '.png');
      FilePathNameFig = strcat(FilePath, FileName, '.fig');
      saveas(gcf,strcat(FilePath, FileName),'pdf')
      saveas(gcf,FilePathNamePng)
      savefig(FilePathNameFig);
  end
  iFigure=iFigure+1;
  
  
  figure(iFigure)
  histogram(EintIni,50,'Normalization','probability')
  hold on
  histogram(EintFin,50,'Normalization','probability')
  set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
  set(gcf, 'PaperPositionMode', 'auto');
  xlabel('E_{Int} [Eh]')
  ylabel('PDF')
  if SaveFigs == 1
      FileName = strcat(MoleculesName(1,:),'-Dissociation');
      FilePathNamePng = strcat(FilePath, FileName, '.png');
      FilePathNameFig = strcat(FilePath, FileName, '.fig');
      saveas(gcf,strcat(FilePath, FileName),'pdf')
      saveas(gcf,FilePathNamePng)
      savefig(FilePathNameFig);
  end
  iFigure=iFigure+1;
  
  
  figure(iFigure)
  histogram(RbondIni,50,'Normalization','probability')
  hold on
  histogram(RbondFin,50,'Normalization','probability')
  set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
  set(gcf, 'PaperPositionMode', 'auto');
  xlabel('R_{Bond} [a_{0}]')
  ylabel('PDF')
  if SaveFigs == 1
      FileName = strcat(MoleculesName(1,:),'-Dissociation');
      FilePathNamePng = strcat(FilePath, FileName, '.png');
      FilePathNameFig = strcat(FilePath, FileName, '.fig');
      saveas(gcf,strcat(FilePath, FileName),'pdf')
      saveas(gcf,FilePathNamePng)
      savefig(FilePathNameFig);
  end
  iFigure=iFigure+1;
  
  
  figure(iFigure)
  histogram(AngMomIni,50,'Normalization','probability')
  hold on
  histogram(AngMomFin,50,'Normalization','probability')
  set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
  set(gcf, 'PaperPositionMode', 'auto');
  xlabel('Angolar Momentum]')
  ylabel('PDF')
  if SaveFigs == 1
      FileName = strcat(MoleculesName(1,:),'-Dissociation');
      FilePathNamePng = strcat(FilePath, FileName, '.png');
      FilePathNameFig = strcat(FilePath, FileName, '.fig');
      saveas(gcf,strcat(FilePath, FileName),'pdf')
      saveas(gcf,FilePathNamePng)
      savefig(FilePathNameFig);
  end
  iFigure=iFigure+1;
  
end