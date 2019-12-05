% -- MATLAB --
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

function [iFigure] = PlotTrajectory(iFigure, iTint)

  global iTraj tMin tMax StartBin FinalBin iNode iProc AtomColor AtomSize AllMoleculesName PairColor NTint T0_Vec
  
  DimMin    = -30.d0;
  DimMax    =  30.d0;
  iStepInt  = 100;
  DeltaStep = 3; 

  for iBins = StartBin:FinalBin
    
    f = figure;
    fig = gcf;
    screensize = get( groot, 'Screensize' );
    fig.Position=screensize;
    set(gcf,'color','w');
    pause(15)

    filename = strcat('../Test/T_',num2str(T0_Vec(iTint)),'_',num2str(T0_Vec(iTint)),'/Bins_',num2str(iBins),'_0/Node_',num2str(iNode),'/Proc_',num2str(iProc),'/trajectories.out')
    startRow = 2;
    formatSpec = '%8f%7f%14f%14f%14f%17f%17f%17f%17f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    TrajIdxVec = dataArray{:, 1};
    iPES = dataArray{:, 2};
    bMaxVec = dataArray{:, 3};
    bVec = dataArray{:, 4};
    jqnIniVec = dataArray{:, 5};
    vqnIniVec = dataArray{:, 6};
    ArrIniVec = dataArray{:, 7};
    jqnFinVec = dataArray{:, 8};
    vqnFinVec = dataArray{:, 9};
    ArrFinVec = dataArray{:, 10};
    clearvars filename startRow formatSpec fileID dataArray ans;
    for jTraj = 1:length(TrajIdxVec)
      if TrajIdxVec(jTraj) == iTraj
        bMax    = bMaxVec(jTraj);
        b       = bVec(jTraj);
        jqnIni  = round(jqnIniVec(jTraj)-0.5d0);
        vqnIni  = round(vqnIniVec(jTraj)-0.5d0);
        ArrIni  = fix(round(ArrIniVec(jTraj)-0.5d0)/16);
        jqnFin  = round(jqnFinVec(jTraj)-0.5d0);
        vqnFin  = round(vqnFinVec(jTraj)-0.5d0);
        ArrFin  = fix(round(ArrFinVec(jTraj)-0.5d0)/16);
      end
    end


    filename = strcat('../Test/T_',num2str(T0_Vec(iTint)),'_',num2str(T0_Vec(iTint)),'/Bins_',num2str(iBins),'_0/Node_',num2str(iNode),'/Proc_',num2str(iProc),'/XXEvo-',num2str(iTraj),'.out')
    startRow = 2;
    formatSpec = '%20f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    t = dataArray{:, 1};
    clearvars formatSpec fileID dataArray ans;
    formatSpec = '%*20s%18f%18f%18f%18f%18f%18f%18f%18f%18f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    XTemp = [dataArray{1:end-1}];
    clearvars formatSpec fileID dataArray ans;
    formatSpec = '%*182s%18f%18f%18f%18f%18f%18f%18f%18f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    VTemp = [dataArray{1:end-1}];
    clearvars filename startRow formatSpec fileID dataArray ans;

    for iStep = 1:length(t)
      k=1;
      for iAtom=1:3
        for iDir=1:3
            XMatrix(iDir,iAtom,iStep) = XTemp(iStep,k);
            VMatrix(iDir,iAtom,iStep) = VTemp(iStep,k);
          k=k+1;
        end
      end
    end


    RMatrix  = zeros(length(t),3);
    dVMatrix = zeros(length(t),3);
    filename = strcat('../Test/T_',num2str(T0_Vec(iTint)),'_',num2str(T0_Vec(iTint)),'/Bins_',num2str(iBins),'_0/Node_',num2str(iNode),'/Proc_',num2str(iProc),'/PESEvo-',num2str(iTraj),'.out');
    startRow = 2;
    formatSpec = '%20f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    timePES = dataArray{:, 1};
    clearvars formatSpec dataArray ans;
    formatSpec = '%*74s%18f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    VTemp = dataArray{:, 1};
    clearvars formatSpec fileID dataArray ans;
    fileID = fopen(filename,'r');
    formatSpec = '%*20s%18f%18f%18f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    RMatrixTemp = [dataArray{1:end-1}];
    clearvars formatSpec fileID dataArray ans;
    formatSpec = '%*92s%18f%18f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    dVMatrixTemp = [dataArray{1:end-1}];
    clearvars filename startRow formatSpec fileID dataArray ans;
    MaxStepPES = length(timePES);
    RMatrix(1,:)  = RMatrixTemp(1,:);
    V(1)          = VTemp(1);
    dVMatrix(1,:) = dVMatrixTemp(1,:);
    iStepPES  = 1;
    for iStep = 2:length(t)
       while timePES(iStepPES) < t(iStep) && iStepPES < MaxStepPES
        iStepPES = iStepPES + 1;
       end
       iStepPES = iStepPES-1;
       RMatrix(iStep,:)  = RMatrixTemp(iStepPES,:);
       V(iStep)          = VTemp(iStepPES);
       dVMatrix(iStep,:) = dVMatrixTemp(iStepPES,:);
    end


    x = linspace(DimMin,DimMax,2);
    y = linspace(DimMin,DimMax,2);
    [X,Y] = meshgrid(x,y);
    Z     = 0.d0*X;

    p = uipanel('Parent',f,'BorderType','none','BackgroundColor','white'); 
%     if ArrFin ~= 0 
%       p.Title = strcat(AllMoleculesName(ArrIni,:),': (v,J) = (', num2str(vqnIni), ',', num2str(jqnIni), ') -> ',AllMoleculesName(ArrFin,:), ': (v,J) = (', num2str(vqnFin), ',', num2str(jqnFin),')'); 
%     else
%       p.Title = strcat(AllMoleculesName(ArrIni,:),': (v,J) = (', num2str(vqnIni), ',', num2str(jqnIni), ') -> Dissociation'); 
%     end
%     p.TitlePosition = 'centertop'; 
%     p.FontSize = 20;
    set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter','latex','color','w');
    %p.FontWeight = 'bold';
    iStep=1
    while t(iStep) < tMin
      iStep = iStep + 1;
    end
    iStepIni = iStep
    iStep=iStepIni
    while t(iStep) < tMax && iStep < length(t)
      iStep = iStep + 1;
    end
    iStepFin = min(iStep,length(t))
    for iStep = iStepIni:DeltaStep:iStepFin


      subplot(2,3,[1 2 4 5],'Parent',p)
      s=surf(X,Y,Z,'FaceColor',[0.94,0.94,0.94]);
      alpha(s,'z')
      hold on
      scatter3(XMatrix(1,:,iStep),XMatrix(2,:,iStep),XMatrix(3,:,iStep),AtomSize,AtomColor,'filled')
      q1 = quiver3(XMatrix(1,:,iStep),XMatrix(2,:,iStep),XMatrix(3,:,iStep),VMatrix(1,:,iStep),VMatrix(2,:,iStep),VMatrix(3,:,iStep),0.3);
      % c = q1.Color;
      % q1.Color = 'k';
      % set(gca,'Color','c')
      hold off
      grid off
      set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex','color','w');
      set(gcf, 'PaperPositionMode', 'auto','color','w');
      %title(['Time = ', num2str(t(iStep)), ' a.u.']);
      view(125,10);
      xlim([DimMin, DimMax]);
      ylim([DimMin, DimMax]);
      zlim([DimMin, DimMax]);
      xlabel(['X [bohr]']);
      ylabel(['Y [bohr]']);
      zlabel(['Z [bohr]']);


      subplot(2,3,3,'Parent',p)
      plot(t(max(iStep-iStepInt,1):iStep,1),RMatrix(max(iStep-iStepInt,1):iStep,1),'Color',PairColor(3,:))
      hold on 
      plot(t(max(iStep-iStepInt,1):iStep,1),RMatrix(max(iStep-iStepInt,1):iStep,2),'Color',PairColor(2,:))
      plot(t(max(iStep-iStepInt,1):iStep,1),RMatrix(max(iStep-iStepInt,1):iStep,3),'Color',PairColor(1,:))
      %hold off
      set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex','color','w');
      set(gcf, 'PaperPositionMode', 'auto','color','w');
      title(['Time = ', num2str(t(iStep)), ' a.u.']);
      %xlim([DimMin, DimMax]);
      %ylim([DimMin, DimMax]);
      xlabel(['time [a.u.]']);
      ylabel(['R [bohr]']);


      subplot(2,3,6,'Parent',p)
      plot(t(max(iStep-iStepInt,1):iStep),dVMatrix(max(iStep-iStepInt,1):iStep,1),'Color',PairColor(3,:))
      hold on 
      plot(t(max(iStep-iStepInt,1):iStep),dVMatrix(max(iStep-iStepInt,1):iStep,2),'Color',PairColor(2,:))
      plot(t(max(iStep-iStepInt,1):iStep),dVMatrix(max(iStep-iStepInt,1):iStep,3),'Color',PairColor(1,:))
      % hold off
      set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex','color','w');
      set(gcf, 'PaperPositionMode', 'auto','color','w');
      %title(['Time = ', num2str(t(iStep)), ' a.u.']);
      %xlim([DimMin, DimMax]);
      %ylim([DimMin, DimMax]);
      xlabel(['time [a.u.]']);
      ylabel(['dV/dR [Eh/bohr]']);
      
      %pause
      pause(0.0001)

    end

    pause

  end
      
end
