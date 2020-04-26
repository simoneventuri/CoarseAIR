%% The Function plots the Molecules' Diatomic Potentials and their Energy Levels
%
%  Input Arguments:  - ...
%
%  Input Global Var: - PlotDiatPotFlg:          Flag 1/2; if =1, the Diat. Potential is computed by means of MATLAB Functions; if =2, it is obtained from Fortran Code.
%                    - PlotPairUser:            Vectors of Flags 0/1 containing the pairs to plot
%                    - OutputFileName:          Name of the File that contains the Output Quantity
%                    - ForceFlg:                Flag 0/1; if =0, Output Quantity is a Diat. Potential; if =1, Output Quantity is Diat. Force.
%

function [iFigure] = PlotDiatPot(iFigure, NLevels, Leveljqn, Levelvqn, LevelEh, LevelEeV, rMin, rMax, VMin, VMax, rIn, rOut, LevToBin, DeltaEintDiss, LevelEeVVib0, LevelEeVRot)

  
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
  
  global PlotDiatPotFlg PlotPairUser OutputFileName ForceFlg ParaViewFlg
  global PathToOutput PlotPair Pair_To_Molecule MoleculesName AtomsName NMolecules
  global Pair_to_Atoms
  global AtomMass
  global SaveFigs FigDirPath AxisFontSz AxisFontNm LegendFontSz AxisLabelSz AxisLabelNm LegendFontNm XLimPlot YLimPlot PlotPairUservqnColor 
  global RCVec BCVec GCVec KCVec OCVec PCVec WCVec JCVec YCVec CCVec MCVec
  global B_To_Ang Kcm_To_Hart KcmAng_To_HartB Hart_To_eV EhToeV EnergyUM 

  if PlotDiatPotFlg == 1

    rMinb = 1.55d0;
    rMaxx = 10.d0;
    
%     iMol                = Pair_To_Molecule(1);%%% for assinging the same                                  
%     rMin = min(rIn(:,iMol))                   %%%   rmin and rmax as the energy levels
%     rMax = max(rOut(:,iMol))                  %%%
    Nr   = 250;
    jqn  = 0;

    yConv = 1;
    if EnergyUM == 'eV'
      yConv = EhToeV;
    end
    yConv
    

    for iMol = 1:NMolecules 

      for iLevels = 1:NLevels(iMol)
        VMaxJ(Leveljqn(iLevels,iMol)+1,iMol) = VMax(iLevels,iMol) ./ EhToeV;
      end


      if ParaViewFlg == 1
       
        FolderName         = strcat('./', MoleculesName(iMol,:), '_DiatPot')
        [status,msg,msgID] = mkdir(FolderName)

        fileName = strcat(FolderName,'/CentrifugalBarrier.csv');
        fileID   = fopen(fileName,'w');
        fprintf(fileID,'Variables = "rMax", "JMax", "EMax"\n');
        jj       = 0;
        iLevel   = 1;
        while iLevel <= NLevels(iMol)
          if Leveljqn(iLevel,iMol) == jj
            if rMax(iLevel,iMol) <= 30.d0
              fprintf(fileID,'%15.6e, %15.6e, %15.6e\n',rMax(iLevel,iMol), jj, VMax(iLevel,iMol));  
            end
            jj = jj + 1;
          end
          iLevel = iLevel + 1;
        end
        fclose(fileID);


%         fileName1 = strcat(FolderName,'/SurfacePoints.csv');
%         fileID1   = fopen(fileName1,'w');
%         fprintf(fileID1,'Variables = "r", "J", "E", "dEdr"\n');
%         jj       = 0;
%         iLevel   = 1;
%         while iLevel <= NLevels(iMol)
%           if Leveljqn(iLevel,iMol) == jj
%             MaxLevel = iLevel;
%             for jLevel = 1:NLevels(iMol)
%               if LevelEeV(jLevel,iMol) > LevelEeV(MaxLevel,iMol) && Leveljqn(jLevel,iMol) == jj
%                 MaxLevel = jLevel;
%               end
%             end
%             rvec = linspace(rMin(iLevel,iMol),min(rMax(iLevel,iMol)+5.d0,30.d0),Nr);
%             %rvec = linspace(rIn(MaxLevel,iMol),min(rMax(iLevel,iMol)+5.d0,30.d0),Nr);
%             for ir=1:Nr
%               [Ve, dVe] = DiatPot(rvec(ir), jj, iMol);
%               %Ve = Ve;%- VMaxJ(iJqn+1,iMol);
%               fprintf(fileID1,'%15.6e, %15.6e, %15.6e, %15.6e\n',rvec(ir), jj, Ve.*yConv, dVe.*yConv);
%             end
%             jj = jj + 1;
%           end
%           iLevel = iLevel + 1;
%         end
%         fclose(fileID1);
        
        FileName = strcat(FolderName,'/InelLevels.csv');
        fileID = fopen(FileName,'w');
        fprintf(fileID,'id,v,J,EeV,rIn,rOut,EeVVib,EeVRot,dCentBarr\n');
        for i = 1:NLevels(iMol)
          fprintf(fileID,'%i,%i,%i,%e,%e,%e,%e,%e,%e\n', i, Levelvqn(i,iMol), Leveljqn(i,iMol), LevelEeV(i,iMol), rIn(i,iMol), rOut(i,iMol), LevelEeVVib0(i,iMol), LevelEeVRot(i,iMol), DeltaEintDiss(i,iMol));
        end
        fclose(fileID);
        
        
      end
        
    end


    for iP = 1:3
      if PlotPair(iP) == 1
        MoleculeName(1,1:2) = MoleculesName(Pair_To_Molecule(iP),1:2);
        iMol                = Pair_To_Molecule(iP);

%         figure(iFigure)
%   %     %   for jqn = 1:2:200
%   %     %     rvec = linspace(rMin,rMax,Nr);
%   %     %     for ir = 1:Nr
%   %     %       [Vv(ir)]  = LeRoy(rvec(ir));
%   %     %       [Vc(ir)]  = CentPot(rvec(ir), jqn);
%   %     %     end
%   %     %     Ve =  Vv + Vc;
%   %     %     plot3(rvec,jqn*ones(Nr,1),Ve,'k');
%   %     %     hold on
%   %     %   end
%   % 
%          clear cmap
%          cmap = colormap(lines(max(LevToBin(:,iMol))+1));
%         for i = 1:NLevels(iMol)
%           if mod(Leveljqn(i,iMol),2)==0
%             plot3([rIn(i,iMol); rOut(i,iMol)], [Leveljqn(i,iMol); Leveljqn(i,iMol)], [LevelEeV(i,iMol); LevelEeV(i,iMol)], '-','Color',cmap(LevToBin(i),1:3));
%             hold on
%           end
%         end
%         set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
%         set(gcf, 'PaperPositionMode', 'auto');
%         if SaveFigs == 1
%           FileName = strcat(Molecules(1,:),'-',Pair_Name(iP,:),'-',TempChar(iT,5:end));
%           FilePathNamePng = strcat(FilePath, FileName, '.png');
%           FilePathNameFig = strcat(FilePath, FileName, '.fig');
%           saveas(gcf,strcat(FilePath, FileName),'pdf')
%           saveas(gcf,FilePathNamePng)
%           savefig(FilePathNameFig);
%         end
%         iFigure=iFigure+1;

%         cmap = colormap(lines(max(LevToBin(:,iMol))+1));
%         iTot = 1;
%         for i = 1:NLevels(iMol)
%           %if mod(Leveljqn(i,iMol),2)==0
%           if rIn(i,iMol) >= 20
%             i;
%             pause
%           end
%             X(iTot,1:2)    = [rIn(i,iMol); rOut(i,iMol)]; 
%             Y(iTot,1:2)    = [Leveljqn(i,iMol); Leveljqn(i,iMol)];
%             Z(iTot,1:2)    = [LevelEeV(i,iMol); LevelEeV(i,iMol)];%.*yConv; %%% already in Ev
%             ColorVec(iTot) = DeltaEintDiss(i,1);%LevToBin(i,iMol);
%             iTot           = iTot + 1;
%           %end
%         end
%         FolderStr = strcat('./', MoleculeName(1,1:2), '_DiatPot/');
%         [status, msg, msgID] = mkdir(FolderStr);
%         Write3DLinesToVTK(FolderStr, 'LevToBin', X, Y, Z, ColorVec, cmap)
%         clear X Y Z ColorVec

      end 
    end
    
  elseif PlotDiatPotFlg == 2
    
    for iP = 1:3
      
      if PlotPairUser(iP) == 1 
    
        OutputFile = strcat(PathToOutput,'/PlotPES/',OutputFileName)
        delimiter = ',';
        startRow = 2;
        formatSpec = '%f%f%f%[^\n\r]';
        fileID = fopen(OutputFile,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
        fclose(fileID);
        R1(:) = dataArray{:, 1};
        R2(:) = dataArray{:, 2};
        V(:) = dataArray{:, 3};
        clearvars filename delimiter startRow formatSpec fileID dataArray ans;

        if ForceFlg == 0

          for i=1:2:length(R1(:))
            RTemp(ceil(i/2))=R1(i); 
            VTemp(ceil(i/2))=V(i); 
          end

          h = RTemp(2) - RTemp(1);
          dV = zeros(length(RTemp),1);
          Vint = zeros(length(RTemp),1);

          for i=3:length(RTemp)-2
            dV(i) = ( -VTemp(i+2) + 8.d0*VTemp(i+1) - 8.d0*VTemp(i-1) + VTemp(i-2) ) / (12.d0*h);
            Vint(i)=trapz(RTemp(1:i),dV(1:i)) + VTemp(3);
          end

        end

        Exts_x      = [0.8, 5];
        iP
        str_x       = ['$r(',AtomsName(Pair_to_Atoms(iP,1)),'_A - ',AtomsName(Pair_to_Atoms(iP,2)),'_B) [Bohr]$'];
        if ForceFlg == 0
          str_y       = '$\Delta E [Hartree]$';
          Exts_y      = [-1, 8];
        else
          str_y       = '$Force [Hartree/Bohr]$';
          Exts_y      = [-15, 2];
        end
        LineStr     = '-';
        LineSize    = 3;
        LegendFlg   = 0;
        LegendStr   = {''};
        FolderStr   = '/DiatPot/';
        if ForceFlg == 0
          FileStr   = 'DiatPot';
          LineColor   = GCVec;
          PlotFigure(iFigure, RTemp(3:end-2), dV(3:end-2), Exts_x, Exts_y, str_x, str_y, LineStr, LineColor, LineSize, LegendFlg, LegendStr, FolderStr, FileStr );
        else
          FileStr   = 'Force';
        end
        LineColor   = BCVec;
        PlotFigure(iFigure, R1(:), V(:)-V(end), Exts_x, Exts_y, str_x, str_y, LineStr, LineColor, LineSize, LegendFlg, LegendStr, FolderStr, FileStr );
        iFigure     = iFigure + 1;
        
      end
      
    end
    
  end

end