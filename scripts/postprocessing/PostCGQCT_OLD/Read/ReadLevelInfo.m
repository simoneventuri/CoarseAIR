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

function [NLevels, Levelvqn, Leveljqn, LevelEh, LevelEeV, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelg, LevelToBin, LevelQ, vToLevel, DeltaEintDiss, DeltaEintDepth, rIn, rOut, rMin, rMax, VMin, VMax, Tau, Egam, dVIn, ddVIn, dVOut, ddVOut, dVdJIn, dVdJOut, rEquival] = ReadLevelInfo(iT, EeV);
      
    
  %% (METHOD DEPENDENT)

  global T0_Vec System SystemPath 
  global NBins BinnedMolName NBinnedMol KeV KinMthd EhToeV MoleculesName AtomMass Pair_to_Atoms ParaViewFlg
  
  
  RecomputeMinMaxFlg               = false%true
  RecomputeInOutFlg                = false%true
  ComputerEquivalFlg               = false%true
  ComputeDerivativesAtTournPntsFlg = false
  ReadLevelsInfoFlg                = true
  WriteFileFlg                     = false
  

  RejectedLevels = [];
  for iBinnedMol=1:NBinnedMol

    if sum(KinMthd(iBinnedMol,:) == 'CGM')==3 || sum(KinMthd(iBinnedMol,:) == 'VIB')==3
      filename = strcat(SystemPath,'/',BinnedMolName(iBinnedMol,:),'/',BinnedMolName(iBinnedMol,:),'_',num2str(NBins(iBinnedMol)),'/qnsEnBin.dat')
    elseif sum(KinMthd(iBinnedMol,:) =='STS')==3
      filename = strcat(SystemPath,'/',BinnedMolName(iBinnedMol,:),'/',BinnedMolName(iBinnedMol,:),'_',num2str(NBins(iBinnedMol)),'/qnsEnBin.dat')
    end
    startRow = 2;
%     formatSpec = '%8f%7f%7f%15f%15f%f%[^\n\r]';
%     fileID = fopen(filename,'r');
%     dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%     formatSpec = '%8f%7f%7f%15f%15f%f%[^\n\r]';
%     fileID = fopen(filename,'r');
%     dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');        
    formatSpec = '%8f%7f%7f%15f%15f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    NLevels(iBinnedMol)=size(dataArray{:,1},1)
    Levelvqn(1:NLevels(iBinnedMol),iBinnedMol)   = dataArray{:, 2};
    Leveljqn(1:NLevels(iBinnedMol),iBinnedMol)   = dataArray{:, 3};
    LevelEeV(1:NLevels(iBinnedMol),iBinnedMol)   = dataArray{:, 4};
    Levelg(1:NLevels(iBinnedMol),iBinnedMol)     = dataArray{:, 5};
    LevelToBin(1:NLevels(iBinnedMol),iBinnedMol) = dataArray{:, 6};
    clearvars filename startRow formatSpec fileID dataArray ans;
    
    
    
    %% Read Levels Original Data File
    filename = strcat(SystemPath,'/',BinnedMolName(iBinnedMol,:),'/levels_cut.inp')
    startRow = 16;
    formatSpec = '%6f%5f%15f%15f%15f%15f%15f%15f%15f%15f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    vqnVec  = dataArray{:, 1};
    jqnVec  = dataArray{:, 2};
    EintVec = dataArray{:, 3};
    Egam(1:NLevels(iBinnedMol),iBinnedMol) = max(dataArray{:, 4},1.d-300);
    rMin(1:NLevels(iBinnedMol),iBinnedMol) = dataArray{:, 5};
    rMax(1:NLevels(iBinnedMol),iBinnedMol) = dataArray{:, 6};
    VMin(1:NLevels(iBinnedMol),iBinnedMol) = dataArray{:, 7}.*EhToeV;
    VMax(1:NLevels(iBinnedMol),iBinnedMol) = dataArray{:, 8}.*EhToeV;
    Tau(1:NLevels(iBinnedMol),iBinnedMol)  = dataArray{:, 9};
    rIn(1:NLevels(iBinnedMol),iBinnedMol)  = dataArray{:, 10};
    rOut(1:NLevels(iBinnedMol),iBinnedMol) = dataArray{:, 11};
    clearvars filename startRow formatSpec fileID dataArray ans;

    LevelEeV(1:NLevels(iBinnedMol),iBinnedMol)   = LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) - VMax(1,iBinnedMol);
    LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol)  = LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) - min(LevelEeV(1:NLevels(iBinnedMol),iBinnedMol));
    LevelEh(1:NLevels(iBinnedMol),iBinnedMol)    = LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) ./ EhToeV;
    VMax(1:NLevels(iBinnedMol),iBinnedMol)       = VMax(1:NLevels(iBinnedMol),iBinnedMol) - VMax(1,iBinnedMol);
    VMin(1:NLevels(iBinnedMol),iBinnedMol)       = VMin(1:NLevels(iBinnedMol),iBinnedMol) - VMax(1,iBinnedMol);

    FileName = strcat('./',BinnedMolName(iBinnedMol,:),'_LevelsInfo_FromDataFile.csv');
    fileID   = fopen(FileName,'w');
    fprintf(fileID,'#J,E,rIn,rOut,rMin,rMax,VMin,VMax\n');
    for iLevels = 1:length(LevelEeV)
      fprintf(fileID,'%i,%e,%e,%e,%e,%e,%e,%e\n', Leveljqn(iLevels), LevelEeV(iLevels,iBinnedMol), rIn(iLevels,iBinnedMol), rOut(iLevels,iBinnedMol), rMin(iLevels,iBinnedMol), rMax(iLevels,iBinnedMol), VMin(iLevels,iBinnedMol), VMax(iLevels,iBinnedMol));
    end
    fclose(fileID);


    
    %% (Re-)Computing Properties
    
    if RecomputeMinMaxFlg
      rMinOld = 1.5;
      rMaxOld = 50;
      [rMin0, VMin0, rMax0, VMax0, VInfTemp] = FindMinMaxDiatPot(0, iBinnedMol, rMinOld, rMaxOld);
      for jqn = 1:max(Leveljqn(:,iBinnedMol))
        jqn
        [rMinTemp(jqn), VMinTemp(jqn), rMaxTemp(jqn), VMaxTemp(jqn), VInfTemp] = FindMinMaxDiatPot(jqn, iBinnedMol, rMinOld, rMaxOld);
        rMinOld = rMinTemp(jqn);
        rMaxOld = rMaxTemp(jqn);
        %[rMinTemp(jqn), VMinTemp(jqn), rMaxTemp(jqn), VMaxTemp(jqn)]
      end
    end
    
    % Opening File for Writing Final Properties 
    if WriteFileFlg
      FileName = strcat('./',BinnedMolName(iBinnedMol,:),'_LevelsInfo_Final.csv');
      fileID   = fopen(FileName,'w');
      fprintf(fileID,'Variables = "J", "v", "rIn", "rOut", "rEquival", "E", "DeltaEDiss", "ERot", "EVib", "g", "dVIn", "dVOut", "rMin", "VMin", "rMax", "VMax", "Tau", "Egam"\n');
    end
    
    jLevels = 0;
    for iLevels = 1:NLevels(iBinnedMol)
      %%iLevels
      
      
      %% Compute rMin & rMax, VMin & VMax
      if RecomputeMinMaxFlg
        if Leveljqn(iLevels,iBinnedMol) ~= 0
          rMin(iLevels,iBinnedMol) = rMinTemp(Leveljqn(iLevels,iBinnedMol)); 
          rMax(iLevels,iBinnedMol) = rMaxTemp(Leveljqn(iLevels,iBinnedMol));
          VMin(iLevels,iBinnedMol) = VMinTemp(Leveljqn(iLevels,iBinnedMol));%.*EhToeV; 
          VMax(iLevels,iBinnedMol) = VMaxTemp(Leveljqn(iLevels,iBinnedMol));%.*EhToeV;
        else
          rMin(iLevels,iBinnedMol) = rMin0; 
          rMax(iLevels,iBinnedMol) = rMax0; 
          VMin(iLevels,iBinnedMol) = VMin0;%.*EhToeV; 
          VMax(iLevels,iBinnedMol) = 0.d0;%VMax0.*EhToeV;
        end
      end
      
      
      %% Compute rIn & rOut
      if RecomputeInOutFlg
        if (LevelEeV(iLevels,iBinnedMol) <= VMax(iLevels,iBinnedMol))
          [rIn(iLevels,iBinnedMol), rOut(iLevels,iBinnedMol), ErrorFlg] = TurningPoints( [1.50d0, rMin(iLevels,iBinnedMol)], [rMin(iLevels,iBinnedMol), rMax(iLevels,iBinnedMol)], LevelEeV(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), iBinnedMol);
        else
          jLevels=jLevels+1;
          RejectedLevels(jLevels,iBinnedMol) = iLevels;
          fprintf('%i-th Molecule, %i-th Level has Energy (%e eV) Larger than the Dissociation Barrier (%e eV)!\n', iBinnedMol, iLevels, LevelEeV(iLevels,iBinnedMol), VMax(iLevels,iBinnedMol))
        end
        % Check V @ rIn and rOut
        %ErrorRx = abs(LevelEeV(iLevels,iBinnedMol) - DiatPot(rIn(iLevels,iBinnedMol),  Leveljqn(iLevels,iBinnedMol), 1));
        %ErrorLx = abs(LevelEeV(iLevels,iBinnedMol) - DiatPot(rOut(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), 1));
      end
      
      
      %% Compute DeltaEintDiss (Distance from Centrifugal Barrier)
      DeltaEintDiss(1:NLevels(iBinnedMol),iBinnedMol)  = LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) - VMax(1:NLevels(iBinnedMol),iBinnedMol); 
      DeltaEintDepth(1:NLevels(iBinnedMol),iBinnedMol) = LevelEeV(1:NLevels(iBinnedMol),iBinnedMol) - VMin(1:NLevels(iBinnedMol),iBinnedMol); 
    
      
      %% Compute rEquival (Equivalent Radius)
      if ComputerEquivalFlg
        NHalfPeriods = 2;
        NPoints      = 100000;
        f      = @(t,y) [y(2); -DiatPotdV(y(1),Leveljqn(iLevels,iBinnedMol),1) ./ 27.2113839712790 ./ (AtomMass(1)*2) ];
        tSpan  = [0.0, Tau(iLevels,iBinnedMol)*NHalfPeriods];
        xInit  = [rIn(iLevels,iBinnedMol), 0.00];
        Opts   = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);
        ODESol = ode45(f, tSpan, xInit, Opts);
        t      = linspace(0, Tau(iLevels,iBinnedMol)*NHalfPeriods, NPoints);
        xx     = deval(ODESol, t);
        rEquival(iLevels,iBinnedMol) = trapz(t', xx(1,:)) ./ (Tau(iLevels,iBinnedMol)*NHalfPeriods);
        %figure
        %plot(ODESol.x, ODESol.y(1,:))
        %hold on 
        %plot(t,yy)
      else
        rEquival(iLevels,iBinnedMol) = 0.d0;
      end
      
      
      %% Compute Derivatives at Tourning Points (Equivalent Radius)
      if ComputeDerivativesAtTournPntsFlg
        [dVIn(iLevels,iBinnedMol), ddVIn(iLevels,iBinnedMol), dVdJIn(iLevels,iBinnedMol)]    = FindDerivativeDiatPot(rIn(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), iBinnedMol, Leveljqn);
        [dVOut(iLevels,iBinnedMol), ddVOut(iLevels,iBinnedMol), dVdJOut(iLevels,iBinnedMol)] = FindDerivativeDiatPot(rOut(iLevels,iBinnedMol), Leveljqn(iLevels,iBinnedMol), iBinnedMol, Leveljqn);
      else
        dVIn(iLevels,iBinnedMol)    = 0.d0;
        ddVIn(iLevels,iBinnedMol)   = 0.d0;
        dVOut(iLevels,iBinnedMol)   = 0.d0;
        ddVOut(iLevels,iBinnedMol)  = 0.d0;
        dVdJIn(iLevels,iBinnedMol)  = 0.d0;
        dVdJOut(iLevels,iBinnedMol) = 0.d0;
      end
      
      
      %% Splotting Energy Contributions
      if (Leveljqn(iLevels,iBinnedMol) == 0) 
        vToLevel(Levelvqn(iLevels,iBinnedMol)+1,iBinnedMol) = iLevels; 
        vEeVVib(Levelvqn(iLevels,iBinnedMol)+1,iBinnedMol)  = LevelEeV(iLevels,iBinnedMol); 
        vEeVVib0(Levelvqn(iLevels,iBinnedMol)+1,iBinnedMol) = LevelEeV0(iLevels,iBinnedMol);         
      end
      LevelEeVVib0(iLevels,iBinnedMol) = LevelEeV0(vToLevel(Levelvqn(iLevels,iBinnedMol)+1,iBinnedMol),iBinnedMol);
      LevelEeVRot(iLevels,iBinnedMol)  = LevelEeV0(iLevels,iBinnedMol) - LevelEeVVib0(iLevels,iBinnedMol);  
      
      
      %% Write Levels Info
      if WriteFileFlg
        fprintf(fileID,'%i,%i,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n',      Leveljqn(iLevels,iBinnedMol), ...
                                                                                       Levelvqn(iLevels,iBinnedMol), ...    
                                                                                       rIn(iLevels,iBinnedMol), ...
                                                                                       rOut(iLevels,iBinnedMol),  ...
                                                                                       rEquival(iLevels,iBinnedMol), ...
                                                                                       LevelEeV(iLevels,iBinnedMol), ...
                                                                                       DeltaEintDiss(iLevels,iBinnedMol), ...
                                                                                       LevelEeVRot(iLevels,iBinnedMol), ...
                                                                                       LevelEeVVib0(iLevels,iBinnedMol), ...
                                                                                       Levelg(iLevels,iBinnedMol), ...
                                                                                       dVIn(iLevels,iBinnedMol), ...
                                                                                       dVOut(iLevels,iBinnedMol), ...
                                                                                       rMin(iLevels,iBinnedMol), ...
                                                                                       VMin(iLevels,iBinnedMol), ...
                                                                                       rMax(iLevels,iBinnedMol), ...
                                                                                       VMax(iLevels,iBinnedMol), ...
                                                                                       Tau(iLevels,iBinnedMol), ...
                                                                                       Egam(iLevels,iBinnedMol) ...
                                                                                       );
      end
      
      
    end
    if WriteFileFlg
      fclose(fileID);
    end
    fprintf('iMol = %i, List of Rejected Levels: \n', iBinnedMol)
    RejectedLevels
    
    
    %% Read Levels Info
    if ReadLevelsInfoFlg
      filename = strcat('./',BinnedMolName(iBinnedMol,:),'_LevelsInfo_Final.csv')
      delimiter = ',';
      startRow = 2;
      formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
      fileID = fopen(filename,'r');
      dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
      fclose(fileID);
      Leveljqn   = dataArray{:, 1};
      Levelvqn   = dataArray{:, 2};
      rIn        = dataArray{:, 3};
      rOut       = dataArray{:, 4};
      rEquival   = dataArray{:, 5};
      LevelEeV   = dataArray{:, 6};
      DeltaEDiss = dataArray{:, 7};
      ERot       = dataArray{:, 8};
      EVib       = dataArray{:, 9};
      Levelg     = dataArray{:, 10};
      dVIn       = dataArray{:, 11};
      dVOut      = dataArray{:, 12};
      rMin       = dataArray{:, 13};
      VMin       = dataArray{:, 14};
      rMax       = dataArray{:, 15};
      VMax       = dataArray{:, 16};
      Tau        = dataArray{:, 17};
      Egam       = dataArray{:, 18};
      clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    end
    
    
    %% Computing Partition Function
    if sum(KinMthd(iBinnedMol,:) == 'CGM')==3 || sum(KinMthd(iBinnedMol,:) == 'VIB')==3
      LevelQ(1:NLevels(iBinnedMol),iBinnedMol) = Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* exp( - (LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol)) ./ (KeV*T0_Vec(iT)) );
    elseif sum(KinMthd(iBinnedMol,:) =='STS')==3
      %LevelQ(1:NLevels(iBinnedMol),iBinnedMol) = Levelg(1:NLevels(iBinnedMol),iBinnedMol);
      LevelQ(1:NLevels(iBinnedMol),iBinnedMol) = Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* exp( - (LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol)) ./ (KeV*T0_Vec(iT)) );
    end 
    LevelQ = LevelQ ./ sum(LevelQ,1);

    
  end
    
end
