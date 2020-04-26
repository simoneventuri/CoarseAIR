%% The Function reads the rates (K(i,j)) and computes the Processes (Dissociation, Exchange 1,2 and 3) overall rates (K(i,Process))
%
%  Input Arguments:  - iT:            Index for the current Translational Temperature
%                    - iBinsStart:    (Optional) Initial Level/Bin to read. If not given, it is assumed equal to be equal to 1
%                    - iBinsEnd:      (Optional) Final Level/Bin to read. If not given, it is assumed equal to be equal to NBins
%
%  Input Global Var: - T0_Vec:        Vector of Translational Temperatures (e.g.: [10000])
%                    - RatesPath:     The path to the output folder (e.g.: ../Test/N3/N2/Rates)
%                    - MoleculesName: A vector of strings containing the name of the molecules present in the system
%                    - NBins          Nb of Levels/Bins
%

%function [RatesSigma, RatesMatrix, RatesSigmaMatrix, DissRates, DissRatesSigma, ProcessesRates, ProcessesRatesSigma] = ReadRates(iT, RatesSigma, RatesMatrix, RatesSigmaMatrix, DissRates, DissRatesSigma, ProcessesRates, ProcessesRatesSigma, iBinsStart, iBinsEnd)    
function [RatesMatrix, DissRates, ProcessesRates] = ReadRatesSingleBin(iT, RatesMatrix, DissRates, ProcessesRates, iBinsStart, iBinsEnd)    

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

  global NBins T0_Vec RatesPath ProcToLevIP   
  

  for iBins = iBinsStart:iBinsEnd      

    filenameBinsRates = strcat(RatesPath,'/T_',num2str(T0_Vec(iT)),'_',num2str(T0_Vec(iT)),'/Bin',num2str(iBins),'.dat')

    if exist(filenameBinsRates, 'file') == 2

      startRow = 6;
      formatSpec = '%*24s%16f%20f%[^\n\r]';
      fileID = fopen(filenameBinsRates,'r');
      dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

  %       formatSpec = '%*24s%16f%20f%20f%[^\n\r]';
  %       fileID = fopen(filenameBinsRates,'r');
  %       dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

  %       formatSpec = '%*19s%21f%20f%[^\n\r]';
  %       fileID = fopen(filenameBinsRates,'r');
  %       dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

      fclose(fileID);
      Process    = dataArray{:, 1};
      RateTemp   = dataArray{:, 2};
      SigmaTemp  = zeros(length(RateTemp),1);
      %SigmaTemp  = dataArray{:, 3};

      if length(Process) > 0
        if Process(1) ~= 1
          for ii=1:length(Process)
            RatesMatrix(iBins,ProcToLevIP(Process(ii),1),ProcToLevIP(Process(ii),2),iT)      = RateTemp(ii);
            %RatesSigmaMatrix(iBins,ProcToLevIP(Process(ii),1),ProcToLevIP(Process(ii),2),iT) = SigmaTemp(ii);
          end
        else
          ProcessesRates(iBins,1,iT)      = RateTemp(1);
          ProcessesRatesSigma(iBins,1,iT) = SigmaTemp(1);  
          DissRates(iBins,iT)             = RateTemp(1);  
          DissRatesSigma(iBins,iT)        = SigmaTemp(1);  
          for ii=2:length(Process)
            RatesMatrix(iBins,ProcToLevIP(Process(ii),1),ProcToLevIP(Process(ii),2),iT)      = RateTemp(ii); 
            %RatesSigmaMatrix(iBins,ProcToLevIP(Process(ii),1),ProcToLevIP(Process(ii),2),iT) = SigmaTemp(ii); 
          end
        end 
        for iP=1:3
          ProcessesRates(iBins,iP+1,iT)      = sum(RatesMatrix(iBins,:,iP,iT));
          %ProcessesRatesSigma(iBins,iP+1,iT) = sqrt(sum(RatesSigmaMatrix(iBins,:,iP,iT).^2));
        end
      end
      clearvars RateTemp Process filename startRow formatSpec fileID dataArray ans;

    end

  end
      
end