%% The Function reads the rates (K(i,j)) and computes the Processes (Dissociation, Exchange 1,2 and 3) overall rates (K(i,Process))
%
%  Input Arguments:  - iT:               Index for the current Translational Temperature
%                    - iBinsStart:       (Optional) Initial Level/Bin to read. If not given, it is assumed equal to be equal to 1
%                    - iBinsEnd:         (Optional) Final Level/Bin to read. If not given, it is assumed equal to be equal to NBins
%
%  Input Global Var: - ReadAllRatesFlg   Flag = -1/0/1/2; if == -1, rates are taken from the .mat file;
%                                                         if ==  0, each of the Levels/Bins has its own rates file for each of the temperatures; 
%                                                         if ==  1, each temperature has its own rates file containing all the Levels/Bins; 
%                                                         if ==  2, there is one rates file for all the Levels/Bins at all the temperatures.
%                    - SaveRatesFlg      Flag = 0/1;      if ==  1, after being uploaded, rates are saved in a .mat file.
%                    - T0_Vec:           Vector of Translational Temperatures (e.g.: [10000])
%                    - RatesPath:        The path to the output folder (e.g.: ../Test/N3/N2/Rates)
%                    - MoleculesName:    A vector of strings containing the name of the molecules present in the system
%                    - NBins             Nb of Levels/Bins
%

%function [RatesSigma, RatesMatrix, RatesSigmaMatrix, DissRates, DissRatesSigma, ProcessesRates, ProcessesRatesSigma] = ReadRates(iT, RatesSigma, RatesMatrix, RatesSigmaMatrix, DissRates, DissRatesSigma, ProcessesRates, ProcessesRatesSigma, iBinsStart, iBinsEnd)    
function [RatesMatrix, DissRates, ProcessesRates, KRec, KRecOverg] = ReadRates(iT, RatesMatrix, DissRates, ProcessesRates, iBinsStart, iBinsEnd, DissCorrectionFactor, NLevels, Levelg, LevelEeV)

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

  global NBins T0_Vec RatesPath ProcToLevIP ReadAllRatesFlg ProduceMatFlg ComponentMass DSWtoKg AvN UKb Ue Plnck RxLxIdx ComponentDeg NComp
  
  
  if (~exist('iBinsStart', 'var'))
    iBinsStart = ReadRatesSingleBin
  end
  if (~exist('iBinsStart', 'var'))
    iBinsEnd   = NBins(1)
  end
  
  
  if ReadAllRatesFlg == -1
    
    filenameRates = strcat(RatesPath,'/T_',num2str(T0_Vec(iT)),'_',num2str(T0_Vec(iT)),'/Rates.mat')
    %load(filenameRates, 'RatesSigma', 'RatesMatrix', 'RatesSigmaMatrix', 'DissRates', 'DissRatesSigma', 'ProcessesRates', 'ProcessesRatesSigma')
    load(filenameRates, 'RatesMatrix', 'DissRates', 'ProcessesRates','T0_VecTemp')
    %if T0_VecTemp ~= T0_Vec
    %  pause
    %end
    
  else

    if ReadAllRatesFlg == 0

      [RatesMatrix, DissRates, ProcessesRates] = ReadRatesSingleBin(iT, RatesMatrix, DissRates, ProcessesRates, iBinsStart, iBinsEnd);

    elseif ReadAllRatesFlg == 1

      [RatesMatrix, DissRates, ProcessesRates] = ReadRatesAll(iT, RatesMatrix, DissRates, ProcessesRates, iBinsStart, iBinsEnd);

    elseif ReadAllRatesFlg == 2

      [RatesMatrix, DissRates, ProcessesRates] = ReadRatesAllTint(RatesMatrix, DissRates, ProcessesRates, iBinsStart, iBinsEnd);

    end
    
    DissRates           = DissRates           .* DissCorrectionFactor;
    ProcessesRates(:,1) = ProcessesRates(:,1) .* DissCorrectionFactor;
    
    if ProduceMatFlg == 1
      
      T0_VecTemp = T0_Vec;
      
      filenameRates = strcat(RatesPath,'/T_',num2str(T0_Vec(iT)),'_',num2str(T0_Vec(iT)),'/Rates')
      %save(filenameRates,'RatesSigma', 'RatesMatrix', 'RatesSigmaMatrix', 'DissRates', 'DissRatesSigma', 'ProcessesRates', 'ProcessesRatesSigma','-v7.3');
      save(filenameRates,'RatesMatrix', 'DissRates', 'ProcessesRates','T0_VecTemp','-v7.3');
      
    end

  end
  
    iBinnedMol=1;
  
    ExpVec(1:NLevels(iBinnedMol),1)  = Levelg(1:NLevels(iBinnedMol),1) .* exp( - LevelEeV(1:NLevels(iBinnedMol),1) .* Ue ./ (T0_Vec(1) .* UKb) );
    %ExpVec                           = ExpVec ./ sum(ExpVec);

    iReac = 1;
    for iComp = 1:NComp
      QTran(iComp) = (2.d0 .* pi .* ComponentMass(iComp) .* DSWtoKg ./ AvN .* UKb .* T0_Vec(iT) ./ Plnck.^2).^(3/2 .* RxLxIdx(iComp));
    end
    IntDeg(:) = ComponentDeg(:).^(RxLxIdx(:));
    EqVec     = prod(QTran) .* prod(IntDeg) .* ExpVec(1:NLevels(iBinnedMol),1);


    KEq(1:NLevels(iBinnedMol),1)       = ( EqVec(1:NLevels(iBinnedMol)) ).^(-1);
    KRec(1:NLevels(iBinnedMol),1)      = 0.0;%DissRates(1:NLevels(iBinnedMol),1) ./ KEq(1:NLevels(iBinnedMol),1);
    KRecOverg(1:NLevels(iBinnedMol),1) = 0.0;%KRec(1:NLevels(iBinnedMol),1)      ./ Levelg(1:NLevels(iBinnedMol),1);
  
end