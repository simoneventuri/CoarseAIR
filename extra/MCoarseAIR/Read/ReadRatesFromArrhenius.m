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
function [RatesMatrix, DissRates, ProcessesRates] = ReadRatesFromArrhenius(iT, RatesMatrix, DissRates, ProcessesRates, iBinsStart, iBinsEnd)    

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

  global NBins T0_Vec RatesPath ProcToLevIP FldrStr
  

  opts = delimitedTextImportOptions("NumVariables", 7);
  opts.DataLines = [1, Inf];
  opts.Delimiter = ["(", ")", ",", ":"];
  opts.VariableNames = ["Var1", "VarName2", "Var3", "E12", "Var5", "Var6", "Var7"];
  opts.SelectedVariableNames = ["VarName2", "E12"];
  opts.VariableTypes = ["string", "double", "string", "double", "string", "string", "string"];
  opts = setvaropts(opts, [1, 3, 5, 6, 7], "WhitespaceRule", "preserve");
  opts = setvaropts(opts, [1, 3, 5, 6, 7], "EmptyFieldRule", "auto");
  opts.ExtraColumnsRule = "ignore";
  opts.EmptyLineRule = "read";
  tbl = readtable(strcat("/home/venturi/WORKSPACE/Mars_Database/Run_0D/database/kinetics/O3_UMN",FldrStr,"/T",num2str((T0_Vec(iT))),"K/Diss_Corrected.dat"), opts)
  Idx      = tbl.VarName2;
  DissTemp = tbl.E12;
  clear opts tbl
  ProcessesRates(:,1,iT) = DissTemp(Idx(:));
  DissRates              = ProcessesRates(:,1,iT);
  
  
  opts = delimitedTextImportOptions("NumVariables", 9);
  opts.DataLines = [1, Inf];
  opts.Delimiter = ["(", ")", ",", ":"];
  opts.VariableNames = ["Var1", "VarName2", "Var3", "VarName4", "Var5", "e11", "Var7", "Var8", "Var9"];
  opts.SelectedVariableNames = ["VarName2", "VarName4", "e11"];
  opts.VariableTypes = ["string", "double", "string", "double", "string", "double", "string", "string", "string"];
  opts = setvaropts(opts, [1, 3, 5, 7, 8, 9], "WhitespaceRule", "preserve");
  opts = setvaropts(opts, [1, 3, 5, 7, 8, 9], "EmptyFieldRule", "auto");
  opts.ExtraColumnsRule = "ignore";
  opts.EmptyLineRule = "read";
  tbl = readtable(strcat("/home/venturi/WORKSPACE/Mars_Database/Run_0D/database/kinetics/O3_UMN",FldrStr,"/T",num2str((T0_Vec(iT))),"K/Inel.dat"), opts);
  iIdx = tbl.VarName2;
  jIdx = tbl.VarName4;
  Rate = tbl.e11;
  for iLine = 1:length(Rate)
    RatesMatrix(iIdx(iLine),jIdx(iLine),1,iT) = Rate(iLine);
  end
  clear opts tbl iIdx jIdx Rate
  
  
  opts = delimitedTextImportOptions("NumVariables", 9);
  opts.DataLines = [1, Inf];
  opts.Delimiter = ["(", ")", ",", ":"];
  opts.VariableNames = ["Var1", "VarName2", "Var3", "VarName4", "Var5", "e11", "Var7", "Var8", "Var9"];
  opts.SelectedVariableNames = ["VarName2", "VarName4", "e11"];
  opts.VariableTypes = ["string", "double", "string", "double", "string", "double", "string", "string", "string"];
  opts = setvaropts(opts, [1, 3, 5, 7, 8, 9], "WhitespaceRule", "preserve");
  opts = setvaropts(opts, [1, 3, 5, 7, 8, 9], "EmptyFieldRule", "auto");
  opts.ExtraColumnsRule = "ignore";
  opts.EmptyLineRule = "read";
  tbl = readtable(strcat("/home/venturi/WORKSPACE/Mars_Database/Run_0D/database/kinetics/O3_UMN",FldrStr,"/T",num2str((T0_Vec(iT))),"K/Exch_Type1.dat"), opts);
  iIdx = tbl.VarName2;
  jIdx = tbl.VarName4;
  Rate = tbl.e11;
  for iLine = 1:length(Rate)
    RatesMatrix(iIdx(iLine),jIdx(iLine),2,iT) = Rate(iLine)./2.0;
    RatesMatrix(iIdx(iLine),jIdx(iLine),3,iT) = Rate(iLine)./2.0;
  end
  clear opts tbl iIdx jIdx Rate
  ProcessesRates(:,2,iT) = sum(RatesMatrix(:,:,2,iT),2);
  ProcessesRates(:,3,iT) = sum(RatesMatrix(:,:,3,iT),2);
  
  
end