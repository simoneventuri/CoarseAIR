%% The Function loads the physical constants and the parameters for plotting
%
%  Input Global Var: - System:        The chemical system of interest (e.g.: N3 / CO2 / O3 / etc.)
%                    - PathToOutput:  The path to the output folder (e.g.: ../Test/ )
%                    - MoleculesName: A vector of strings containing the name of the molecules present in the system
%                    - RunKONIGDir:   (Optional for PostCQCT.m) Path to KONIG's run directory
%

function [iFigure] = UpdateParameters()

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
  
  global System PathToOutput MoleculesName RunKONIGDir Pair_to_Atoms
  
  global Plnck UKb Ue KeV AvN AMUToKg EhToeV DSWtoKg ATMToPa
  
  global SystemPath RatesPath DatabasePath OutputPath SaveFigs linS linST FigDirPath AxisFontSz AxisFontNm AxisLabelSz AxisLabelNm LegendFontSz LegendFontNm  
  
  global RCVec BCVec GCVec KCVec OCVec PCVec WCVec JCVec YCVec CCVec MCVec
  
  Plnck    = 6.62607004d-34;
  UKb      = 1.380658e-23;
  Ue       = 1.602191e-19;
  KeV      = 8.617330e-05;
  AvN      = 6.0221409e+23;
  AMUToKg  = 1.d0/AvN*1.d-3;
  EhToeV   = 27.2114;
  DSWtoKg  = 1.d-3/1.8208e+03;
  ATMToPa  = 101325.d0;

  Pair_to_Atoms = [1,2;1,3;2,3];
  
  System(1,:)
  SystemPath   = strcat(PathToOutput,'/',System(1,:));                                               % For qnsEnBin.dat (Levels Info)
  RatesPath    = strcat(PathToOutput,'/',System(1,:),'/',MoleculesName(1,:),'/Rates/');              % For Levels / Bins Rates
  DatabasePath = strcat(PathToOutput,'/',RunKONIGDir,'/database/');                                  % For thermo (Bins Info)
  OutputPath   = strcat(PathToOutput,'/',RunKONIGDir,'/output/');                                    % For box.dat, Tint.dat, 
  
  iFigure  = 1;
  SaveFigs = 0;
  
  linS     = {'--','-.',':','-'};
  linST    = {'-','-.',':','--'};
  
%   AxisFontSz = 36;
%   AxisFontNm = 'Arial';
%   
%   AxisLabelSz = 44;
%   AxisLabelNm = 'Arial';
%   
%   LegendFontSz = 40;
%   LegendFontNm = 'Arial';
  
  AxisFontSz = 26;
  AxisFontNm = 'Times';
  
  AxisLabelSz = 30;
  AxisLabelNm = 'Times';
  
  LegendFontSz = 25;
  LegendFontNm = 'Times';
  
  RCVec = [255  50  20] ./ 255;
  BCVec = [  0  70 200] ./ 255;
  GCVec = [  0 140  50] ./ 255;
  KCVec = [  0   0   0] ./ 255;
  OCVec = [255 105  45] ./ 255;
  PCVec = [155  45 175] ./ 255;
  WCVec = [  1   1   1] ./ 255;
  JCVec = [100 100 100] ./ 255;
  YCVec = [255 255   0] ./ 255;
  CCVec = [205 205 205] ./ 255;
  MCVec = [100  25  15] ./ 255;
  
end