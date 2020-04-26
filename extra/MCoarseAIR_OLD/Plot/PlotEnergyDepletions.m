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

function [iFigure] = PlotEnergyDepletions(iFigure, t, CDInt, CDVib, CDRot, iQSS_Start, iQSS, iQSS_End)    

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
  
  global T0_Vec DissFlg InelFlg ExchFlg
  
  t_QSS     = t(iQSS);
  CDInt_QSS = CDInt(iQSS);
  CDVib_QSS = CDVib(iQSS);
  CDRot_QSS = CDRot(iQSS);

  
  figure(iFigure)
  semilogx(t, CDInt,'k')
  hold on
  semilogx(t_QSS, CDInt_QSS,'ko')
  semilogx(t, CDVib,'b')
  semilogx(t_QSS, CDVib_QSS,'bo')
  semilogx(t, CDRot,'r')
  semilogx(t_QSS, CDRot_QSS,'ro')
  iFigure = iFigure + 1;
  

  FileName = strcat('./CDEnergies_',num2str(DissFlg),'_',num2str(InelFlg),'_',num2str(ExchFlg),'_',num2str(ExchFlg),'.csv');
  if exist(FileName, 'file')
    fileID1  = fopen(FileName,'a');
  else
    fileID1  = fopen(FileName,'w');
    fprintf(fileID1,'# T [K], CDInt_Eq, CDVib_Eq, CDRot_Eq, CDInt_QSS, CDVib_QSS, CDRot_QSS\n');
  end
  fprintf(fileID1,'%e,%e,%e,%e,%e,%e,%e\n', T0_Vec(1), CDInt(end), CDVib(end), CDRot(end), CDInt_QSS, CDVib_QSS, CDRot_QSS )
  fclose(fileID1);
  
  
end