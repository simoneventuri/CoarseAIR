%% The Function Loads the Variables of the Chemical System
%                       
function Initialize_ChemicalSyst()

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
  
    global Syst Input

    
    fprintf('= Initialize_ChemicalSyst ==========================\n')
    fprintf('====================================================\n')
    fprintf(['Loading Variables for Chemical System:' ' ' Input.SystNameLong '\n'])
    
    
    if strcmp(Input.SystNameLong, 'CO2_NASA')
        
        %%% System
        Syst.Name              = 'CO2';
        Syst.NameLong_Opposite = 'O2C_NASA';

        Syst.NProc = 4; %(Diss+Inel+Exch+Exch)
        
        
        %%% Atoms
        Syst.NAtoms        = 3;
        
        Syst.Atom(1).Name  = 'C';
        Syst.Atom(2).Name  = 'O';
        Syst.Atom(3).Name  = 'O';

        Syst.Atom(1).Color = [0, 0, 0];
        Syst.Atom(2).Color = [0, 0, 1];
        Syst.Atom(3).Color = [0, 0, 1];

        Syst.Atom(1).Size  = 100;
        Syst.Atom(2).Size  = 200;
        Syst.Atom(3).Size  = 200;

        Syst.Atom(1).Mass  = 21868.661757;
        Syst.Atom(2).Mass  = 29148.94559;
        Syst.Atom(3).Mass  = 29148.94559;
        
        
        %%% Molecules
        Syst.NMolecules                   = 2;
        
        Syst.Molecule(1).Name             = 'CO';
        Syst.Molecule(2).Name             = 'O2';

        Syst.Molecule(1).DissEn           = 0.0;
        Syst.Molecule(2).DissEn           = 0.0;

        Syst.Molecule(1).DegeneracyFactor = [  1,   1];
        Syst.Molecule(2).DegeneracyFactor = [1/2, 1/2];

        Syst.Molecule(1).Mu               = 28.0104e-3;
        Syst.Molecule(2).Mu               = 31.9988e-3;

        Syst.Molecule(1).NLevelsOrig      = 13521;
        Syst.Molecule(2).NLevelsOrig      = 6078;
    
        Syst.Molecule(1).ToAtoms          = [1,2];
        Syst.Molecule(2).ToAtoms          = [2,3];
        
        Syst.Molecule(1).DiatPot          = 'CO_NASA';
        Syst.Molecule(2).DiatPot          = 'O2_NASA';
                
        %%% Pairs
        Syst.Pair(1).Name  = 'CO';
        Syst.Pair(2).Name  = 'CO';
        Syst.Pair(3).Name  = 'O2';

        Syst.Pair(1).ToMol = 1;
        Syst.Pair(2).ToMol = 1;
        Syst.Pair(3).ToMol = 2;

        Syst.Pair(1).Color = [17, 17, 17] ./ 256;
        Syst.Pair(2).Color = [17, 17, 17] ./ 256;
        Syst.Pair(3).Color = [0, 0, 256]  ./ 256;

        
        %% CFD Components (For PLATO and KONIG)
        Syst.NComp             =  4;
        
        Syst.CFDComp(1).Name   = 'C';
        Syst.CFDComp(2).Name   = 'O';
        Syst.CFDComp(3).Name   = 'CO';
        Syst.CFDComp(4).Name   = 'O2';

        Syst.CFDComp(1).ToMol   = 0;
        Syst.CFDComp(2).ToMol   = 0;
        Syst.CFDComp(3).ToMol   = 1;
        Syst.CFDComp(4).ToMol   = 2;

        Syst.CFDComp(1).Mass    = Syst.Atom(1).Mass;
        Syst.CFDComp(2).Mass    = Syst.Atom(2).Mass;
        Syst.CFDComp(3).Mass    = Syst.Atom(1).Mass + Syst.Atom(2).Mass;
        Syst.CFDComp(4).Mass    = 2.0 * Syst.Atom(3).Mass;

        Syst.CFDComp(1).Qe      = 1;
        Syst.CFDComp(2).Qe      = 9;
        Syst.CFDComp(3).Qe      = 1;
        Syst.CFDComp(4).Qe      = 3;

        Syst.CFDComp(1).Color   = [ 102, 102, 102] ./ 256;
        Syst.CFDComp(2).Color   = [   0,   0,   0] ./ 256;
        Syst.CFDComp(3).Color   = [ 204,   0,   0] ./ 256;
        Syst.CFDComp(4).Color   = [   0,   0, 234] ./ 256;

        Syst.CFDComp(1).LineStyle = ':';
        Syst.CFDComp(2).LineStyle = '-.';
        Syst.CFDComp(3).LineStyle = '-';
        Syst.CFDComp(4).LineStyle = '--';
        
        Syst.RxLxIdx = [-1, -1, 1,  0;   % Diss
                         0,  0, 0,  0;   % Inel 
                         0,  0, 0,  0;   % Exch1
                        -1,  1, 1, -1];  % Exch2
        
        Syst.MolToCFDComp       = [3, 4];

        
        %% Exchange Properties
        Syst.ExchToMol          = [1; 2];
        Syst.ExchToAtom         = [3; 1];
        
        Syst.PairToExch         = [1; 2];
        
        Syst.ToOtherExch        = [0; 1];
        
        Syst.ColPartToComp      = 2; 
        
        
    elseif strcmp(Input.SystNameLong, 'O2C_NASA')
        
        %%% System
        Syst.Name              = 'O2C';
        Syst.NameLong_Opposite = 'CO2_NASA';

        Syst.NProc = 3; %(Diss+Inel+Exch)
        
        
        %%% Atoms
        Syst.NAtoms        = 3;
        
        Syst.Atom(1).Name  = 'O';
        Syst.Atom(2).Name  = 'O';
        Syst.Atom(3).Name  = 'C';

        Syst.Atom(1).Color = [0, 0, 1];
        Syst.Atom(2).Color = [0, 0, 1];
        Syst.Atom(3).Color = [0, 0, 0];

        Syst.Atom(1).Size  = 200;
        Syst.Atom(2).Size  = 200;
        Syst.Atom(3).Size  = 100;

        Syst.Atom(1).Mass  = 29148.94559;
        Syst.Atom(2).Mass  = 29148.94559;
        Syst.Atom(3).Mass  = 21868.661757;
        
        
        %%% Molecules
        Syst.NMolecules                   = 2;
        
        Syst.Molecule(1).Name             = 'O2';
        Syst.Molecule(2).Name             = 'CO';

        Syst.Molecule(1).DissEn           = 0.0;
        Syst.Molecule(2).DissEn           = 0.0;

        Syst.Molecule(1).DegeneracyFactor = [1/2, 1/2];
        Syst.Molecule(2).DegeneracyFactor = [  1,   1];

        Syst.Molecule(1).Mu               = 31.9988e-3;
        Syst.Molecule(2).Mu               = 28.0104e-3;

        Syst.Molecule(1).NLevelsOrig      = 6078;
        Syst.Molecule(2).NLevelsOrig      = 13521;
    
        Syst.Molecule(1).ToAtoms          = [1,2];
        Syst.Molecule(2).ToAtoms          = [2,3];
        
        Syst.Molecule(1).DiatPot          = 'O2_NASA';
        Syst.Molecule(2).DiatPot          = 'CO_NASA';
        
        
        %%% Pairs
        Syst.Pair(1).Name  = 'O2';
        Syst.Pair(2).Name  = 'CO';
        Syst.Pair(3).Name  = 'CO';

        Syst.Pair(1).ToMol = 1;
        Syst.Pair(2).ToMol = 2;
        Syst.Pair(3).ToMol = 2;

        Syst.Pair(1).Color = [0, 0, 256]  ./ 256;
        Syst.Pair(2).Color = [17, 17, 17] ./ 256;
        Syst.Pair(3).Color = [17, 17, 17] ./ 256;

        
        %% CFD Components (For PLATO and KONIG)
        Syst.NComp             =  4;
        
        Syst.CFDComp(1).Name   = 'C';
        Syst.CFDComp(2).Name   = 'O';
        Syst.CFDComp(3).Name   = 'CO';
        Syst.CFDComp(4).Name   = 'O2';

        Syst.CFDComp(1).ToMol   = 0;
        Syst.CFDComp(2).ToMol   = 0;
        Syst.CFDComp(3).ToMol   = 2;
        Syst.CFDComp(4).ToMol   = 1;

        Syst.CFDComp(1).Mass    = Syst.Atom(3).Mass;
        Syst.CFDComp(2).Mass    = Syst.Atom(1).Mass;
        Syst.CFDComp(3).Mass    = Syst.Atom(1).Mass + Syst.Atom(3).Mass;
        Syst.CFDComp(4).Mass    = 2.0 * Syst.Atom(1).Mass;

        Syst.CFDComp(1).Qe      = 1;
        Syst.CFDComp(2).Qe      = 9;
        Syst.CFDComp(3).Qe      = 1;
        Syst.CFDComp(4).Qe      = 3;

        Syst.CFDComp(1).Color   = [ 102, 102, 102] ./ 256;
        Syst.CFDComp(2).Color   = [   0,   0,   0] ./ 256;
        Syst.CFDComp(3).Color   = [ 204,   0,   0] ./ 256;
        Syst.CFDComp(4).Color   = [   0,   0, 234] ./ 256;

        Syst.CFDComp(1).LineStyle = ':';
        Syst.CFDComp(2).LineStyle = '-.';
        Syst.CFDComp(3).LineStyle = '-';
        Syst.CFDComp(4).LineStyle = '--';

        Syst.RxLxIdx = [ 0, -2,  0,  1;   % Diss
                         0,  0,  0,  0;   % Inel 
                         1, -1, -1,  1];  % Exch1
        
        Syst.MolToCFDComp       = [4, 3];

        
        %% Exchange Properties
        Syst.ExchToMol          = [2];
        Syst.ExchToAtom         = [3];

        Syst.PairToExch         = [2];
        
        Syst.ToOtherExch        = [2];
        
        Syst.ColPartToComp      = 1; 
        
        
    elseif strcmp(Input.SystNameLong, 'N2O_UMN')
        
        %%% System
        Syst.Name              = 'N2O';
        Syst.NameLong_Opposite = 'NON_UMN';

        Syst.NProc = 3; %(Diss+Inel+Exch)
        
        
        %%% Atoms
        Syst.NAtoms        = 3;
        
        Syst.Atom(1).Name  = 'N';
        Syst.Atom(2).Name  = 'N';
        Syst.Atom(3).Name  = 'O';

        Syst.Atom(1).Color = [0, 0, 1];
        Syst.Atom(2).Color = [0, 0, 1];
        Syst.Atom(3).Color = [0, 0, 0];

        Syst.Atom(1).Size  = 200;
        Syst.Atom(2).Size  = 200;
        Syst.Atom(3).Size  = 100;

        Syst.Atom(1).Mass  = 25526.04298;
        Syst.Atom(2).Mass  = 25526.04298;
        Syst.Atom(3).Mass  = 29148.94559;
        
        
        %%% Molecules
        Syst.NMolecules                   = 2;
        
        Syst.Molecule(1).Name             = 'N2';
        Syst.Molecule(2).Name             = 'NO';

        Syst.Molecule(1).DissEn           = 0.0;
        Syst.Molecule(2).DissEn           = 0.0;

        Syst.Molecule(1).DegeneracyFactor = [ 3, 6]; %Odd, Even
        Syst.Molecule(2).DegeneracyFactor = [ 2, 2];

        Syst.Molecule(1).Mu               = 31.9988e-3;
        Syst.Molecule(2).Mu               = 30.0061e-3;

        Syst.Molecule(1).NLevelsOrig      = 9093;
        Syst.Molecule(2).NLevelsOrig      = 6793;
    
        Syst.Molecule(1).ToAtoms          = [1,2];
        Syst.Molecule(2).ToAtoms          = [2,3];
        
        Syst.Molecule(1).DiatPot          = 'N2_UMN';
        Syst.Molecule(2).DiatPot          = 'NO_UMN';
        
        
        %%% Pairs
        Syst.Pair(1).Name  = 'N2';
        Syst.Pair(2).Name  = 'NO';
        Syst.Pair(3).Name  = 'NO';

        Syst.Pair(1).ToMol = 1;
        Syst.Pair(2).ToMol = 2;
        Syst.Pair(3).ToMol = 2;

        Syst.Pair(1).Color = [0, 0, 256]  ./ 256;
        Syst.Pair(2).Color = [17, 17, 17] ./ 256;
        Syst.Pair(3).Color = [17, 17, 17] ./ 256;

        
        %% CFD Components (For PLATO and KONIG)
        Syst.NComp             =  4;
        
        Syst.CFDComp(1).Name   = 'N';
        Syst.CFDComp(2).Name   = 'O';
        Syst.CFDComp(3).Name   = 'N2';
        Syst.CFDComp(4).Name   = 'NO';

        Syst.CFDComp(1).ToMol   = 0;
        Syst.CFDComp(2).ToMol   = 0;
        Syst.CFDComp(3).ToMol   = 1;
        Syst.CFDComp(4).ToMol   = 2;

        Syst.CFDComp(1).Mass    = Syst.Atom(1).Mass;
        Syst.CFDComp(2).Mass    = Syst.Atom(3).Mass;
        Syst.CFDComp(3).Mass    = 2.0 * Syst.Atom(1).Mass;
        Syst.CFDComp(4).Mass    = Syst.Atom(1).Mass + Syst.Atom(3).Mass;

        Syst.CFDComp(1).Qe      = 12;
        Syst.CFDComp(2).Qe      = 9;
        Syst.CFDComp(3).Qe      = 1;
        Syst.CFDComp(4).Qe      = 4;

        Syst.CFDComp(1).Color   = [ 102, 102, 102] ./ 256;
        Syst.CFDComp(2).Color   = [   0,   0,   0] ./ 256;
        Syst.CFDComp(3).Color   = [ 204,   0,   0] ./ 256;
        Syst.CFDComp(4).Color   = [   0,   0, 234] ./ 256;

        Syst.CFDComp(1).LineStyle = ':';
        Syst.CFDComp(2).LineStyle = '-.';
        Syst.CFDComp(3).LineStyle = '-';
        Syst.CFDComp(4).LineStyle = '--';

        Syst.RxLxIdx = [-2,  0,  1,  0;   % Diss
                         0,  0,  0,  0;   % Inel 
                        -1,  1,  1, -1];  % Exch1
        
        Syst.MolToCFDComp       = [3, 4];

        
        %% Exchange Properties
        Syst.ExchToMol          = [2];
        Syst.ExchToAtom         = [2];

        Syst.PairToExch         = [2];
        
        Syst.ToOtherExch        = [1];
        
        Syst.ColPartToComp      = 2; 
        
        
        
    elseif strcmp(Input.SystNameLong, 'NON_UMN')
        
        %%% System
        Syst.Name              = 'NON';
        Syst.NameLong_Opposite = 'N2O_UMN';
        
        Syst.NProc = 4; %(Diss+Inel+Exch+Exch)
        
        
        %%% Atoms
        Syst.NAtoms        = 3;
        
        Syst.Atom(1).Name  = 'N';
        Syst.Atom(2).Name  = 'O';
        Syst.Atom(3).Name  = 'N';

        Syst.Atom(1).Color = [0, 0, 1];
        Syst.Atom(2).Color = [0, 0, 0];
        Syst.Atom(3).Color = [0, 0, 1];

        Syst.Atom(1).Size  = 200;
        Syst.Atom(2).Size  = 100;
        Syst.Atom(3).Size  = 200;

        Syst.Atom(1).Mass  = 25526.04298;
        Syst.Atom(2).Mass  = 29148.94559;
        Syst.Atom(3).Mass  = 25526.04298;
        
        
        %%% Molecules
        Syst.NMolecules                   = 2;
        
        Syst.Molecule(1).Name             = 'NO';
        Syst.Molecule(2).Name             = 'N2';

        Syst.Molecule(1).DissEn           = 0.0;
        Syst.Molecule(2).DissEn           = 0.0;

        Syst.Molecule(1).DegeneracyFactor = [2, 2];
        Syst.Molecule(2).DegeneracyFactor = [3, 6];

        Syst.Molecule(1).Mu               = 30.0061e-3;
        Syst.Molecule(2).Mu               = 31.9988e-3;

        Syst.Molecule(1).NLevelsOrig      = 6793;
        Syst.Molecule(2).NLevelsOrig      = 9093;
    
        Syst.Molecule(1).ToAtoms          = [1,2];
        Syst.Molecule(2).ToAtoms          = [1,3];
        
        Syst.Molecule(1).DiatPot          = 'NO_UMN';
        Syst.Molecule(2).DiatPot          = 'N2_UMN';
        
        
        %%% Pairs
        Syst.Pair(1).Name  = 'NO';
        Syst.Pair(2).Name  = 'N2';
        Syst.Pair(3).Name  = 'NO';

        Syst.Pair(1).ToMol = 1;
        Syst.Pair(2).ToMol = 2;
        Syst.Pair(3).ToMol = 1;

        Syst.Pair(1).Color = [17, 17, 17] ./ 256;
        Syst.Pair(2).Color = [0, 0, 256]  ./ 256;
        Syst.Pair(3).Color = [17, 17, 17] ./ 256;

        
        %% CFD Components (For PLATO and KONIG)
        Syst.NComp             =  4;
        
        Syst.CFDComp(1).Name   = 'N';
        Syst.CFDComp(2).Name   = 'O';
        Syst.CFDComp(3).Name   = 'N2';
        Syst.CFDComp(4).Name   = 'NO';

        Syst.CFDComp(1).ToMol   = 0;
        Syst.CFDComp(2).ToMol   = 0;
        Syst.CFDComp(3).ToMol   = 2;
        Syst.CFDComp(4).ToMol   = 1;

        Syst.CFDComp(1).Mass    = Syst.Atom(1).Mass;
        Syst.CFDComp(2).Mass    = Syst.Atom(2).Mass;
        Syst.CFDComp(3).Mass    = 2.0 * Syst.Atom(3).Mass;
        Syst.CFDComp(4).Mass    = Syst.Atom(1).Mass + Syst.Atom(2).Mass;

        Syst.CFDComp(1).Qe      = 12;
        Syst.CFDComp(2).Qe      = 9;
        Syst.CFDComp(3).Qe      = 1;
        Syst.CFDComp(4).Qe      = 4;

        Syst.CFDComp(1).Color   = [ 102, 102, 102] ./ 256;
        Syst.CFDComp(2).Color   = [   0,   0,   0] ./ 256;
        Syst.CFDComp(3).Color   = [ 204,   0,   0] ./ 256;
        Syst.CFDComp(4).Color   = [   0,   0, 234] ./ 256;

        Syst.CFDComp(1).LineStyle = ':';
        Syst.CFDComp(2).LineStyle = '-.';
        Syst.CFDComp(3).LineStyle = '-';
        Syst.CFDComp(4).LineStyle = '--';

        Syst.RxLxIdx = [-1, -1,  0,  1;   % Diss
                         0,  0,  0,  0;   % Inel 
                         1, -1, -1,  1;   % Exch1
                         0,  0,  0,  0];  % Exch2 

        Syst.MolToCFDComp       = [4, 3];

        
        %% Exchange Properties
        Syst.ExchToMol          = [2, 1];
        Syst.ExchToAtom         = [2, 1];

        Syst.PairToExch         = [2, 3];

        Syst.ToOtherExch        = [0, 1];
        
        Syst.ColPartToComp      = 1; 
        

    elseif strcmp(Input.SystNameLong, 'N4_NASA')

        %%% System
        Syst.Name = 'N4';
        
        Syst.NProc = 3; %(Diss+Inel+Exch)
        
        
        %%% Atoms
        Syst.NAtoms        = 4;
        
        Syst.Atom(1).Name  = 'N';
        Syst.Atom(2).Name  = 'N';
        Syst.Atom(3).Name  = 'N';
        Syst.Atom(4).Name  = 'N';

        Syst.Atom(1).Color = [0, 0, 0];
        Syst.Atom(2).Color = [0, 0, 1];
        Syst.Atom(3).Color = [0, 0, 1];
        Syst.Atom(4).Color = [0, 0, 1];

        Syst.Atom(1).Size  = 150;
        Syst.Atom(2).Size  = 150;
        Syst.Atom(3).Size  = 150;
        Syst.Atom(4).Size  = 150;

        Syst.Atom(1).Mass  = 25526.04298d0;
        Syst.Atom(2).Mass  = 25526.04298d0;
        Syst.Atom(3).Mass  = 25526.04298d0;
        Syst.Atom(4).Mass  = 25526.04298d0;


        %%% Molecules
        Syst.NMolecules                   = 1;
        
        Syst.Molecule(1).Name             = 'N2';
        
        Syst.Molecule(1).DissEn           = 0.0;
        
        Syst.Molecule(1).DegeneracyFactor = [3, 6];
        
        Syst.Molecule(1).Mu               = 28.0134d-3;
        
        Syst.Molecule(1).NLevelsOrig      = 9390;
        
        Syst.Molecule(1).ToAtoms          = [1,2];
        
        Syst.Molecule(1).DiatPot          = 'N2_NASA';
        
        
        %%% Pairs
        Syst.Pair(1).Name  = 'N2';
        Syst.Pair(2).Name  = 'N2';
        Syst.Pair(3).Name  = 'N2';
        Syst.Pair(4).Name  = 'N2';
        Syst.Pair(5).Name  = 'N2';
        Syst.Pair(6).Name  = 'N2';
        
        Syst.Pair(1).ToMol = 1;
        Syst.Pair(2).ToMol = 1;
        Syst.Pair(3).ToMol = 1;
        Syst.Pair(4).ToMol = 1;
        Syst.Pair(5).ToMol = 1;
        Syst.Pair(6).ToMol = 1;
        
        Syst.Pair(1).Color = [17, 17, 17] ./ 256;
        Syst.Pair(2).Color = [17, 17, 17] ./ 256;
        Syst.Pair(3).Color = [17, 17, 17] ./ 256;
        Syst.Pair(4).Color = [17, 17, 17] ./ 256;
        Syst.Pair(5).Color = [17, 17, 17] ./ 256;
        Syst.Pair(6).Color = [17, 17, 17] ./ 256;
        
        
        %% CFD Components (For PLATO and KONIG)
        Syst.NComp             =  2;
        
        Syst.CFDComp(1).Name   = 'N';
        Syst.CFDComp(2).Name   = 'N2';

        Syst.CFDComp(1).ToMol   = 0;
        Syst.CFDComp(2).ToMol   = 1;

        Syst.CFDComp(1).Mass    = Syst.Atom(1).Mass;
        Syst.CFDComp(2).Mass    = 2.0 * Syst.Atom(1).Mass;

        Syst.CFDComp(1).Qe      = 1;
        Syst.CFDComp(2).Qe      = 1;

        Syst.CFDComp(1).Color   = [ 102, 102, 102] ./ 256;
        Syst.CFDComp(2).Color   = [   0,   0,   0] ./ 256;

        Syst.CFDComp(1).LineStyle = ':';
        Syst.CFDComp(2).LineStyle = '-';
        
        Syst.RxLxIdx = [-4, 2;   % Double Diss
                         2, 1;   % Single Diss
                         0, 0;   % Inel
                         0, 0];  % Exch1
        
        
        %% Exchange Properties
        Syst.MolToCFDComp       = [2];

        Syst.ExchToMol          = [1,1];
        Syst.ExchToAtom         = [0,0];

        Syst.ColPartToComp      = 1; 
        
  
    end
    
    Syst
    fprintf('====================================================\n\n')  

end