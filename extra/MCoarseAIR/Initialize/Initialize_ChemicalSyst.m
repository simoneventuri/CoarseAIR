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
        Syst.Name  = 'CO2';
        
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

        Syst.Molecule(1).DegeneracyFactor = 1;
        Syst.Molecule(2).DegeneracyFactor = 6;

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

        Syst.CFDComp(1).Deg     = 1;
        Syst.CFDComp(2).Deg     = 9;
        Syst.CFDComp(3).Deg     = 1;
        Syst.CFDComp(4).Deg     = 3;

        Syst.CFDComp(1).Color   = [ 102, 102, 102] ./ 256;
        Syst.CFDComp(2).Color   = [   0,   0,   0] ./ 256;
        Syst.CFDComp(3).Color   = [ 204,   0,   0] ./ 256;
        Syst.CFDComp(4).Color   = [   0,   0, 234] ./ 256;

        Syst.CFDComp(1).LineStyle = ':';
        Syst.CFDComp(2).LineStyle = '-.';
        Syst.CFDComp(3).LineStyle = '-';
        Syst.CFDComp(4).LineStyle = '--';

        Syst.CFDComp(1).RxLxIdx = -1;
        Syst.CFDComp(2).RxLxIdx = -1;
        Syst.CFDComp(3).RxLxIdx =  1;
        Syst.CFDComp(4).RxLxIdx =  0;
        
        Syst.MolToCFDComp       = [3, 4];

        
        %% Exchange Properties
        Syst.ExchToMol          = [1; 2];
        Syst.ExchToAtom         = [3; 1];

        Syst.PairToExch         = [1; 2];
        
        Syst.ColPartToComp      = 2; 
        

%     elseif strcmp(Input.SystNameLong, 'O2C')
%       Input.SystNameLong = 'CO2';
% 
%       NTint  = length(T0_Vec);
% 
%       Syst.NAtoms    = 3;
%       AtomsName = [  'O',   'O',   'C'];
%       AtomColor = [0 0 1; 0 0 1; 0 0 0];
%       AtomSize  = [  200,   200,   100];
%       AtomMass  = [  29148.94559d0, 29148.94559d0, 21868.661757d0]
% 
% 
%       MoleculesName    = ['O2'; 'CO'];
%       MoleculedDissEn  = [ 0.0;  0.0];
%       NMolecules       = length(MoleculesName);
%       DegeneracyFactor = [6; 1];
%       MoleculeMu       = [31.9988d-3, 28.0104e-3];
% 
%       AllMoleculesName = ['O2'; 'CO'; 'CO'];
% 
%       PairColor = [0 0 256; 17 17 17; 17 17 17]./256;
% 
%       BinnedMolName   = ['O2'; 'CO'];
%       NBinnedMol      = length(BinnedMolName);
%       BinnedMolToComp = [   4;    3];
% 
%       NComp = 4;
%         CompNames=['  ';'  ';'  ';'  '];
%         CompNames(1,1:2) = ' C';
%         CompNames(2,1:2) = ' O';
%         CompNames(3,1:2) = 'CO';
%         CompNames(4,1:2) = 'O2';
% 
%         ComponentMass = [ AtomMass(3), AtomMass(1), AtomMass(1)+AtomMass(3), 2.d0.*AtomMass(1)];
%         ComponentDeg  = [           1,           9,                       1,                 1];
% 
%       ColPartToComp = 2;
% 
%         CompColor(1,:) = [ 102, 102, 102];
%         CompColor(2,:) = [   0, 153, 102];
%         CompColor(3,:) = [ 204,   0,   0];
%         CompColor(4,:) = [   0,   0, 234];
%         CompColor      = CompColor ./ 256;
% 
%         ColorVec = CompColor
% 
%       RxLxIdx          = [-1, -1, 0, 1];
%       
%       
%     elseif strcmp(Input.SystNameLong, 'N2O')
% 
%       NTint  = length(T0_Vec);
% 
%       Syst.NAtoms    = 3;
%       AtomsName = [  'N',   'N',   'O'];
%       AtomColor = [0 0 0; 0 0 1; 0 0 1];
%       AtomSize  = [  200,   200,   200];
%       AtomMass  = [  25526.04298d0, 25526.04298d0, 229148.94559d0]
% 
%       MoleculesName    = ['N2_UMN'; 'NO_UMN'];
%       MoleculedDissEn  = [ 0.0;  0.0];
%       NMolecules       = length(MoleculesName);
%       DegeneracyFactor = [1; 6];
%       MoleculeMu       = [28.0104e-3, 31.9988d-3];
% 
%       AllMoleculesName = ['N2_UMN'; 'NO_UMN'; 'NO_UMN'];
% 
%       PairColor = [17 17 17;0 0 256; 0 0 256]./256;
% 
%       BinnedMolName   = ['N2_UMN'; 'NO_UMN'];
%       NBinnedMol      = length(BinnedMolName);
%       BinnedMolToComp = [   3;    4];
% 
%       NComp = 4;
%         CompNames=['  ';'  ';'  ';'  '];
%         CompNames(1,1:2) = ' N';
%         CompNames(2,1:2) = ' O';
%         CompNames(3,1:2) = 'N2';
%         CompNames(4,1:2) = 'NO';
% 
%         ComponentMass = [ AtomMass(1), AtomMass(2), AtomMass(1)+AtomMass(2), AtomMass(1)+AtomMass(3)];
%         ComponentDeg  = [           1,           9,                       1,                 1]
% 
%       ColPartToComp = 2;
% 
%         CompColor(1,:) = [ 102, 102, 102];
%         CompColor(2,:) = [   0, 153, 102];
%         CompColor(3,:) = [ 204,   0,   0];
%         CompColor(4,:) = [   0,   0, 234];
%         CompColor      = CompColor ./ 256;
% 
%       ColorVec = CompColor
% 
%       RxLxIdx          = [-1, -1, 1, 0];
% 
% 
%     elseif strcmp(Input.SystNameLong, 'N3_UMN')
% 
%       NTint  = length(T0_Vec);
% 
%       Syst.NAtoms    = 3;
%       AtomsName = [    'N',     'N',    'N'];
%       AtomColor = [0 0.3 0; 0 0.3 0; 0 0.3 0];
%       AtomSize  = [    150,     150,     150];
%       AtomMass  = [  25526.04298d0, 25526.04298d0, 25526.04298d0]
% 
% 
%       MoleculesName    = ['N2_UMN'];
%       MoleculedDissEn  = [ 0.0];
%       NMolecules       = size(MoleculesName,1);
%       DegeneracyFactor = (2);
%       MoleculeMu       = [28.0134d-3];
% 
% 
%       AllMoleculesName = ['N2_UMN'; 'N2_UMN'; 'N2_UMN'];
% 
%       PairColor = [0 0 256; 0 256 0; 256 0 0]./256;
% 
%       BinnedMolName   = ['N2_UMN'];
%       NBinnedMol      = size(BinnedMolName,1);
%       BinnedMolToComp = [   2];
% 
%       NComp = 2;
%         CompNames=['  ';'  '];
%         CompNames(1,:) = ' N';
%         CompNames(2,:) = 'N2';
% 
%         ComponentMass = [AtomMass(1), 2.d0.*AtomMass(1)];
%         ComponentDeg  = [         4,                 1]
% 
% 
%       ColPartToComp = 1;
% 
%         CompColor(1,:) = [ 102, 102, 102];
%         CompColor(2,:) = [ 204,   0,   0];
%         CompColor      = CompColor ./ 256;
% 
%         ColorVec = CompColor
% 
%       RxLxIdx          = [-1, 1];
%       
%     elseif strcmp(Input.SystNameLong, 'N3_NASA')
% 
%       NTint  = length(T0_Vec);
% 
%       Syst.NAtoms    = 3;
%       AtomsName = [    'N',     'N',    'N'];
%       AtomColor = [0 0.3 0; 0 0.3 0; 0 0.3 0];
%       AtomSize  = [    150,     150,     150];
%       AtomMass  = [  25526.04298d0, 25526.04298d0, 25526.04298d0]
% 
% 
%       MoleculesName    = ['N2_NASA'];
%       MoleculedDissEn  = [ 0.0];
%       NMolecules       = size(MoleculesName,1);
%       DegeneracyFactor = (2);
%       MoleculeMu       = [28.0134d-3];
% 
% 
%       AllMoleculesName = ['N2_NASA'; 'N2_NASA'; 'N2_NASA'];
% 
%       PairColor = [0 0 256; 0 256 0; 256 0 0]./256;
% 
%       BinnedMolName   = ['N2_NASA'];
%       NBinnedMol      = size(BinnedMolName,1);
%       BinnedMolToComp = [   2];
% 
%       NComp = 2;
%         CompNames=['  ';'  '];
%         CompNames(1,:) = ' N';
%         CompNames(2,:) = 'N2';
% 
%         ComponentMass = [AtomMass(1), 2.d0.*AtomMass(1)];
%         ComponentDeg  = [         4,                 1]
% 
% 
%       ColPartToComp = 1;
% 
%         CompColor(1,:) = [ 102, 102, 102];
%         CompColor(2,:) = [ 204,   0,   0];
%         CompColor      = CompColor ./ 256;
% 
%         ColorVec = CompColor
% 
%       RxLxIdx          = [-1, 1];
% 
    elseif strcmp(Input.SystNameLong, 'O3_UMN')
        
        %%% System
        Syst.Name  = 'O3';
        
        Syst.NProc = 3; %(Diss+Inel+Exch)
        
        
        %%% Atoms
        Syst.NAtoms        = 3;
        
        Syst.Atom(1).Name  = 'O';
        Syst.Atom(2).Name  = 'O';
        Syst.Atom(3).Name  = 'O';

        Syst.Atom(1).Color = [0, 0, 1];
        Syst.Atom(2).Color = [0, 0, 1];
        Syst.Atom(3).Color = [0, 0, 1];

        Syst.Atom(1).Size  = 200;
        Syst.Atom(2).Size  = 200;
        Syst.Atom(3).Size  = 200;

        Syst.Atom(1).Mass  = 29148.94559;
        Syst.Atom(2).Mass  = 29148.94559;
        Syst.Atom(3).Mass  = 29148.94559;
        
        
        %%% Molecules
        Syst.NMolecules                   = 1;
        
        Syst.Molecule(1).Name             = 'O2';

        Syst.Molecule(1).DissEn           = 5.2112618288711969;

        Syst.Molecule(1).DegeneracyFactor = 3;

        Syst.Molecule(1).Mu               = 31.9988e-3;

        Syst.Molecule(1).NLevelsOrig      = 6115;

        Syst.Molecule(1).ToAtoms          = [1,2];

        Syst.Molecule(1).DiatPot          = 'O2_UMN';

        
        %%% Pairs
        Syst.Pair(1).Name  = 'O2';
        Syst.Pair(2).Name  = 'O2';
        Syst.Pair(3).Name  = 'O2';

        Syst.Pair(1).ToMol = 1;
        Syst.Pair(2).ToMol = 1;
        Syst.Pair(3).ToMol = 1;

        Syst.Pair(1).Color = [0, 0, 256]  ./ 256;
        Syst.Pair(2).Color = [0, 0, 256]  ./ 256;
        Syst.Pair(3).Color = [0, 0, 256]  ./ 256;

        
        %% CFD Components (For PLATO and KONIG)
        Syst.NComp             =  2;
        
        Syst.CFDComp(1).Name   = 'O';
        Syst.CFDComp(2).Name   = 'O2';

        Syst.CFDComp(1).ToMol   = 0;
        Syst.CFDComp(2).ToMol   = 1;

        Syst.CFDComp(1).Mass    = Syst.Atom(1).Mass;
        Syst.CFDComp(2).Mass    = Syst.Atom(1).Mass + Syst.Atom(2).Mass;

        Syst.CFDComp(1).Deg     = 9;
        Syst.CFDComp(2).Deg     = 3;

        Syst.CFDComp(1).Color   = [   0,   0,   0] ./ 256;
        Syst.CFDComp(2).Color   = [   0,   0, 234] ./ 256;

        Syst.CFDComp(1).LineStyle = ':';
        Syst.CFDComp(2).LineStyle = '-';

        Syst.CFDComp(1).RxLxIdx = -1;
        Syst.CFDComp(2).RxLxIdx =  1;
        
        Syst.MolToCFDComp       = [2];

        
        %% Exchange Properties
        Syst.ExchToMol          = [1];
        Syst.ExchToAtom         = [3];

        Syst.PairToExch         = [1; 1];
        
        Syst.ColPartToComp      = 2; 

      
%     elseif strcmp(Input.SystNameLong, 'HNC')
%       
%       NTint  = length(T0_Vec);
% 
%       Syst.NAtoms    = 3;
%       AtomsName = [             'H',            'N',            'C'];
%       AtomColor = [           1 0 0;          0 1 0;          0 0 1];
%       AtomSize  = [             100,            200,            150];
%       AtomMass  = [  1835.0397616d0, 25519.042285d0, 21868.661757d0];
% 
% 
%       MoleculesName    = ['HN'; 'CH'; 'CN'];
%       MoleculedDissEn  = [ 0.0;  0.0; 0.0];
%       NMolecules       = length(MoleculesName);
%       DegeneracyFactor = [3; 2; 2];
%       MoleculeMu       = [15.01454e-3, 13.01854e-3, 26.0174e-3];
% 
%       AllMoleculesName = ['HN'; 'CH'; 'CN'];
% 
%       PairColor = [0 0 256; 0 256 0; 256 0 0]./256;
% 
%       BinnedMolName   = ['HN'; 'CH'; 'CN'];
%       NBinnedMol      = length(BinnedMolName);
%       BinnedMolToComp = [   6;    5;    4];
% 
%       NComp = 6;
%         CompNames=['  ';'  ';'  ';'  ';'  ';'  '];
%         CompNames(1,1:2) = ' C';
%         CompNames(2,1:2) = ' H';
%         CompNames(3,1:2) = ' N';
%         CompNames(4,1:2) = 'CN';
%         CompNames(5,1:2) = 'CH';
%         CompNames(6,1:2) = 'HN';
% 
%         ComponentMass = [ AtomMass(1), AtomMass(2), AtomMass(3), AtomMass(1)+AtomMass(2), AtomMass(1)+AtomMass(3), AtomMass(2)+AtomMass(3)];
%         ComponentDeg  = [           1,           1,           1,                       1,                       1,                       1];
% 
%       ColPartToComp = 4;
% 
%         CompColor(1,:) = [ 102, 102, 102];
%         CompColor(2,:) = [   0, 153, 102];
%         CompColor(3,:) = [ 204,   0,   0];
%         CompColor(4,:) = [   0,   0, 234];
%         CompColor(5,:) = [   0, 204,   0];
%         CompColor(6,:) = [ 100,   0, 100];
%         CompColor      = CompColor ./ 256;
% 
%       ColorVec = CompColor;
% 
%       RxLxIdx          = [-1, -1, -1, 1, 1, 1];
%       
%       Input.SystNameLong(1:3)      = 'CHN'
%     
%       
%     elseif strcmp(Input.SystNameLong, 'CNH')
%       
%       NTint  = length(T0_Vec);
% 
%       Syst.NAtoms    = 3;
%       AtomsName = [             'C',            'N',       'H'];
%       AtomColor = [           0 0 1;          0 1 0;     1 0 0];
%       AtomSize  = [             150,            200,       100];
%       AtomMass  = [  21868.661757d0, 25519.042285d0, 1835.0397616d0];
% 
% 
%       MoleculesName    = ['CN'; 'CH'; 'NH'];
%       MoleculedDissEn  = [ 0.0;  0.0;  0.0];
%       NMolecules       = length(MoleculesName);
%       DegeneracyFactor = [2; 2; 3];
%       MoleculeMu       = [26.0174e-3, 13.01854e-3, 15.01454e-3];
% 
%       AllMoleculesName = ['CN'; 'CH'; 'NH'];
% 
%       PairColor = [256 0 0; 0 256 0; 0 0 256]./256;
% 
%       BinnedMolName   = ['CN'; 'CH'; 'NH'];
%       NBinnedMol      = length(BinnedMolName);
%       BinnedMolToComp = [   4;    5;    6];
% 
%       NComp = 6;
%         CompNames=['  ';'  ';'  ';'  ';'  ';'  '];
%         CompNames(1,1:2) = ' C';
%         CompNames(2,1:2) = ' H';
%         CompNames(3,1:2) = ' N';
%         CompNames(4,1:2) = 'CN';
%         CompNames(5,1:2) = 'CH';
%         CompNames(6,1:2) = 'HN';
% 
%         ComponentMass = [ AtomMass(1), AtomMass(2), AtomMass(3), AtomMass(1)+AtomMass(2), AtomMass(1)+AtomMass(3), AtomMass(2)+AtomMass(3)];
%         ComponentDeg  = [           1,           1,           1,                       1,                       1,                       1];
% 
%       ColPartToComp = 4;
% 
%         CompColor(1,:) = [ 102, 102, 102];
%         CompColor(2,:) = [   0, 153, 102];
%         CompColor(3,:) = [ 204,   0,   0];
%         CompColor(4,:) = [   0,   0, 234];
%         CompColor(5,:) = [   0, 204,   0];
%         CompColor(6,:) = [ 100,   0, 100];
%         CompColor      = CompColor ./ 256;
% 
%       ColorVec = CompColor;
% 
%       RxLxIdx          = [-1, -1, -1, 1, 1, 1];
%       
%       Input.SystNameLong(1:3)      = 'CHN'
%       
%       
%     elseif strcmp(Input.SystNameLong, 'CHN')
%       
%       NTint  = length(T0_Vec);
% 
%       Syst.NAtoms    = 3;
%       AtomsName = [             'C',            'H',           'N'];
%       AtomColor = [           0 0 1;          1 0 0;          0 1 0];
%       AtomSize  = [             150,            100,            200];
%       AtomMass  = [  21868.661757d0, 1835.0397616d0, 25519.042285d0];
% 
% 
%       MoleculesName    = ['CH'; 'CN'; 'HN'];
%       MoleculedDissEn  = [ 0.0;  0.0;  0.0];
%       NMolecules       = length(MoleculesName);
%       DegeneracyFactor = [2; 2; 3];
%       MoleculeMu       = [15.01454e-3, 13.01854e-3, 26.0174e-3];
% 
%       AllMoleculesName = ['CH'; 'CN'; 'HN'];
% 
%       PairColor = [0 0 256; 0 256 0; 256 0 0]./256;
% 
%       BinnedMolName   = ['CH'; 'CN'; 'HN'];
%       NBinnedMol      = length(BinnedMolName);
%       BinnedMolToComp = [   5;    4;    6];
% 
%       NComp = 6;
%         CompNames=['  ';'  ';'  ';'  ';'  ';'  '];
%         CompNames(1,1:2) = ' C';
%         CompNames(2,1:2) = ' H';
%         CompNames(3,1:2) = ' N';
%         CompNames(4,1:2) = 'CN';
%         CompNames(5,1:2) = 'CH';
%         CompNames(6,1:2) = 'HN';
% 
%         ComponentMass = [ AtomMass(1), AtomMass(2), AtomMass(3), AtomMass(1)+AtomMass(2), AtomMass(1)+AtomMass(3), AtomMass(2)+AtomMass(3)];
%         ComponentDeg  = [           1,           1,           1,                       1,                       1,                       1];
% 
%       ColPartToComp = 4;
% 
%         CompColor(1,:) = [ 102, 102, 102];
%         CompColor(2,:) = [   0, 153, 102];
%         CompColor(3,:) = [ 204,   0,   0];
%         CompColor(4,:) = [   0,   0, 234];
%         CompColor(5,:) = [   0, 204,   0];
%         CompColor(6,:) = [ 100,   0, 100];
%         CompColor      = CompColor ./ 256;
% 
%       ColorVec = CompColor;
% 
%       RxLxIdx          = [-1, -1, -1, 1, 1, 1];
%       
%     elseif strcmp(Input.SystNameLong, 'H3')
% 
%       NTint  = length(T0_Vec);
% 
%       Syst.NAtoms    = 3;
%       AtomsName = [    'H',     'H',    'H'];
%       AtomColor = [0 0.3 0; 0 0.3 0; 0 0.3 0];
%       AtomSize  = [    150,     150,     150];
%       AtomMass  = [  1835.0397616d0, 1835.0397616d0, 1835.0397616d0]
% 
% 
%       MoleculesName    = ['H2'];
%       MoleculedDissEn  = [ 0.0];
%       NMolecules       = size(MoleculesName,1);
%       DegeneracyFactor = (2);
%       MoleculeMu       = [2.01588];
% 
% 
%       AllMoleculesName = ['H2'; 'H2'; 'H2'];
% 
%       PairColor = [0 0 256; 0 256 0; 256 0 0]./256;
% 
%       BinnedMolName   = ['H2'];
%       NBinnedMol      = size(BinnedMolName,1);
%       BinnedMolToComp = [   2];
% 
%       NComp = 2;
%         CompNames=['  ';'  '];
%         CompNames(1,:) = ' H';
%         CompNames(2,:) = 'H2';
% 
%         ComponentMass = [AtomMass(1), 2.d0.*AtomMass(1)];
%         ComponentDeg  = [         4,                 1]
% 
% 
%       ColPartToComp = 1;
% 
%         CompColor(1,:) = [ 102, 102, 102];
%         CompColor(2,:) = [ 204,   0,   0];
%         CompColor      = CompColor ./ 256;
% 
%         ColorVec = CompColor
% 
%       RxLxIdx          = [-1, 1];
%       
%       
%    elseif strcmp(Input.SystNameLong, 'NaNbNcNd_NASA')
% 
%       NTint  = length(T0_Vec);
% 
%       Syst.NAtoms    = 4;
%       AtomsName = [    'N',     'N',     'N',      'N'];
%       AtomColor = [0 0.3 0; 0 0.3 0; 0 0.3 0;  0 0.3 0];
%       AtomSize  = [    150,     150,     150,      150];
%       AtomMass  = [  25526.04298d0, 25526.04298d0, 25526.04298d0, 25526.04298d0]
% 
% 
%       MoleculesName    = ['NaNb'];
%       MoleculedDissEn  = [ 0.0];
%       NMolecules       = size(MoleculesName,1);
%       DegeneracyFactor = (2);
%       MoleculeMu       = [28.0134d-3];
% 
% 
%       AllMoleculesName = ['N2_NASA'; 'N2_NASA'; 'N2_NASA'; 'N2_NASA'; 'N2_NASA'; 'N2_NASA'];
% 
%       PairColor = [0 0 256; 0 256 0; 256 0 0]./256;
% 
%       BinnedMolName   = ['N2_NASA'];
%       NBinnedMol      = size(BinnedMolName,1);
%       BinnedMolToComp = [   2];
% 
%       NComp = 2;
%         CompNames=['  ';'  '];
%         CompNames(1,:) = ' N';
%         CompNames(2,:) = 'N2';
% 
%         ComponentMass = [AtomMass(1), 2.d0.*AtomMass(1)];
%         ComponentDeg  = [          4,                 1]
% 
% 
%       ColPartToComp = 1;
% 
%         CompColor(1,:) = [ 102, 102, 102];
%         CompColor(2,:) = [ 204,   0,   0];
%         CompColor      = CompColor ./ 256;
% 
%         ColorVec = CompColor
% 
%       RxLxIdx          = [-1, 1];
% 
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
        
        Syst.Molecule(1).DegeneracyFactor = 1;
        
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

        Syst.CFDComp(1).Deg     = 1;
        Syst.CFDComp(2).Deg     = 1;

        Syst.CFDComp(1).Color   = [ 102, 102, 102] ./ 256;
        Syst.CFDComp(2).Color   = [   0,   0,   0] ./ 256;

        Syst.CFDComp(1).LineStyle = ':';
        Syst.CFDComp(2).LineStyle = '-';
        
        Syst.CFDComp(1).RxLxIdx = -1;
        Syst.CFDComp(2).RxLxIdx =  1;
        
        
        %% Exchange Properties
        Syst.MolToCFDComp       = [2];

        Syst.ExchToMol          = [1,1];
        Syst.ExchToAtom         = [0,0];

        Syst.ColPartToComp      = 1; 
        
   
%     elseif strcmp(Input.SystNameLong, 'N4_UMN')
% 
%       NTint  = length(T0_Vec);
% 
%       Syst.NAtoms    = 4;
%       AtomsName = [    'N',     'N',     'N',      'N'];
%       AtomColor = [0 0.3 0; 0 0.3 0; 0 0.3 0;  0 0.3 0];
%       AtomSize  = [    150,     150,     150,      150];
%       AtomMass  = [  25526.04298d0, 25526.04298d0, 25526.04298d0, 25526.04298d0]
% 
% 
%       MoleculesName    = ['N2'];
%       MoleculedDissEn  = [ 0.0];
%       NMolecules       = size(MoleculesName,1);
%       DegeneracyFactor = (2);
%       MoleculeMu       = [28.0134d-3];
% 
% 
%       AllMoleculesName = ['N2_NASA'; 'N2_NASA'; 'N2_NASA'; 'N2_NASA'; 'N2_NASA'; 'N2_NASA'];
% 
%       PairColor = [0 0 256; 0 256 0; 256 0 0]./256;
% 
%       BinnedMolName   = ['N2'];
%       NBinnedMol      = size(BinnedMolName,1);
%       BinnedMolToComp = [   2];
% 
%       NComp = 2;
%         CompNames=['  ';'  '];
%         CompNames(1,:) = ' N';
%         CompNames(2,:) = 'N2';
% 
%         ComponentMass = [AtomMass(1), 2.d0.*AtomMass(1)];
%         ComponentDeg  = [          4,                 1]
% 
% 
%       ColPartToComp = 1;
% 
%         CompColor(1,:) = [ 102, 102, 102];
%         CompColor(2,:) = [ 204,   0,   0];
%         CompColor      = CompColor ./ 256;
% 
%         ColorVec = CompColor
% 
%       RxLxIdx          = [-1, 1];

    end
    
    
    Syst
    fprintf('====================================================\n\n')  

end