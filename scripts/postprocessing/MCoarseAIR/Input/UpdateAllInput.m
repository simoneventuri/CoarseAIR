%% The Function loads the parameters of the chemical system
%
%  Input Global Var: - System: The chemical system of interest (e.g.: N3 / CO2 / O3 / etc.)
%                    - T0_Vec: Vector of Translational Temperatures (e.g.: [10000])
%                       

function UpdateAllInput()    

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
  
  global System T0_Vec
  global NBins NTint NAtoms AtomsName MoleculesName NMolecules DegeneracyFactor BinnedMolName ColPartToComp NBinnedMol BinnedMolToComp NComp ...
         CompNames CompColor AtomColor AtomSize AllMoleculesName PairColor AtomMass ComponentMass ColorVec ComponentDeg MoleculeMu RxLxIdx MoleculedDissEn
  global BCVec
  
    if strcmp(System, 'CO2')

      NTint  = length(T0_Vec);

      NAtoms    = 3;
      AtomsName = [  'C',   'O',   'O'];
      AtomColor = [0 0 0; 0 0 1; 0 0 1];
      AtomSize  = [  100,   200,   200];
      AtomMass  = [  21868.661757d0, 29148.94559d0, 29148.94559d0]

      MoleculesName    = ['CO'; 'O2'];
      MoleculedDissEn  = [ 0.0;  0.0];
      NMolecules       = length(MoleculesName);
      DegeneracyFactor = [1; 6];
      MoleculeMu       = [28.0104e-3, 31.9988d-3];

      AllMoleculesName = ['CO'; 'CO'; 'O2'];

      PairColor = [17 17 17; 17 17 17; 0 0 256]./256;

      BinnedMolName   = ['CO'; 'O2'];
      NBinnedMol      = length(BinnedMolName);
      BinnedMolToComp = [   3;    4];

      NComp = 4;
        CompNames=['  ';'  ';'  ';'  '];
        CompNames(1,1:2) = ' C';
        CompNames(2,1:2) = ' O';
        CompNames(3,1:2) = 'CO';
        CompNames(4,1:2) = 'O2';

        ComponentMass = [ AtomMass(1), AtomMass(2), AtomMass(1)+AtomMass(2), 2.d0.*AtomMass(2)];
        ComponentDeg  = [           1,           9,                       1,                 1]

      ColPartToComp = 2;

        CompColor(1,:) = [ 102, 102, 102];
        CompColor(2,:) = [   0, 153, 102];
        CompColor(3,:) = [ 204,   0,   0];
        CompColor(4,:) = [   0,   0, 234];
        CompColor      = CompColor ./ 256;

      ColorVec = CompColor

      RxLxIdx          = [-1, -1, 1, 0];


    elseif strcmp(System, 'O2C')
      System = 'CO2';

      NTint  = length(T0_Vec);

      NAtoms    = 3;
      AtomsName = [  'O',   'O',   'C'];
      AtomColor = [0 0 1; 0 0 1; 0 0 0];
      AtomSize  = [  200,   200,   100];
      AtomMass  = [  29148.94559d0, 29148.94559d0, 21868.661757d0]


      MoleculesName    = ['O2'; 'CO'];
      MoleculedDissEn  = [ 0.0;  0.0];
      NMolecules       = length(MoleculesName);
      DegeneracyFactor = [6; 1];
      MoleculeMu       = [31.9988d-3, 28.0104e-3];

      AllMoleculesName = ['O2'; 'CO'; 'CO'];

      PairColor = [0 0 256; 17 17 17; 17 17 17]./256;

      BinnedMolName   = ['O2'; 'CO'];
      NBinnedMol      = length(BinnedMolName);
      BinnedMolToComp = [   4;    3];

      NComp = 4;
        CompNames=['  ';'  ';'  ';'  '];
        CompNames(1,1:2) = ' C';
        CompNames(2,1:2) = ' O';
        CompNames(3,1:2) = 'CO';
        CompNames(4,1:2) = 'O2';

        ComponentMass = [ AtomMass(3), AtomMass(1), AtomMass(1)+AtomMass(3), 2.d0.*AtomMass(1)];
        ComponentDeg  = [           1,           9,                       1,                 1];

      ColPartToComp = 2;

        CompColor(1,:) = [ 102, 102, 102];
        CompColor(2,:) = [   0, 153, 102];
        CompColor(3,:) = [ 204,   0,   0];
        CompColor(4,:) = [   0,   0, 234];
        CompColor      = CompColor ./ 256;

        ColorVec = CompColor

      RxLxIdx          = [-1, -1, 0, 1];
      
      
    elseif strcmp(System, 'N2O')

      NTint  = length(T0_Vec);

      NAtoms    = 3;
      AtomsName = [  'N',   'N',   'O'];
      AtomColor = [0 0 0; 0 0 1; 0 0 1];
      AtomSize  = [  200,   200,   200];
      AtomMass  = [  25526.04298d0, 25526.04298d0, 229148.94559d0]

      MoleculesName    = ['N2_UMN'; 'NO_UMN'];
      MoleculedDissEn  = [ 0.0;  0.0];
      NMolecules       = length(MoleculesName);
      DegeneracyFactor = [1; 6];
      MoleculeMu       = [28.0104e-3, 31.9988d-3];

      AllMoleculesName = ['N2_UMN'; 'NO_UMN'; 'NO_UMN'];

      PairColor = [17 17 17;0 0 256; 0 0 256]./256;

      BinnedMolName   = ['N2_UMN'; 'NO_UMN'];
      NBinnedMol      = length(BinnedMolName);
      BinnedMolToComp = [   3;    4];

      NComp = 4;
        CompNames=['  ';'  ';'  ';'  '];
        CompNames(1,1:2) = ' N';
        CompNames(2,1:2) = ' O';
        CompNames(3,1:2) = 'N2';
        CompNames(4,1:2) = 'NO';

        ComponentMass = [ AtomMass(1), AtomMass(2), AtomMass(1)+AtomMass(2), AtomMass(1)+AtomMass(3)];
        ComponentDeg  = [           1,           9,                       1,                 1]

      ColPartToComp = 2;

        CompColor(1,:) = [ 102, 102, 102];
        CompColor(2,:) = [   0, 153, 102];
        CompColor(3,:) = [ 204,   0,   0];
        CompColor(4,:) = [   0,   0, 234];
        CompColor      = CompColor ./ 256;

      ColorVec = CompColor

      RxLxIdx          = [-1, -1, 1, 0];


    elseif strcmp(System, 'N3_UMN')

      NTint  = length(T0_Vec);

      NAtoms    = 3;
      AtomsName = [    'N',     'N',    'N'];
      AtomColor = [0 0.3 0; 0 0.3 0; 0 0.3 0];
      AtomSize  = [    150,     150,     150];
      AtomMass  = [  25526.04298d0, 25526.04298d0, 25526.04298d0]


      MoleculesName    = ['N2_UMN'];
      MoleculedDissEn  = [ 0.0];
      NMolecules       = size(MoleculesName,1);
      DegeneracyFactor = [1];
      MoleculeMu       = [28.0134d-3];


      AllMoleculesName = ['N2_UMN'; 'N2_UMN'; 'N2_UMN'];

      PairColor = [0 0 256; 0 256 0; 256 0 0]./256;

      BinnedMolName   = ['N2_UMN'];
      NBinnedMol      = size(BinnedMolName,1);
      BinnedMolToComp = [   2];

      NComp = 2;
        CompNames=['  ';'  '];
        CompNames(1,:) = ' N';
        CompNames(2,:) = 'N2';

        ComponentMass = [AtomMass(1), 2.d0.*AtomMass(1)];
        ComponentDeg  = [         4,                 1]


      ColPartToComp = 1;

        CompColor(1,:) = [ 102, 102, 102];
        CompColor(2,:) = [ 204,   0,   0];
        CompColor      = CompColor ./ 256;

        ColorVec = CompColor

      RxLxIdx          = [-1, 1];
      
    elseif strcmp(System, 'N3_NASA')

      NTint  = length(T0_Vec);

      NAtoms    = 3;
      AtomsName = [    'N',     'N',    'N'];
      AtomColor = [0 0.3 0; 0 0.3 0; 0 0.3 0];
      AtomSize  = [    150,     150,     150];
      AtomMass  = [  25526.04298d0, 25526.04298d0, 25526.04298d0]


      MoleculesName    = ['N2_NASA'];
      MoleculedDissEn  = [ 0.0];
      NMolecules       = size(MoleculesName,1);
      DegeneracyFactor = [1];
      MoleculeMu       = [28.0134d-3];


      AllMoleculesName = ['N2_NASA'; 'N2_NASA'; 'N2_NASA'];

      PairColor = [0 0 256; 0 256 0; 256 0 0]./256;

      BinnedMolName   = ['N2_NASA'];
      NBinnedMol      = size(BinnedMolName,1);
      BinnedMolToComp = [   2];

      NComp = 2;
        CompNames=['  ';'  '];
        CompNames(1,:) = ' N';
        CompNames(2,:) = 'N2';

        ComponentMass = [AtomMass(1), 2.d0.*AtomMass(1)];
        ComponentDeg  = [         4,                 1]


      ColPartToComp = 1;

        CompColor(1,:) = [ 102, 102, 102];
        CompColor(2,:) = [ 204,   0,   0];
        CompColor      = CompColor ./ 256;

        ColorVec = CompColor

      RxLxIdx          = [-1, 1];

    elseif strcmp(System, 'O3')

      NTint  = length(T0_Vec);

      NAtoms     = 3;
      AtomsName  = [    'O',     'O',    'O'];
      AtomColor  = [0    0.2745    0.7843; 0    0.2745    0.7843; 0    0.2745    0.7843];
      AtomSize   = [    150,     150,     150];
      AtomMass   = [  29148.94559d0, 29148.94559d0, 29148.94559d0]

      MoleculesName    = ['O2'];
      MoleculedDissEn  = [5.2112618288711969];
      NMolecules       = size(MoleculesName,1);
      DegeneracyFactor = [1];
      MoleculeMu       = [31.9988d-3];

      AllMoleculesName = ['O2'; 'O2'; 'O2'];

      PairColor = [0 0 256; 0 256 0; 256 0 0]./256;

      BinnedMolName   = ['O2'];
      NBinnedMol      = size(BinnedMolName,1);
      BinnedMolToComp = [   2];

      NComp = 2;
        CompNames=['  ';'  '];
        CompNames(1,:) = ' O';
        CompNames(2,:) = 'O2';

        ComponentMass = [AtomMass(1), 2.d0.*AtomMass(1)];
        ComponentDeg  = [          9,                 1]


      ColPartToComp = 1;

        CompColor(1,:) = [ 102, 102, 102];
        CompColor(2,:) = [ 204,   0,   0];
        CompColor      = CompColor ./ 256;

        ColorVec = [CompColor; [0,   204,   0]./ 256; [0,   0,   204]./ 256];

      RxLxIdx          = [-1, 1];
      
      
    elseif strcmp(System, 'HNC')
      
      NTint  = length(T0_Vec);

      NAtoms    = 3;
      AtomsName = [             'H',            'N',            'C'];
      AtomColor = [           1 0 0;          0 1 0;          0 0 1];
      AtomSize  = [             100,            200,            150];
      AtomMass  = [  1835.0397616d0, 25519.042285d0, 21868.661757d0];


      MoleculesName    = ['HN'; 'CH'; 'CN'];
      MoleculedDissEn  = [ 0.0;  0.0; 0.0];
      NMolecules       = length(MoleculesName);
      DegeneracyFactor = [3; 2; 2];
      MoleculeMu       = [15.01454e-3, 13.01854e-3, 26.0174e-3];

      AllMoleculesName = ['HN'; 'CH'; 'CN'];

      PairColor = [0 0 256; 0 256 0; 256 0 0]./256;

      BinnedMolName   = ['HN'; 'CH'; 'CN'];
      NBinnedMol      = length(BinnedMolName);
      BinnedMolToComp = [   6;    5;    4];

      NComp = 6;
        CompNames=['  ';'  ';'  ';'  ';'  ';'  '];
        CompNames(1,1:2) = ' C';
        CompNames(2,1:2) = ' H';
        CompNames(3,1:2) = ' N';
        CompNames(4,1:2) = 'CN';
        CompNames(5,1:2) = 'CH';
        CompNames(6,1:2) = 'HN';

        ComponentMass = [ AtomMass(1), AtomMass(2), AtomMass(3), AtomMass(1)+AtomMass(2), AtomMass(1)+AtomMass(3), AtomMass(2)+AtomMass(3)];
        ComponentDeg  = [           1,           1,           1,                       1,                       1,                       1];

      ColPartToComp = 4;

        CompColor(1,:) = [ 102, 102, 102];
        CompColor(2,:) = [   0, 153, 102];
        CompColor(3,:) = [ 204,   0,   0];
        CompColor(4,:) = [   0,   0, 234];
        CompColor(5,:) = [   0, 204,   0];
        CompColor(6,:) = [ 100,   0, 100];
        CompColor      = CompColor ./ 256;

      ColorVec = CompColor;

      RxLxIdx          = [-1, -1, -1, 1, 1, 1];
      
      System(1:3)      = 'CHN'
    
      
    elseif strcmp(System, 'CNH')
      
      NTint  = length(T0_Vec);

      NAtoms    = 3;
      AtomsName = [             'C',            'N',       'H'];
      AtomColor = [           0 0 1;          0 1 0;     1 0 0];
      AtomSize  = [             150,            200,       100];
      AtomMass  = [  21868.661757d0, 25519.042285d0, 1835.0397616d0];


      MoleculesName    = ['CN'; 'CH'; 'NH'];
      MoleculedDissEn  = [ 0.0;  0.0;  0.0];
      NMolecules       = length(MoleculesName);
      DegeneracyFactor = [2; 2; 3];
      MoleculeMu       = [26.0174e-3, 13.01854e-3, 15.01454e-3];

      AllMoleculesName = ['CN'; 'CH'; 'NH'];

      PairColor = [256 0 0; 0 256 0; 0 0 256]./256;

      BinnedMolName   = ['CN'; 'CH'; 'NH'];
      NBinnedMol      = length(BinnedMolName);
      BinnedMolToComp = [   4;    5;    6];

      NComp = 6;
        CompNames=['  ';'  ';'  ';'  ';'  ';'  '];
        CompNames(1,1:2) = ' C';
        CompNames(2,1:2) = ' H';
        CompNames(3,1:2) = ' N';
        CompNames(4,1:2) = 'CN';
        CompNames(5,1:2) = 'CH';
        CompNames(6,1:2) = 'HN';

        ComponentMass = [ AtomMass(1), AtomMass(2), AtomMass(3), AtomMass(1)+AtomMass(2), AtomMass(1)+AtomMass(3), AtomMass(2)+AtomMass(3)];
        ComponentDeg  = [           1,           1,           1,                       1,                       1,                       1];

      ColPartToComp = 4;

        CompColor(1,:) = [ 102, 102, 102];
        CompColor(2,:) = [   0, 153, 102];
        CompColor(3,:) = [ 204,   0,   0];
        CompColor(4,:) = [   0,   0, 234];
        CompColor(5,:) = [   0, 204,   0];
        CompColor(6,:) = [ 100,   0, 100];
        CompColor      = CompColor ./ 256;

      ColorVec = CompColor;

      RxLxIdx          = [-1, -1, -1, 1, 1, 1];
      
      System(1:3)      = 'CHN'
      
      
    elseif strcmp(System, 'CHN')
      
      NTint  = length(T0_Vec);

      NAtoms    = 3;
      AtomsName = [             'C',            'H',           'N'];
      AtomColor = [           0 0 1;          1 0 0;          0 1 0];
      AtomSize  = [             150,            100,            200];
      AtomMass  = [  21868.661757d0, 1835.0397616d0, 25519.042285d0];


      MoleculesName    = ['CH'; 'CN'; 'HN'];
      MoleculedDissEn  = [ 0.0;  0.0;  0.0];
      NMolecules       = length(MoleculesName);
      DegeneracyFactor = [2; 2; 3];
      MoleculeMu       = [15.01454e-3, 13.01854e-3, 26.0174e-3];

      AllMoleculesName = ['CH'; 'CN'; 'HN'];

      PairColor = [0 0 256; 0 256 0; 256 0 0]./256;

      BinnedMolName   = ['CH'; 'CN'; 'HN'];
      NBinnedMol      = length(BinnedMolName);
      BinnedMolToComp = [   5;    4;    6];

      NComp = 6;
        CompNames=['  ';'  ';'  ';'  ';'  ';'  '];
        CompNames(1,1:2) = ' C';
        CompNames(2,1:2) = ' H';
        CompNames(3,1:2) = ' N';
        CompNames(4,1:2) = 'CN';
        CompNames(5,1:2) = 'CH';
        CompNames(6,1:2) = 'HN';

        ComponentMass = [ AtomMass(1), AtomMass(2), AtomMass(3), AtomMass(1)+AtomMass(2), AtomMass(1)+AtomMass(3), AtomMass(2)+AtomMass(3)];
        ComponentDeg  = [           1,           1,           1,                       1,                       1,                       1];

      ColPartToComp = 4;

        CompColor(1,:) = [ 102, 102, 102];
        CompColor(2,:) = [   0, 153, 102];
        CompColor(3,:) = [ 204,   0,   0];
        CompColor(4,:) = [   0,   0, 234];
        CompColor(5,:) = [   0, 204,   0];
        CompColor(6,:) = [ 100,   0, 100];
        CompColor      = CompColor ./ 256;

      ColorVec = CompColor;

      RxLxIdx          = [-1, -1, -1, 1, 1, 1];
      
    elseif strcmp(System, 'H3')

      NTint  = length(T0_Vec);

      NAtoms    = 3;
      AtomsName = [    'H',     'H',    'H'];
      AtomColor = [0 0.3 0; 0 0.3 0; 0 0.3 0];
      AtomSize  = [    150,     150,     150];
      AtomMass  = [  1835.0397616d0, 1835.0397616d0, 1835.0397616d0]


      MoleculesName    = ['H2'];
      MoleculedDissEn  = [ 0.0];
      NMolecules       = size(MoleculesName,1);
      DegeneracyFactor = [1];
      MoleculeMu       = [2.01588];


      AllMoleculesName = ['H2'; 'H2'; 'H2'];

      PairColor = [0 0 256; 0 256 0; 256 0 0]./256;

      BinnedMolName   = ['H2'];
      NBinnedMol      = size(BinnedMolName,1);
      BinnedMolToComp = [   2];

      NComp = 2;
        CompNames=['  ';'  '];
        CompNames(1,:) = ' H';
        CompNames(2,:) = 'H2';

        ComponentMass = [AtomMass(1), 2.d0.*AtomMass(1)];
        ComponentDeg  = [         4,                 1]


      ColPartToComp = 1;

        CompColor(1,:) = [ 102, 102, 102];
        CompColor(2,:) = [ 204,   0,   0];
        CompColor      = CompColor ./ 256;

        ColorVec = CompColor

      RxLxIdx          = [-1, 1];
      

    end

end