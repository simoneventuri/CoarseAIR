%% Computing Turning Points for the Diatomic Potential
%
function Syst = LevelProperties( Syst, iMol )

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

    global Param
    
    fprintf('\n  = LevelProperties ==================================\n')
    fprintf('  ====================================================\n')
    
    fprintf('  ReComputing Level Properties for Molecule Nb %i \n',  iMol )


    clear( Syst.Molecule(iMol).DiatPot )   
    [Syst.Molecule(iMol).EeVDiss, dVv] = feval(Syst.Molecule(iMol).DiatPot, 100.0);
    
    Syst.Molecule(iMol).LevelEeV  = Syst.Molecule(iMol).LevelEeV - Syst.Molecule(iMol).EeVDiss;
    
    rMin = zeros(Syst.Molecule(iMol).Njqn+1,1);
    VMin = zeros(Syst.Molecule(iMol).Njqn+1,1);
    rMax = zeros(Syst.Molecule(iMol).Njqn+1,1);
    VMax = zeros(Syst.Molecule(iMol).Njqn+1,1);
    rMinOld = 1.2;
    rMaxOld = 100.0;

    for ijqn=0:Syst.Molecule(iMol).Njqn
        if (mod(ijqn,20) == 0)
            fprintf("    Computing Diatomic Potential's Extrema for jqn = %i \n", ijqn)
        end
        [rMin(ijqn+1), VMin(ijqn+1), rMax(ijqn+1), VMax(ijqn+1)] = MinMax(rMinOld, rMaxOld, double(ijqn), Syst.Molecule(iMol).iMol);
        rMinOld = rMin(ijqn+1);
        rMaxOld = rMax(ijqn+1);
    end

    for iLevel=1:Syst.Molecule(iMol).NLevels
        ijqn = Syst.Molecule(iMol).Leveljqn(iLevel);
        if (mod(iLevel,100) == 0)
            fprintf("    Computing Levels' Turning Points for iLevel = %i \n", iLevel)
        end
        Syst.Molecule(iMol).LevelrMin(iLevel) = rMin(ijqn+1);
        Syst.Molecule(iMol).LevelVMin(iLevel) = VMin(ijqn+1);
        Syst.Molecule(iMol).LevelrMax(iLevel) = rMax(ijqn+1);
        Syst.Molecule(iMol).LevelVMax(iLevel) = VMax(ijqn+1);
        %[Syst.Molecule(iMol).LevelrIn(iLevel), Syst.Molecule(iMol).LevelrOut(iLevel)] = TurningPoints( rMin(ijqn+1)-0.001, rMin(ijqn+1)+0.001, Syst.Molecule(iMol).LevelEeV(iLevel), double(ijqn), Syst.Molecule(iMol).iMol);
    end
        
    Syst.Molecule(iMol).EEhDiss   = Syst.Molecule(iMol).EeVDiss  ./ Param.EhToeV;
    Syst.Molecule(iMol).LevelEEh  = Syst.Molecule(iMol).LevelEeV ./ Param.EhToeV;
     
end