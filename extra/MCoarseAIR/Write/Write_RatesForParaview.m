%% The Function plots the Ro-Vibrational Populations at Given Time Steps
%
function Write_RatesForParaview(Controls)    
    
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

    global Input Kin Param Syst Temp Rates

    fprintf('= Write_RatesForParaview =============== T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    
    
    for iMol = 1:1
        fprintf(['Molecule Nb ' num2str(iMol) ', ' Syst.Molecule(iMol).Name '\n'] );

        [status,msg,msgID] = mkdir(Input.Paths.SaveDataFldr);
        FileName           = strcat(Input.Paths.SaveDataFldr, '/DiatPot_and_DissRates.csv');
        fileID = fopen(FileName,'w');
        fprintf(fileID,'id,v,J,EeV,rIn,rOut,EeVVib,EeVRot,dCentBarr,KDiss,KDiss_iP1,KDiss_iP2,KDiss_iP3\n');
        for iLevel = 1:Syst.Molecule(iMol).NLevels
            fprintf(fileID,'%i,%i,%i,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n', iLevel, 	                                 ...
                                                                       Syst.Molecule(iMol).Levelvqn(iLevel),     ...
                                                                       Syst.Molecule(iMol).Leveljqn(iLevel),     ...
                                                                       Syst.Molecule(iMol).LevelEeV(iLevel),     ...
                                                                       Syst.Molecule(iMol).LevelrIn(iLevel),     ...
                                                                       Syst.Molecule(iMol).LevelrOut(iLevel),    ...
                                                                       Syst.Molecule(iMol).LevelEeVVib0(iLevel), ...
                                                                       Syst.Molecule(iMol).LevelEeVRot(iLevel),  ...
                                                                       Syst.Molecule(iMol).LevelECB(iLevel),     ...
                                                                       Rates.T(Temp.iT).Diss(iLevel,1),          ...
                                                                       Rates.T(Temp.iT).Diss(iLevel,2),          ...
                                                                       Rates.T(Temp.iT).Diss(iLevel,3),          ...
                                                                       Rates.T(Temp.iT).Diss(iLevel,4)           ...
                                                                       );
        end
        fclose(fileID);
        
    end
    
    
    fprintf('====================================================\n\n')
    
end