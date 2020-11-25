%% The Function Initializes the Reaming Global Variables 
%
%  Required Variables: - Syst.NMolecules
%
function Initialize_Input()

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

    
    global Input Syst Temp OtherSyst

    
    fprintf('= Initialize_Input =================================\n')
    fprintf('====================================================\n')
    fprintf('Initializing Syst and Temp Objects based on Input Object\n' )
    fprintf('====================================================\n\n')  
    
    
    Syst.NameLong   = Input.SystNameLong;
    
    Syst.ThCollPart = false;
    if exist('Input.Kin.ThCollPart', 'var')
        Syst.ThCollPart = Input.Kin.ThCollPart;
    end
    
    for iMol = 1:size(Input.Kin.MolResolutionIn,2)
        Syst.Molecule(iMol).KinMthdIn    = Input.Kin.MolResolutionIn(iMol);
        Syst.Molecule(iMol).EqNStatesIn  = Input.Kin.EqNStatesIn(iMol);
        
        Syst.Molecule(iMol).MinStateIn   = Input.Kin.MinStateIn(iMol);
        Syst.Molecule(iMol).MaxStateIn   = Input.Kin.MaxStateIn(iMol);
        
        Syst.Molecule(iMol).NGroupsIn    = Input.Kin.NGroupsIn(iMol);
        
        clear( Syst.Molecule(iMol).DiatPot );
        Syst.Molecule(iMol).DissEn       = feval(Syst.Molecule(iMol).DiatPot, fminsearch(Syst.Molecule(iMol).DiatPot, 1.5));
    end
    
    
    Temp.TranVec = Input.TranVec;
    Temp.NTran   = length(Temp.TranVec);

    Temp.IntVec  = Input.TranVec;    
    Temp.NInt    = length(Temp.TranVec);
    
    
    TempStr = '';
    if (Input.Kin.NBinsSuffix > 0)
        TempStr = strcat('_', num2str(Input.Kin.NBinsSuffix), 'Bins');
    end
    Input.Kin.Proc.OverallFlg = strcat(num2str(Input.Kin.Proc.DissFlg),'_',num2str(Input.Kin.Proc.InelFlg),'_',num2str(Input.Kin.Proc.ExchFlg1),'_',num2str(Input.Kin.Proc.ExchFlg2), TempStr);
    
    Syst.iPES = ''
    if (Input.iPES > 0)
        Syst.iPES = strcat('_PES', num2str(Input.iPES));
    end
    Syst.Suffix              = Input.Suffix;

    Input.Paths.SaveFigsFldr = strcat(Input.Paths.SaveFigsFldr, '/', Syst.NameLong, Input.RunSuffix, Syst.iPES, '/');
    Input.Paths.SaveDataFldr = strcat(Input.Paths.SaveDataFldr, '/', Syst.NameLong, Input.RunSuffix, Syst.iPES, '/');

    Syst.HDF5_File           = strcat(Input.Paths.ToHDF5Fldr, Syst.NameLong, Input.Suffix, Syst.iPES, '.hdf5');
    

    
    for iSyst = 1:size(Syst.OtherSyst_NameLong,1)
        OtherSyst(iSyst).Syst.NameLong  = Syst.OtherSyst_NameLong;
        OtherSyst(iSyst).Syst           = Initialize_ChemicalSyst(OtherSyst(iSyst).Syst);
        OtherSyst(iSyst).Syst.HDF5_File = strcat(Input.Paths.ToHDF5Fldr, OtherSyst(iSyst).Syst.NameLong, Input.Suffix, Syst.iPES, '.hdf5');
    
        for iMol = 1:OtherSyst(iSyst).Syst.NMolecules
            OtherSyst(iSyst).Syst.Molecule(iMol).KinMthdIn    = {'StS'};
            OtherSyst(iSyst).Syst.Molecule(iMol).MinStateIn   = 1;
            OtherSyst(iSyst).Syst.Molecule(iMol).MaxStateIn   = 100000;
            OtherSyst(iSyst).Syst.Molecule(iMol).NGroupsIn    = 1;
        end
    end

  
end