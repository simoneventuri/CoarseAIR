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

    
    global Input Syst Temp

    
    fprintf('= Initialize_Input =================================\n')
    fprintf('====================================================\n')
    fprintf('Initializing Syst and Temp Objects based on Input Object\n' )
    fprintf('====================================================\n\n')  
    
    
    Syst.NameLong   = Input.SystNameLong;
    
    for iMol = 1:size(Input.Kin.MolResolutionIn,1)
        Syst.Molecule(iMol).KinMthdIn(:) = Input.Kin.MolResolutionIn(iMol,:);

        Syst.Molecule(iMol).MinStateIn   = Input.Kin.MinStateIn(iMol);
        Syst.Molecule(iMol).MaxStateIn   = Input.Kin.MaxStateIn(iMol);
        
        Syst.Molecule(iMol).NGroupsIn    = Input.Kin.NGroupsIn(iMol);
    end
    
    
    Temp.TranVec = Input.TranVec;
    Temp.NTran   = length(Temp.TranVec);

    Temp.IntVec  = Input.TranVec;    
    Temp.NInt    = length(Temp.TranVec);
    
    TempStr = '';
    if (Input.Kin.NBinsSuffix > 0)
        TempStr = strcat('_', num2str(Input.Kin.NBinsSuffix), 'Bins');
    end
    if (Syst.NAtoms == 3)
        Input.Kin.Proc.OverallFlg = strcat(num2str(Input.Kin.Proc.DissFlg),'_',num2str(Input.Kin.Proc.InelFlg),'_',num2str(Input.Kin.Proc.ExchFlg1),'_',num2str(Input.Kin.Proc.ExchFlg2), TempStr);
    else
        Input.Kin.Proc.OverallFlg = strcat(num2str(Input.Kin.Proc.DissFlg),'_',num2str(Input.Kin.Proc.DissInelFlg),'_',num2str(Input.Kin.Proc.InelFlg),'_',num2str(Input.Kin.Proc.ExchFlg1),'_',num2str(Input.Kin.Proc.ExchFlg2), TempStr);
    end
    
    
    Syst.iPES = ''
    if (Input.iPES > 0)
        Syst.iPES = strcat('_PES', num2str(Input.iPES));
    end
    Syst.Suffix              = Input.Suffix;
    Syst.HDF5_File           = strcat(Input.Paths.ToHDF5Fldr, Syst.NameLong,          Input.Suffix, Syst.iPES, '.hdf5');
    Syst.HDF5_File_OtherExch = strcat(Input.Paths.ToHDF5Fldr, Syst.NameLong_Opposite, Input.Suffix, Syst.iPES, '.hdf5');
    
    Input.Paths.SaveFigsFldr = strcat(Input.Paths.SaveFigsFldr, '/', Syst.NameLong, Input.Suffix, Syst.iPES, '/');

    Input.Paths.SaveDataFldr = strcat(Input.Paths.SaveDataFldr, '/', Syst.NameLong, Input.Suffix, Syst.iPES, '/');

    
%     filename = strcat(Input.Paths.ToQCTFldr,'/InputForBash.inp');
%     delimiter = ' ';
%     formatSpec = '%s%[^\n\r]';
%     fileID = fopen(filename,'r');
%     dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
%     fclose(fileID);
%     InputForBash = [dataArray{1:end-1}];
%     clearvars filename delimiter formatSpec fileID dataArray ans;
% 
%     iTot = 1;
%     if exist('Input.SystNameLong','var')
%         Syst.Name = Input.SystNameLong
%     else
%         Syst.Name = char(InputForBash(iTot,:))
%     end
%     
%     iTot = iTot + 1;
%     Temp.TranFlg = str2num(char(InputForBash(iTot,:)))
% 
%     iTot = iTot + 1;
%     Temp.NTran  = str2num(char(InputForBash(iTot,:)))
%     for iTtra = 1: Temp.NTran
%     iTot = iTot + 1;
%     Temp.TranVec(iTtra) = str2num(char(InputForBash(iTot,:)))
%     end
% 
%     iTot = iTot + 1;
%     Temp.NInt = str2num(char(InputForBash(iTot,:)))
%     if Temp.NInt > 1 
%         for iTtra = 1:Temp.NInt
%             iTot = iTot + 1;
%             Temp.IntVec(iTint) = str2num(char(InputForBash(iTot,:)))
%         end
%     end
% 
%     iTot = iTot + 1;
%     Syst.NMolecules = str2num(char(InputForBash(iTot,:)))
%     for iMol = 1:Syst.NMolecules
%         iTot = iTot + 1;
%         Syst.Molecule(iMol).Name(:) = char(InputForBash(iTot,:))
%         iTot = iTot + 1;
%         SortMthdTemp          = char(InputForBash(iTot,:))
%         if SortMthdTemp(1:3) == 'Sta'
%             Syst.Molecule(iMol).MolResolutionIn(:) = 'STS'
%             iTot = iTot + 1;
%             Syst.Molecule(iMol).MinStateIn = str2num(char(InputForBash(iTot,:)))
%             iTot = iTot + 1;
%             Syst.Molecule(iMol).MaxStateIn = str2num(char(InputForBash(iTot,:)))
%             iTot = iTot + 1;
%         elseif SortMthdTemp(1:3) == 'Vib'
%             Syst.Molecule(iMol).MolResolutionIn(:) = 'VSM'
%             iTot = iTot + 1;
%             Syst.Molecule(iMol).MinStateIn = str2num(char(InputForBash(iTot,:)))
%             iTot = iTot + 1;
%             Syst.Molecule(iMol).MaxStateIn = str2num(char(InputForBash(iTot,:)))
%             iTot = iTot + 1;
%         elseif sum(SortMthdTemp(1:3) == 'RoV') == 3 || sum(SortMthdTemp(1:3) == 'Fro') == 3
%             Syst.Molecule(iMol).MolResolutionIn(:) = 'CGM'
%             iTot = iTot + 1;
%             Syst.Molecule(iMol).MinStateIn = str2num(char(InputForBash(iTot,:)))
%             iTot = iTot + 1;
%             Syst.Molecule(iMol).MaxStateIn = str2num(char(InputForBash(iTot,:)))
%             iTot = iTot + 1;
%             Syst.Molecule(iMol).NGroupsIn  = str2num(char(InputForBash(iTot,:)))
%         end
%     end
  
end
