%% The Function reads the Rates from the HD5 File
%
%  Input Global Var: - Temp.TNowChar
%                    - Syst.HDF5_File
%
function Read_Rates_FromHDF5()    

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
    
    global Rates Syst Temp
    
    fprintf('  = Read_Rates_FromHDF5 ================== T = %i K\n', Temp.TNow)
    fprintf('  ====================================================\n')
    fprintf('  Reading Rates in HDF5 Format \n' )
    fprintf(['  Reading from File: ' Syst.HDF5_File '\n'] )

    
    if (Syst.NAtoms == 3)
       
        DissChar       = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/Diss/');
        %h5disp(Syst.HDF5_File, DissChar)
        RatesTemp                 = h5read(Syst.HDF5_File, DissChar);
        Rates.T(Temp.iT).Diss     = permute(RatesTemp, [2,1]);
        fprintf(['  Rates.T(' num2str(Temp.iT) ').Diss, size: (' num2str(size(Rates.T(Temp.iT).Diss)) ') \n'])
    
        InelChar       = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/Inel/');
        %h5disp(Syst.HDF5_File, InelChar)
        RatesTemp                 = h5read(Syst.HDF5_File, InelChar);
        Rates.T(Temp.iT).Inel     = permute(RatesTemp, [2,1]);
        fprintf(['  Rates.T(' num2str(Temp.iT) ').Inel, size: (' num2str(size(Rates.T(Temp.iT).Inel)) ') \n'])
    
        for iExch = 1:Syst.NProc-2
            ExchChar       = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/Exch_', num2str(iExch), '/');
            %h5disp(Syst.HDF5_File, ExchChar)
            RatesTemp                         = h5read(Syst.HDF5_File, ExchChar);
            Rates.T(Temp.iT).ExchType(iExch).Exch = permute(RatesTemp, [2,1]);
            fprintf(['  Rates.T(' num2str(Temp.iT) ').ExchType(' num2str(iExch) ').Exch, size: (' num2str(size(Rates.T(Temp.iT).ExchType(iExch).Exch)) ') \n'])
        end
        
    else
    
        DissChar       = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/Diss/');
        %h5disp(Syst.HDF5_File, DissChar)
        RatesTemp                 = h5read(Syst.HDF5_File, DissChar);
        Rates.T(Temp.iT).Diss     = permute(RatesTemp, [3,2,1]);
        fprintf(['  Rates.T(' num2str(Temp.iT) ').Diss, size: (' num2str(size(Rates.T(Temp.iT).Diss)) ') \n'])

        DissCharInel   = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/DissInel/');
        %h5disp(Syst.HDF5_File, DissCharInel)
        RatesTemp                 = h5read(Syst.HDF5_File, DissCharInel);
        Rates.T(Temp.iT).DissInel = permute(RatesTemp, [4,3,2,1]);
        fprintf(['  Rates.T(' num2str(Temp.iT) ').DissInel, size: (' num2str(size(Rates.T(Temp.iT).DissInel)) ') \n'])

    end
    
    
    fprintf('  ====================================================\n\n')
    
end