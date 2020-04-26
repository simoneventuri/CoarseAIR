%% The Function reads the Thermodynamics Variables Outputted by KONIG
%
%  Input Global Var: - Syst.NAtoms
%                    - Input.FigureFormat ( 'PrePrint' / 'RePrint' / 'Presentation' )
%
function Read_KONIGBox()    

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
    
    global Input Syst Temp Kin
  
    
    %% Reading Simulation Time
    filename = strcat(Input.Paths.ToKinRunFldr,'/box.dat')
    startRow = 1;
    formatSpec = '%14f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    Kin.T(Temp.iT).t = dataArray{:, 1};
    clearvars startRow formatSpec fileID dataArray ans;


    %% Reading Mole Fractions
    startRow = 1;
    formatSpec = '%*14s';
    for iComp=1:Syst.NComp
        formatSpec=strcat(formatSpec,'%14f');
    end
    formatSpec=strcat(formatSpec,'%[^\n\r]');
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    Kin.T(Temp.iT).MolFracs = [dataArray{1:end-1}];

    
    %% Reading Thermodynamic Variables
    startRow = 1;
    formatSpec = '%*14s';
    for iComp=1:Syst.NComp
        formatSpec=strcat(formatSpec,'%*14f');
    end
    formatSpec=strcat(formatSpec,'%14f%14f%14f%14f%f%[^\n\r]');
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    Kin.T(Temp.iT).TTrans = dataArray{:, 1};
    Kin.T(Temp.iT).rho    = dataArray{:, 2};
    Kin.T(Temp.iT).P      = dataArray{:, 3};
    Kin.T(Temp.iT).nd     = dataArray{:, 4};
    Kin.T(Temp.iT).E      = dataArray{:, 5};


end