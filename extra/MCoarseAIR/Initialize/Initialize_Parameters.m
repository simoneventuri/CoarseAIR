%% The Function loads the physical constants and the parameters for plotting
%
%  Input Global Var: - Syst.NAtoms
%                    - Input.FigureFormat ( 'PrePrint' / 'RePrint' / 'Presentation' )
%
function Initialize_Parameters

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
  
    global Input Syst Param
    
    
    fprintf('= Initialize_Parameters ============================\n')
    fprintf('====================================================\n')
    fprintf('Initializing Param Object\n' )
    fprintf('====================================================\n\n')      

  
    Param.Plnck    = 6.62607004d-34;
    Param.UKb      = 1.380658e-23;
    Param.Ue       = 1.602191e-19;
    Param.KJK      = 1.380649e-23
    Param.KeV      = 8.617330e-05;
    Param.AvN      = 6.0221409e+23;
    Param.AMUToKg  = 1.d0/Param.AvN*1.d-3;
    Param.EhToeV   = 27.2113839712790;
    Param.DSWtoKg  = 1.d-3/1.8208e+03;
    Param.ATMToPa  = 101325.d0;

    if Syst.NAtoms == 3
        Param.Pair_to_Atoms = [1,2;1,3;2,3];
    else
        Param.Pair_to_Atoms = [1,2;1,3;1,4;2,3;2,4;3,4];
        Param.iPOpp         = [6,5,4,3,2,1];
    end
    
    Param.iFigure  = 1;
    Param.SaveFigs = 0;
  
    if strcmp(Input.FigureFormat, 'PrePrint')
        Param.LineWidth  = 2;

        Param.AxisFontSz = 26;
        Param.AxisFontNm = 'Times';

        Param.AxisLabelSz = 30;
        Param.AxisLabelNm = 'Times';
        
        Param.LegendFontSz = 24;
        Param.LegendFontNm = 'Times';
        
    elseif strcmp(Input.FigureFormat, 'RePrint')
        Param.LineWidth  = 3;

        Param.AxisFontSz = 36;
        Param.AxisFontNm = 'Times';

        Param.AxisLabelSz = 40;
        Param.AxisLabelNm = 'Times';
        
        Param.LegendFontSz = 32;
        Param.LegendFontNm = 'Times';
        
    elseif strcmp(Input.FigureFormat, 'Presentation')
        Param.LineWidth  = 4;
        
        Param.AxisFontSz = 30;
        Param.AxisFontNm = 'Times';

        Param.AxisLabelSz = 34;
        Param.AxisLabelNm = 'Times';
        
        Param.LegendFontSz = 30;
        Param.LegendFontNm = 'Times';
        
    end

    Param.KCVec = [  0   0   0] ./ 255;
    Param.RCVec = [235, 70, 50] ./ 255;
    Param.BCVec = [55, 80, 165] ./ 255;
    Param.PCVec = [190 140 140] ./ 255;
    Param.GCVec = [65, 140, 70] ./ 255;
    Param.MCVec = [105  65 155] ./ 255;
    Param.YCVec = [245 165  50] ./ 255;
    Param.OCVec = [255 105  45] ./ 255;
    Param.JCVec = [100 100 100] ./ 255;
    Param.CCVec = [205 205 205] ./ 255;
    
    Param.CMat  = [Param.KCVec; Param.RCVec; Param.BCVec; 
                   Param.PCVec; Param.GCVec; Param.MCVec; 
                   Param.YCVec; Param.OCVec; Param.JCVec; 
                   Param.CCVec]; 

    Param.linS  = {'-',':','-.','--'};

end