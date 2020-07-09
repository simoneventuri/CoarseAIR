%% NO Diatomic Potential from UMN (..., 2016) for N2O2 PES
%
function [V, dV] = N2_UMN_ForN2O2_Fun(R)    

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
    
% N2 UMN Min @ 2.075 (V=-9.904361)

    cs   = [ 2.71405774451e0, 1.32757649829e-1, 2.66756890408e-1, 1.95350725241e-1, -4.08663480982e-1, 3.92451705557e-1, 1.13006674877e0]; 
    red  = 1.098d0;
    red4 = 1.4534810048d0;
    de   = 228.4d0;
    VRef = 0.0d0;%0.191619504727d0;


    RAng    = R * Param.BToAng;

    RAng3   = RAng^3;
    RAng4   = RAng3*RAng;
    TempSum = (RAng4 + red4);
    y       = (RAng4 - red4) / TempSum;
    y2      = y^2;
    y3      = y2 * y;
    y4      = y2 * y2;
    y5      = y3 * y2;
    y6      = y3 * y3;

    fy      =   cs(1) + cs(2)*y + cs(3)*y2 + cs(4)*y3 + cs(5)*y4 + cs(6)*y5 + cs(7)*y6;
    u       =   exp(-fy * (RAng-red));
    minu    =   1.d0 - u;

    dfdy    =   cs(2) + 2.0d0*cs(3)*y + 3.0d0*cs(4)*y2 + 4.0d0*cs(5)*y3 + 5.0d0*cs(6)*y4 + 6.0d0*cs(7)*y5;

    dydr    =   8.0d0 * RAng3 * red4 / TempSum^2;
    dfdr    =   dfdy * dydr;

    V       =   de * minu^2 - de;
    dV      =   2.0d0 * de * minu * u * (dfdr * (RAng-red) + fy);


    V  = (V' * Param.KcmToEh + VRef) * Param.EhToeV;
    dV = dV' * Param.KcmToEh * Param.EhToeV * Param.BToAng;

end