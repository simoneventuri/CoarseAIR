%% NO Diatomic Potential from UMN (..., 2016)
%
function [V, dV] = NO_UMN_Fun(R)    

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
    
    
% NO UMN Min @ 2.17661 (V=-6.617426)

    cs   = [ 0.322338e0, 5.878590e0, -12.790761e0, 13.320811e0, -7.516309e0, 1.875839e0, -0.052723e0, -0.037783e0, 0.48294e0, 1.98697e0]; 
    red  = 1.1508;
    de   = 152.6;
    VRef = 0.0;%0.191619504727d0;

    RAng = R * Param.BToAng;

    u    = exp( -(RAng-red) / cs(9) - (RAng-red) ^2.0 / (cs(10)) );

    dfdr =  ( -2.0 * (RAng-red)/cs(10) - 1.0/cs(9) );

    V    =  - de*(cs(1) * u + cs(2) * u^2 + cs(3) * u^3 + cs(4) * u^4 + cs(5) * u^5.0 + cs(6) * u^6 + cs(7) * u^7 + cs(8) * u^8);
    V    = (V' * Param.KcmToEh + VRef) * Param.EhToeV;

    dV   =  - de * (cs(1) *          dfdr * u        + ...
                   cs(2) * 2.0d0 * dfdr * u^2     + ...
                   cs(3) * 3.0d0 * dfdr * u^3     + ...
                   cs(4) * 4.0d0 * dfdr * u^4     + ...
                   cs(5) * 5.0d0 * dfdr * u^5     + ...
                   cs(6) * 6.0d0 * dfdr * u^6     + ...
                   cs(7) * 7.0d0 * dfdr * u^7     + ...
                   cs(8) * 8.0d0 * dfdr * u^8);

    dV = dV' * Param.KcmToEh * Param.EhToeV * Param.BToAng;

end