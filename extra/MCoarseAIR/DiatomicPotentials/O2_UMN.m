%% O2 Diatomic Potential from UMN (Varga et al., 2018)
%
function [VVec, dVVec] = O2_UMN(RVec)    

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
    
    
    VRef       = 0.0;%0.1915103559;

    a          = [-1.488979427684798e3, 1.881435846488955e4, -1.053475425838226e5, 2.755135591229064e5, -4.277588997761775e5, 4.404104009614092e5, -2.946204062950765e5, 1.176861219078620e5 ];   
    alpha      = 9.439784362354936e-1;
    beta       = 1.262242998506810;

    r2r4Scalar = 2.59361680;                                                                                 
    rs6        = 0.5299;
    rs8        = 2.20;
    c6         = 12.8;

    VVec  = [];
    dVVec = [];
    for iSample = 1:length(RVec)

    %######################################################################################################     <--- Compute_Vd_dVd_O2
    %### 
    R    = RVec(iSample);
    RAng = R * Param.BToAng;


    %##############################################################################################     <--- Ev2gm2_Grad
    %###
    V = 0.0;
    for k=1:8
        V = V + a(k) * exp(-alpha * beta^(k-1) * RAng^2);
    end
    V = V*1.e-3;

    dV = 0.0;
    for k=1:8
        dV = dV - 2.0 * a(k) * alpha * beta^(k-1) * RAng * exp(-alpha * beta^(k-1) * RAng^2);
    end
    dV = dV*1.e-3;
    %######################################################################################     <--- d3disp_Grad
    %###


    %##############################################################################       <--- edisp_Grad
    %###
    c8Step = 3.0 * c6 * r2r4Scalar^2;

    tmp = sqrt(c8Step / c6); 
    e6  = c6     / (R^6 + (rs6*tmp + rs8)^6);
    e8  = c8Step / (R^8 + (rs6*tmp + rs8)^8);

    e6dr =     c6 * (-6.0 * R^5) / (R^6 + (rs6*tmp + rs8)^6)^2;
    e8dr = c8Step * (-8.0 * R^7) / (R^8 + (rs6*tmp + rs8)^8)^2;
    %##############################################################################       ---> edisp_Grad


    VDisp  = (-e6   -2.0*e8);
    dVDisp = (-e6dr -2.0*e8dr) / Param.BToAng;
    %######################################################################################     ---> d3disp_Grad


    VDiat  = VDisp  + V;  
    dVDiat = dVDisp + dV;  
    %##############################################################################################     ---> Ev2gm2_Grad


    V  = (VDiat + VRef)          * Param.EhToeV;   
    dV = (dVDiat * Param.BToAng) * Param.EhToeV;                                                                                  
    %######################################################################################################     ---> Compute_Vd_dVd_O2

    VVec  = [VVec;   V];
    dVVec = [dVVec; dV];

    end      

end