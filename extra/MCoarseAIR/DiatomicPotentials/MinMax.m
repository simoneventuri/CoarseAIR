%% Finding the Minimum and the Maximum of the Diatomic Potential for a given jqn
%     
function [rMin, VMin, rMax, VMax] = MinMax(rMinOld, rMaxOld, jqn, iMol)

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

    global Syst Param

    jqnreal = jqn+0.0;
    
    rMin = fzero(@(x) DiatPotdV(x, jqnreal, iMol), rMinOld);
    if isnan(rMin)
        rMin = 100.0;
    end
    VMin = DiatPot(rMin, jqnreal, iMol);
   
    if (jqn > 0)
        rMax = fzero(@(x) DiatPotdV(x, jqnreal, iMol), rMaxOld);
        if isnan(rMax)
            rMax = 1.0;
        end
        VMax = DiatPot(rMax, jqnreal, iMol);
    else
        rMax = 100.0;
        VMax = 0.0;
    end

end