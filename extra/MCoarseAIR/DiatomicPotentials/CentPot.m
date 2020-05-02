%% Centrifugal Contribution to Effective Diatomic Potential
%
function [Vc, dVc] = CentPot(r, jqn, iMol)

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

    global Syst
    
    iA    = Syst.Molecule(iMol).ToAtoms(1);
    jA    = Syst.Molecule(iMol).ToAtoms(2);
    
    mass  = [Syst.Atom(iA).Mass, Syst.Atom(jA).Mass];
    mu    = mass(1) * mass(2) / ( mass(1) + mass(2) );
    
    Vc_R2 = 1.d0 / mu * 0.5d0 * (jqn+0.5d0)^2;

    Vc    = Vc_R2 ./ r.^2;  
    dVc   = - 2.d0 .* Vc ./ r;  
  
end