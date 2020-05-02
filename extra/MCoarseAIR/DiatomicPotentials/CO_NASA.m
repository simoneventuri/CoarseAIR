%% CO Diatomic Potential from NASA Ames (Scwhenke et al., 2016)
%
function [V, dV] = CO_NASA(r)
   
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
    
    Zero = 0.d0;
    One  = 1.d0;
    Six  = 6.d0;

    nf     = 12;
    coef   = [ -4.5761074811943886e-2, 0.22111734303145603d0, -0.34188485707539029d0, 0.21975372852913833d0, -0.18415829369278044d0, 0.19162121009025507d0, ...
             -0.14681584359973315d0, 7.3117032359438006e-2, -2.2588593079026137e-2, 4.1361679444378470e-3, -4.0747646620877898e-4, 1.6857235957160763e-5];
    c6     =  284.57619366543872d0;
    cinf   = -112.90087267991528d0;
    alpha  = -2.2343545089620513d0;
    damp6  =  418.89225320735858d0;
    npow   =  9;
    a      =  4.18604651162790730d0;
    frel   =  0.12114229643464522d0;
    R0     =  2.14999999999999999d0;

    nfm    = nf - 1;
    freli  = One ./ frel;

    vsr  = 48.0d0 .* exp(alpha .* r) ./ r;
    dvsr = (alpha - (One ./ r)) .* vsr;
    vlr  = -c6 ./ (r.^6 + damp6);
    dvlr = Six .*  r.^5 .* vlr .* vlr ./ c6;

    vint  = coef(nf);
    dvint = Zero;

    for j = nfm:-1:1
        dvint = vint    + dvint .* (r-R0);
        vint  = coef(j) + vint  .* (r-R0);
    end

    dvint = dvint .* r.^npow .* exp(-a .* r) .* freli;
    vint  = vint  .* r.^npow .* exp(-a .* r) .* freli;
    dvint = dvint + ((npow ./ r) - a) .* vint;

    V  = vsr  + vlr  + vint;
    dV = dvsr + dvlr + dvint;
    
    V  = V  .* Param.EhToeV;
    dV = dV .* Param.EhToeV;
end