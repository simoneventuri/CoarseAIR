% -- MATLAB --
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

function [cm, cmdot, xrel, vrel, dist, v, xkin] = ComputeTranslatedCoord(xx, xxdot)

  global Pair_to_Atoms AtomMass
  
  NAtoms  = 3;
  NSpace  = 3;
  NEqtVar = 6;

  
  cm    = zeros(3,2);
  cmdot = zeros(3,2);
  
  for i = 1:NSpace
    for j = 1:NAtoms
      if ( j == Pair_to_Atoms(1,1) ) || ( j == Pair_to_Atoms(1,2) )
        cm(i,1)    = cm(i,1)    + AtomMass(j) * xx(i,j);
        cmdot(i,1) = cmdot(i,1) + AtomMass(j) * xxdot(i,j);
      else
        cm(i,2)    = cm(i,2)    + AtomMass(j) * xx(i,j);
        cmdot(i,2) = cmdot(i,2) + AtomMass(j) * xxdot(i,j);
      end
    end
  end
  
  
  
  xt     =   0.d0;
  xp     =   0.d0;
  for i = 1:NAtoms
    if ( i == Pair_to_Atoms(1,1) ) || ( i == Pair_to_Atoms(1,2) )  
      xp = xp + AtomMass(i);
    else
      xt = xt + AtomMass(i);
    end
  end
  xpi = 1.d0 / xp;
  xti = 1.d0 / xt;

  for i = 1:NSpace
    cm(i,1)    = cm(i,1)    * xpi;
    cm(i,2)    = cm(i,2)    * xti;
    cmdot(i,1) = cmdot(i,1) * xpi;
    cmdot(i,2) = cmdot(i,2) * xti;
  end
  
  for i = 1:NSpace
    xrel(i)    = cm(i,1)    - cm(i,2);
    vrel(i)    = cmdot(i,1) - cmdot(i,2);
  end
  
  
  
  xnum = sum( xrel .* vrel );
  v    = sum( vrel.^2 );
  t    = - xnum ./ v;

  b   =   sqrt( sum( (xrel + t .* vrel ).^2 ) );

  xkin   = 0.d0;
  for i = 1:NSpace
    xkin = xkin + 5.d-1 * xp * cmdot(i,1)^2 + 5.d-1 * xt * cmdot(i,2)^2;
  end
  
  dist  =   sqrt( sum(xrel.^2) );
  v     =   sqrt( v );

  
end