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

function [xx, xxdot, Rbs] = ComputeCoordAndVeloc(PaQ)

  global AtomMass
    
  mMiMn(1:2) =  - AtomMass(1:2) ./ AtomMass(3);
  mMiMn(3)   = 0.d0;

  N       = 3;
  NSpace  = 3;
  NEqtVar = 6;

  
  for iS = 1:NSpace
    for iA = 1:N-1
      xxdot(iS,iA)  =   PaQ(iS+(iA-1)*3);
      xx(iS,iA)     =   PaQ(iS+(iA-1)*3 + NEqtVar);
    end
  end

  for iS = 1:NSpace
    xxdot(iS,N)     =   0.d0;
    xx(iS,N)        =   0.d0;
    for iA = 1:N-1
      xxdot(iS,N)   =   xxdot(iS,N) + mMiMn(iA) * PaQ(iS + (iA-1)*3);
      xx(iS,N)      =   xx(iS,N)    + mMiMn(iA) * PaQ(iS + (iA-1)*3 + NEqtVar);
    end
  end
  
  iP = 0;
  for i = 1:N-1
    for j = i+1:N
      iP = iP + 1;
      Rbs(iP) = ( xx(1,i) - xx(1,j) )^2 + ( xx(2,i) - xx(2,j) )^2 + ( xx(3,i) - xx(3,j) )^2;
    end
  end
  
end