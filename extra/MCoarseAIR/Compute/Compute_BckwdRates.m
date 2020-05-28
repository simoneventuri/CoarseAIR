%% The Function Computes the Backward Rates
%
%  Input Global Var: - Temp.TNowChar
%                    - Syst.HDF5_File
%
function [FinalRateMat] = Compute_BckwdRates(RateMat, iProc, LevelEeV1, Levelq1, LevelEeV2, Levelq2)    

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
        
    global Param Temp Syst
    
    fprintf('= Compute_BckwdRates =================== T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    
    
    RxLxIdx = Syst.RxLxIdx(iProc,:);
    
    Qt = 1.0;
    Qe = 1.0;
    for iComp = 1:Syst.NComp
        Syst.CFDComp(iComp).Qt = Param.Plnck / sqrt( (2.0*pi) * (Syst.CFDComp(iComp).Mass*Param.AMUToKg) * Param.KJK * Temp.TNow );
        
        Qt = Qt * (Syst.CFDComp(iComp).Qt)^RxLxIdx(iComp);
        Qe = Qe * (Syst.CFDComp(iComp).Qe)^RxLxIdx(iComp);
    end

    
    if (iProc == 1)    
        
        FinalRateMat = RateMat .* Levelq1 .* Qe .* Qt;
        
    else

        NLevels1     = size(RateMat,1);
        NLevels2     = size(RateMat,2);
        FinalRateMat = RateMat .* 0.0;

        for iLevel=1:NLevels1
            for jLevel=1:NLevels2
                if LevelEeV1(iLevel) >= LevelEeV2(jLevel) 
                    FinalRateMat(iLevel,jLevel) = RateMat(iLevel,jLevel);
                else
                    FinalRateMat(iLevel,jLevel) = RateMat(jLevel,iLevel) * Levelq1(iLevel) / Levelq2(jLevel) * Qe * Qt;
                end
            end
        end
    
    end
    
    fprintf('====================================================\n\n')
    
end