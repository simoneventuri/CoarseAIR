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

%     NLevels1  = Syst.Molecule(iMol).NLevels;
%     LevelEeV1 = Syst.Molecule(iMol).LevelEeV;
%     Levelq1   = Syst.Molecule(iMol).T(Temp.iT).Levelq;
%     NLevels2  = Syst.Molecule(jMol).NLevels;
%     LevelEeV2 = Syst.Molecule(jMol).LevelEeV;
%     Levelq2   = Syst.Molecule(jMol).T(Temp.iT).Levelq;
% 
%     ExchChar = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/Exch_', num2str(iExch), '/');
%     %h5disp(Syst.HDF5_File, ExchChar)
%     RatesTemp = h5read(Syst.HDF5_File, ExchChar);
%     Rates.T(Temp.iT).ExchType(iExch).Exch = permute(RatesTemp, [2,1]);
%     fprintf(['  Rates.T(' num2str(Temp.iT) ').ExchType(' num2str(iExch) ').Exch, size: (' num2str(size(Rates.T(Temp.iT).ExchType(iExch).Exch)) ') \n'])
% 
%     fprintf(['  Computing Electronic and Translational Partition Function\n'])
%     RxLxIdx = Syst.RxLxIdx(2+iExch,:);
%     Qt      = 1.0;
%     Qe      = 1.0;
%     for iComp = 1:Syst.NComp
%         Syst.CFDComp(iComp).Qt = Param.Plnck / sqrt( (2.0*pi) * (Syst.CFDComp(iComp).Mass*Param.AMUToKg) * Param.KJK * Temp.TNow );
%         Qt = Qt * (Syst.CFDComp(iComp).Qt)^RxLxIdx(iComp);
%         Qe = Qe * (Syst.CFDComp(iComp).Qe)^RxLxIdx(iComp);
%     end
% 
%     ExchChar = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/Exch_', num2str(Syst.ToOtherExch(iExch)), '/');
%     fprintf(['  Finishing Reading Exchange Rates from ', OtherSyst.HDF5_File, ' (Exch. Nb. ', num2str(Syst.ToOtherExch(iExch)), ') \n'])
%     %h5disp(Syst.HDF5_File_OtherExch, ExchChar)
%     RatesTemp  = h5read(OtherSyst.HDF5_File, ExchChar);   %%% Note: We are not permuting it 
%     for jLevel = 1:NLevels2
%         for iLevel = 1:NLevels1
%             if (LevelEeV1(iLevel) < LevelEeV2(jLevel))
%                 Rates.T(Temp.iT).ExchType(iExch).Exch(iLevel,jLevel) = RatesTemp(iLevel,jLevel) / (Levelq1(iLevel) / Levelq2(jLevel) * Qe * Qt);
%             end
%         end
%     end
% 
%     fprintf(['  Saving Merged Rates in HDF5 File \n'])
%     h5create(Syst.HDF5_File, ExchCharMerged, [NLevels2 NLevels1])
%     h5write(Syst.HDF5_File,  ExchCharMerged, Rates.T(Temp.iT).ExchType(iExch).Exch')
    
    
    fprintf('====================================================\n\n')
    
end