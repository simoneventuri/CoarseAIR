%% The Function Groups the Rate Coefficients
%        
function Compute_GroupedRates()
      
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

    
    global Syst Rates Param Temp Input

    fprintf('= Compute_GroupedRates =============================\n')
    fprintf('====================================================\n')



    iMol        = 1;
    NLevels1    = Syst.Molecule(iMol).NLevels;
    LevelToBin1 = Syst.Molecule(iMol).LevelToGroupOut;
    NBins1      = Syst.Molecule(iMol).NGroupsOut;
    ExpVec1     = Syst.Molecule(iMol).T(Temp.iT).Levelq ./ sum(Syst.Molecule(iMol).T(Temp.iT).Levelq);
    ExpMat      = kron(ExpVec1,1.d0./ExpVec1');  



    QBin = zeros(NBins1,1);
    for iLevel=1:NLevels1
       iBin       = LevelToBin1(iLevel); 
       QBin(iBin) = QBin(iBin) + Syst.Molecule(iMol).T(Temp.iT).Levelq(iLevel);
    end


    if (Input.Kin.ReadRatesProc(1))
        KDiss      = Rates.T(Temp.iT).Diss(:,1);
        KDiss_Bins = zeros(NBins1,1);
        for iLevel=1:NLevels1
           iBin             = LevelToBin1(iLevel);
           KDiss_Bins(iBin) = KDiss_Bins(iBin) + KDiss(iLevel) * Syst.Molecule(iMol).T(Temp.iT).Levelq(iLevel) / QBin(iBin); 
        end
        Rates.T(Temp.iT).DissOut = KDiss_Bins;
        fprintf('Grouped Dissociation Rates\n' )
        
        Rates.T(Temp.iT).DissOutRecon = zeros(NLevels1,1);
        for iLevel=1:NLevels1
           iBin             = LevelToBin1(iLevel);
           Rates.T(Temp.iT).DissOutRecon(iLevel) = Rates.T(Temp.iT).DissOut(iBin); 
        end
        fprintf('Reconstructed RoVibrational Specific Rates\n' )
    end

    
    if (Input.Kin.ReadRatesProc(2))
        KInel      = Rates.T(1).Inel;
        KInelExo   = tril(KInel) + tril(KInel .* ExpMat,-1)';
        KInel_Bins = zeros(NBins1,NBins1);
        for iLevel=1:NLevels1
            iBin = LevelToBin1(iLevel);
            for jLevel=1:NLevels1
                jBin = LevelToBin1(jLevel);
                KInel_Bins(iBin,jBin) = KInel_Bins(iBin,jBin) + KInelExo(iLevel, jLevel) * Syst.Molecule(iMol).T(Temp.iT).Levelq(iLevel) / QBin(iBin);
            end
        end
        Rates.T(Temp.iT).InelOut = KInel_Bins;
        clear ExpMat
        fprintf('Grouped Inelastic Rates\n' )
    end

%     for iExch = 1:size(Syst.ExchToMol,1)
%         jMol        = Syst.ExchToMol(iExch);
%         NLevels2    = Syst.Molecule(jMol).NLevels;
%         LevelToBin2 = Syst.Molecule(jMol).LevelToGroupOut;
%         NBins2      = Syst.Molecule(jMol).NGroupsOut;
%         ExpVec2     = Syst.Molecule(jMol).T(Temp.iT).Levelq ./ sum(Syst.Molecule(jMol).T(Temp.iT).Levelq);
%         
%         if (Input.Kin.ReadRatesProc(2+iExch))
%             KExch       = Rates.T(Temp.iT).ExchType(iExch).Exch;
%             ExpMat      = kron(ExpVec1, 1.d0./ExpVec2');  
%             KExchExo    = tril(KExch) + tril(KExch .* ExpMat,-1)';
%             KExch_Bins  = zeros(NBins1,NBins2);
%             for iLevel=1:NLevels1
%                 iBin = LevelToBin1(iLevel);
%                 for jLevel=1:NLevels2
%                     jBin = LevelToBin2(jLevel);
%                     KExch_Bins(iBin,jBin) = KExch_Bins(iBin,jBin) + KExchExo(iLevel, jLevel) * Syst.Molecule(iMol).T(Temp.iT).Levelq(iLevel) / QBin(iBin);
%                 end
%             end
%             Rates.T(Temp.iT).ExchType(iExch).ExchOut = KExch_Bins;
%             fprintf('Grouped Rates for Exchange Nb %i\n', iExch )
%         end
%         
%     end

    
    fprintf('====================================================\n\n')
    
end