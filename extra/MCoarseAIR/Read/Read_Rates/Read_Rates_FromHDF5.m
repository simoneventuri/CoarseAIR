%% The Function reads the Rates from the HD5 File
%
%  Input Global Var: - Temp.TNowChar
%                    - Syst.HDF5_File
%
function Read_Rates_FromHDF5()    

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
    
    global Rates Syst Temp Input Param
    
    fprintf('  = Read_Rates_FromHDF5 ================== T = %i K\n', Temp.TNow)
    fprintf('  ====================================================\n')
    fprintf('  Reading Rates in HDF5 Format \n' )
    fprintf(['  Reading from File: ' Syst.HDF5_File '\n'] )

    
    if (Syst.NAtoms == 3)

        if (Input.Kin.ReadRatesProc(1))
            DissChar       = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/Diss/');
            %h5disp(Syst.HDF5_File, DissChar)
            RatesTemp                 = h5read(Syst.HDF5_File, DissChar);
            Rates.T(Temp.iT).Diss     = permute(RatesTemp, [2,1]);
            fprintf(['  Rates.T(' num2str(Temp.iT) ').Diss, size: (' num2str(size(Rates.T(Temp.iT).Diss)) ') \n'])
        end
        
        if (Input.Kin.ReadRatesProc(2))
            InelChar       = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/Inel/');
            %h5disp(Syst.HDF5_File, InelChar)
            RatesTemp                 = h5read(Syst.HDF5_File, InelChar);
            Rates.T(Temp.iT).Inel     = permute(RatesTemp, [2,1]);
            fprintf(['  Rates.T(' num2str(Temp.iT) ').Inel, size: (' num2str(size(Rates.T(Temp.iT).Inel)) ') \n'])
        end
        
        for iExch = 1:Syst.NProc-2
            if (Input.Kin.ReadRatesProc(2+iExch))
                
                if (Syst.ToOtherExch(iExch) > 0)
                    ExchCharMerged = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/Exch_', num2str(iExch), '_Merged/');

                    if (Input.Kin.OtherExchInHDF5)
                        %h5disp(Syst.HDF5_File, ExchCharMerged)
                        RatesTemp = h5read(Syst.HDF5_File, ExchCharMerged);
                        Rates.T(Temp.iT).ExchType(iExch).Exch = permute(RatesTemp, [2,1]);
                        fprintf(['  Rates.T(' num2str(Temp.iT) ').ExchType(' num2str(iExch) ').Exch, size: (' num2str(size(Rates.T(Temp.iT).ExchType(iExch).Exch)) ') \n'])

                    else
                        iMol      = 1;
                        NLevels1  = Syst.Molecule(iMol).NLevels;
                        LevelEeV1 = Syst.Molecule(iMol).LevelEeV;
                        Levelq1   = Syst.Molecule(iMol).T(Temp.iT).Levelq;
                        jMol      = Syst.ExchToMol(iExch);
                        NLevels2  = Syst.Molecule(jMol).NLevels;
                        LevelEeV2 = Syst.Molecule(jMol).LevelEeV;
                        Levelq2   = Syst.Molecule(jMol).T(Temp.iT).Levelq;

                        ExchChar = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/Exch_', num2str(iExch), '/');
                        %h5disp(Syst.HDF5_File, ExchChar)
                        RatesTemp = h5read(Syst.HDF5_File, ExchChar);
                        Rates.T(Temp.iT).ExchType(iExch).Exch = permute(RatesTemp, [2,1]);
                        fprintf(['  Rates.T(' num2str(Temp.iT) ').ExchType(' num2str(iExch) ').Exch, size: (' num2str(size(Rates.T(Temp.iT).ExchType(iExch).Exch)) ') \n'])

                        fprintf(['  Computing Electronic and Translational Partition Function\n'])
                        RxLxIdx = Syst.RxLxIdx(2+iExch,:);
                        Qt      = 1.0;
                        Qe      = 1.0;
                        for iComp = 1:Syst.NComp
                            Syst.CFDComp(iComp).Qt = Param.Plnck / sqrt( (2.0*pi) * (Syst.CFDComp(iComp).Mass*Param.AMUToKg) * Param.KJK * Temp.TNow );
                            Qt = Qt * (Syst.CFDComp(iComp).Qt)^RxLxIdx(iComp);
                            Qe = Qe * (Syst.CFDComp(iComp).Qe)^RxLxIdx(iComp);
                        end

                        ExchChar = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/Exch_', num2str(Syst.ToOtherExch(iExch)), '/');
                        fprintf(['  Finishing Reading Exchange Rates from ', Syst.HDF5_File_OtherExch, ' (Exch. Nb. ', num2str(Syst.ToOtherExch(iExch)), ') \n'])
                        %h5disp(Syst.HDF5_File_OtherExch, ExchChar)
                        RatesTemp  = h5read(Syst.HDF5_File_OtherExch, ExchChar);   %%% Note: We are not permuting it 
                        for jLevel = 1:NLevels2
                            for iLevel = 1:NLevels1
                                if (LevelEeV1(iLevel) < LevelEeV2(jLevel))
                                    Rates.T(Temp.iT).ExchType(iExch).Exch(iLevel,jLevel) = RatesTemp(iLevel,jLevel) / (Levelq1(iLevel) / Levelq2(jLevel) * Qe * Qt);
                                end
                            end
                        end

                        fprintf(['  Saving Merged Rates in HDF5 File \n'])
                        h5create(Syst.HDF5_File, ExchCharMerged, [NLevels2 NLevels1])
                        h5write(Syst.HDF5_File,  ExchCharMerged, Rates.T(Temp.iT).ExchType(iExch).Exch')

                    end

                end
            
            else
                
                Rates.T(Temp.iT).ExchType(iExch).Exch = 0.0;
                
            end
            
        end
        
    else
    
        DissChar       = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/Diss/');
        %h5disp(Syst.HDF5_File, DissChar)
        RatesTemp                 = h5read(Syst.HDF5_File, DissChar);
        Rates.T(Temp.iT).Diss     = permute(RatesTemp, [3,2,1]);
        fprintf(['  Rates.T(' num2str(Temp.iT) ').Diss, size: (' num2str(size(Rates.T(Temp.iT).Diss)) ') \n'])

        DissCharInel   = strcat('/T_', Temp.TNowChar, '_', Temp.TNowChar, '/Rates/DissInel/');
        %h5disp(Syst.HDF5_File, DissCharInel)
        RatesTemp                 = h5read(Syst.HDF5_File, DissCharInel);
        Rates.T(Temp.iT).DissInel = permute(RatesTemp, [4,3,2,1]);
        fprintf(['  Rates.T(' num2str(Temp.iT) ').DissInel, size: (' num2str(size(Rates.T(Temp.iT).DissInel)) ') \n'])

    end
    
    
    fprintf('  ====================================================\n\n')
    
end