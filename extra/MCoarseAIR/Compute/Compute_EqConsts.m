%% The Function Computes the Equilibrium Constants
%
function Compute_EqConsts()  
    
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
    
    global Syst Param Temp Input
    
    
    fprintf('= Compute_EqConsts ====================== T = %i K\n', Temp.TNow)
    fprintf('====================================================\n')
    
    
    RxLxIdx = Syst.RxLxIdx(1,:);

    Syst.T(Temp.iT).KEqDiss = 1.0;
    for iComp = 1:Syst.NComp
        Syst.CFDComp(iComp).T(Temp.iT).Qt = Param.Plnck / sqrt( (2.0*pi) * (Syst.CFDComp(iComp).Mass*Param.AMUToKg) * Param.KJK * Temp.TNow );

        if Syst.CFDComp(iComp).ToMol > 0
            Syst.CFDComp(iComp).T(Temp.iT).QRot = Syst.Molecule(Syst.CFDComp(iComp).ToMol).T(Temp.iT).QRot;
        else 
            Syst.CFDComp(iComp).T(Temp.iT).QRot = 1.0;
        end

        Syst.T(Temp.iT).KEqDiss = Syst.T(Temp.iT).KEqDiss * (Syst.CFDComp(iComp).T(Temp.iT).Qt * Syst.CFDComp(iComp).T(Temp.iT).QRot * Syst.CFDComp(iComp).Qe)^RxLxIdx(iComp);
    end
    fprintf('KEq for Dissociation/Recombination = %e \n', Syst.T(Temp.iT).KEqDiss );
    
    
%     RxLxIdx = Syst.RxLxIdx(2,:);
% 
%     Syst.T(Temp.iT).KEqInel = 1.0;
%     for iComp = 1:Syst.NComp
%         Syst.CFDComp(iComp).T(Temp.iT).Qt = Param.Plnck / sqrt( (2.0*pi) * (Syst.CFDComp(iComp).Mass*Param.AMUToKg) * Param.KJK * Temp.TNow );
% 
%         if Syst.CFDComp(iComp).ToMol > 0
%             Syst.CFDComp(iComp).T(Temp.iT).QRot = Syst.Molecule(Syst.CFDComp(iComp).ToMol).T(Temp.iT).QRot;
%         else 
%             Syst.CFDComp(iComp).T(Temp.iT).QRot = 1.0;
%         end
% 
%         Syst.T(Temp.iT).KEqInel = Syst.T(Temp.iT).KEqInel * (Syst.CFDComp(iComp).T(Temp.iT).Qt * Syst.CFDComp(iComp).T(Temp.iT).QRot * Syst.CFDComp(iComp).Qe)^RxLxIdx(iComp);
%     end
    
    iExch = 1;
    for iProc = 3:Syst.NProc 
        RxLxIdx = Syst.RxLxIdx(iProc,:);

        Syst.T(Temp.iT).KEqExch(iExch) = 1.0;
        for iComp = 1:Syst.NComp
            Syst.CFDComp(iComp).T(Temp.iT).Qt = Param.Plnck / sqrt( (2.0*pi) * (Syst.CFDComp(iComp).Mass*Param.AMUToKg) * Param.KJK * Temp.TNow );

            if Syst.CFDComp(iComp).ToMol > 0
                Syst.CFDComp(iComp).T(Temp.iT).QRot = Syst.Molecule(Syst.CFDComp(iComp).ToMol).T(Temp.iT).QRot;
            else 
                Syst.CFDComp(iComp).T(Temp.iT).QRot = 1.0;
            end

            Syst.T(Temp.iT).KEqExch(iExch) = Syst.T(Temp.iT).KEqExch(iExch) * (Syst.CFDComp(iComp).T(Temp.iT).Qt * Syst.CFDComp(iComp).T(Temp.iT).QRot * Syst.CFDComp(iComp).Qe)^RxLxIdx(iComp);
        end
        fprintf('KEq for Echange Nb. %i = %e \n', iExch, Syst.T(Temp.iT).KEqExch(iExch) );
        
        iExch = iExch + 1;
    end

    if (Syst.NProc > 2)
        [status,msg,msgID] = mkdir(Input.Paths.SaveDataFldr);
        FileName          = strcat(Input.Paths.SaveDataFldr, '/EqConstants.csv');
        if exist(FileName, 'file')
            fileID1  = fopen(FileName,'a');
        else
            fileID1  = fopen(FileName,'w');
            if (Syst.NProc == 3)
                HeaderStr = strcat('# T [K], KDiss, KExch1\n');
            elseif (Syst.NProc == 4)
                HeaderStr = strcat('# T [K], KDiss, KExch1, KExch2\n');
            end
            fprintf(fileID1,HeaderStr);
        end
        if (Syst.NProc == 3)
            fprintf(fileID1,'%e,%e,%e\n',    Temp.TNow, Syst.T(Temp.iT).KEqDiss, Syst.T(Temp.iT).KEqExch(1) );
        elseif (Syst.NProc == 4)
            fprintf(fileID1,'%e,%e,%e,%e\n', Temp.TNow, Syst.T(Temp.iT).KEqDiss, Syst.T(Temp.iT).KEqExch(1), Syst.T(Temp.iT).KEqExch(2)  );
        end
        fclose(fileID1);
    end
    
    

    fprintf('====================================================\n\n')        
    
end