%% The Function Computes the Equilibrium Constants
%
function Compute_EqConsts()  
    
    global Syst Param Temp
    
    RxLxIdx = Syst.RxLxIdx(Syst.NProc,:);
    
    Syst.KEqExch = 1.0;
    for iComp = 1:Syst.NComp
        Syst.CFDComp(iComp).Qt = Param.Plnck / sqrt( (2.0*pi) * (Syst.CFDComp(iComp).Mass*Param.AMUToKg) * Param.KJK * Temp.TNow );
        
        if Syst.CFDComp(iComp).ToMol > 0
            Syst.CFDComp(iComp).QRot = Syst.Molecule(Syst.CFDComp(iComp).ToMol).T(Temp.iT).QRot;
        else 
            Syst.CFDComp(iComp).QRot = 1.0;
        end
        
        Syst.KEqExch = Syst.KEqExch * (Syst.CFDComp(iComp).Qt * Syst.CFDComp(iComp).QRot * Syst.CFDComp(iComp).Qe)^RxLxIdx(iComp);
    end

end