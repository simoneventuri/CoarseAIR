function [Theta_CM, TOF] = Compute_ScatteringAngle_CM(Vi, Vf, Masses, iProc, Dist)
    

    iTarg = [1,2];
    iProj = [3];
    
    mTarg = sum(Masses(iTarg)); 
    mProj = sum(Masses(iProj));
    mu    = mTarg * mProj / (mTarg + mProj);

    VTargi(:) = Masses(iTarg) * Vi(iTarg,:) ./ mTarg;
    VProji(:) = Masses(iProj) * Vi(iProj,:) ./ mProj;   
    
    vReli_CM  = VProji - VTargi;
    
    Theta_i  = acos( dot(VTargi, VProji) / (vecnorm(VTargi) * vecnorm(VProji)) ) / pi * 180.0;
    
    
    if iProc == 1
        iTarg = [1,2];
        iProj = [3];
        Coeff = 1;
    elseif iProc == 2
        iProj = [2];
        iTarg = [1,3];
        Coeff = -1;
    elseif iProc == 3
        iProj = [1];
        iTarg = [2,3];
        Coeff = -1;
    end
    
    mTarg = sum(Masses(iTarg)); 
    mProj = sum(Masses(iProj));
    mu    = mTarg * mProj / (mTarg + mProj);

    VTargf(:) = Masses(iTarg) * Vf(iTarg,:) ./ mTarg;
    VProjf(:) = Masses(iProj) * Vf(iProj,:) ./ mProj;   
    
    vRelf_CM  = Coeff * (VProjf - VTargf);
    
    Theta_f  = acos( dot(VTargf, VProjf) / (vecnorm(VTargf) * vecnorm(VProjf)) ) / pi * 180.0;
    
    
    
    Theta_CM  = acos( dot(vReli_CM, vRelf_CM) / (vecnorm(vReli_CM) * vecnorm(vRelf_CM)) ) / pi * 180.0;
   
    
    
    TOF       = Dist / norm(VProjf);
    
    
    
end