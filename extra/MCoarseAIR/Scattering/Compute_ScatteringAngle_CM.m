function [Theta_CM] = Compute_ScatteringAngle_CM(Vi, Vf, Masses, iProc)
    
    EhToKcalMol   = 627.5096080305927;

   
    
    iTarg = [1,2];
    iProj = [3];
    
    mTarg = sum(Masses(iTarg)); 
    mProj = sum(Masses(iProj));
    mu    = mTarg * mProj / (mTarg + mProj);

    VTargi(:) = Masses(iTarg) * Vi(iTarg,:) ./ mTarg;
    VProji(:) = Masses(iProj) * Vi(iProj,:) ./ mProj;   
    
    vReli_CM  = VTargi-VProji;
    
    
    
%     if iProc == 1
        iTarg = [1,2];
        iProj = [3];
        Coeff = 1.0;
%     elseif iProc == 2
%         iTarg = [1,3];
%         iProj = [2];
%         Coeff = -1.0;
%     elseif iProc == 3
%         iTarg = [2,3];
%         iProj = [1];
%         Coeff = -1.0;
%     end
    
    mTarg = sum(Masses(iTarg)); 
    mProj = sum(Masses(iProj));
    mu    = mTarg * mProj / (mTarg + mProj);

    VTargf(:) = Masses(iTarg) * Vf(iTarg,:) ./ mTarg;
    VProjf(:) = Masses(iProj) * Vf(iProj,:) ./ mProj;   
    
    vRelf_CM  = Coeff*(VTargf-VProjf);
    
     
    
    Theta_CM  = acos( sum(vReli_CM .* vRelf_CM) / (norm(vReli_CM) * norm(vRelf_CM)) ) / pi * 180.0;

end