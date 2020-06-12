function [EColl_CM] = Compute_CollisionEnergy_CM(V, Masses, iProc)
    
    EhToKcalMol   = 627.5096080305927;

    if iProc == 1
        iTarg = [1,2];
        iProj = [3];
    elseif iProc == 2
        iTarg = [2];
        iProj = [1,3];
    elseif iProc == 3
        iTarg = [1];
        iProj = [2,3];
    end
    
    mTarg = sum(Masses(iTarg)); 
    mProj = sum(Masses(iProj));
    mu    = mTarg * mProj / (mTarg + mProj);

    VTarg(:) = Masses(iTarg) * V(iTarg,:) ./ mTarg;
    VProj(:) = Masses(iProj) * V(iProj,:) ./ mProj;   
    
    vRel_CM    = norm(VTarg-VProj);
    EColl_CM   = 0.5 * mu * vRel_CM^2 * EhToKcalMol;

end