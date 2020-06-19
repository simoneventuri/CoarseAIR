function [VProjf_CM] = Compute_VelocityFlux_CM(Vi, Vf, Masses, iProc)
    
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

    VTargi(:) = Masses(iTarg) * Vi(iTarg,:) ./ mTarg;
    VProji(:) = Masses(iProj) * Vi(iProj,:) ./ mProj;   
     
    
    if iProc == 1
        iTarg = [1,2];
        iProj = [3];
        Coeff = 1;
    elseif iProc == 2
        iProj = [1,3];
        iTarg = [2];
        Coeff = -1;
    elseif iProc == 3
        iProj = [2,3];
        iTarg = [1];
        Coeff = -1;
    end
    
    mTarg = sum(Masses(iTarg)); 
    mProj = sum(Masses(iProj));
    mu    = mTarg * mProj / (mTarg + mProj);

    VTargf(:) = Masses(iTarg) * Vf(iTarg,:) ./ mTarg;
    VProjf(:) = Masses(iProj) * Vf(iProj,:) ./ mProj;    
    
    VProjf_CM = norm(VProjf)*2.18e6;
    
%     VecVec = [norm(VProji), norm(VTargi), norm(VProjf), norm(VTargf)]*2.18e6
    
%     figure(1)
%     %quiver3(0.0,0.0,0.0,VTargi(1)*2.18e6,VTargi(2)*2.18e6,VTargi(3)*2.18e6,'k')
%     quiver3(0.0,0.0,0.0,VProji(1)*2.18e6,VProji(2)*2.18e6,VProji(3)*2.18e6,'r')
%     hold on
%     %quiver3(0.0,0.0,0.0,VTargf(1)*2.18e6,VTargf(2)*2.18e6,VTargf(3)*2.18e6,'b')
%     quiver3(0.0,0.0,0.0,VProjf(1)*2.18e6,VProjf(2)*2.18e6,VProjf(3)*2.18e6,'g')
%     hold off
%     xlim([-4000,4000])
%     ylim([-4000,4000])
%     zlim([-4000,4000])
%     pause
%   
end