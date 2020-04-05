function [V] = CO2_NASA_PES( R_O2, R_CO_1, R_CO_2 )

  roo    = R_O2
  fv1(1) = R_CO_1
  fv1(2) = R_CO_2
  
  

  fc  = xm(1) / (xm(1)+xm(2))
  fo  = xm(2) / (xm(1)+xm(2))
  xmv = xm(1) * xm(2) / (xm(1)+xm(2))
  xmr = xm(3) * (xm(1)+xm(2)) / (xm(1)+xm(2)+xm(3))
  xms = sqrt(xmr/xmv)
  
  
  
  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! This is Simply:    O2_DiaPot%Compute_Vd_dVd_O2( R_O2, V_O2=ham(3,3), dham(3,3,3) )
  vsr   = 64.d0 * exp(alphaco*roo) / roo
  dvsr  = vsr * (alphaco-(One/roo))
  vlr   = -c6o2 / (roo^6+damp6co)
  dvlr  = -Six * (roo^5) * vlr / (roo^6+damp6co)
  fact  = (roo^npowco) * exp(-ao2*roo)
  dfact = fact * ((dfloat(npowco)/roo)-ao2)
  cfcn  = coefo2(nfitco+1)
  dcfcn = Zero
  for k = nfitco:-1:2
        
    dcfcn = cfcn      + dcfcn * (roo-r0co)
    cfcn  = coefo2(k) + cfcn  * (roo-r0co)
    
  end   
  dcfcn   = dfact * cfcn + fact * dcfcn               
  cfcn    = coefo2(1)    + fact * cfcn
  
  V_O2        = cfcn + vsr + vlr %! + 187.888473490235 ????
  ham(3,3)    = V_O2
  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  voo    = 64.d0 * exp(alphax1*roo) / roo
  dvoo   = voo * (alphax1-(One/roo))
  voolr  = -ac6oo / (roo^6+damp6cox)
  dvoolr = -Six * (roo^5) * voolr / (roo^6+damp6cox)
 


  cc     = - (roo^2 - fv1(1)^2 - fv1(2)^2) * Half / (fv1(1)*fv1(2))
  dcc(1) = (One/fv1(2)) - (cc/fv1(1))
  dcc(2) = (One/fv1(1)) - (cc/fv1(2))
  dcc(3) = -roo / (fv1(1)*fv1(2))
      
      
      
  do j = 1:2
  
    vdiatm   = CO_DiaPot%DiatomicPotential( fv1(j) ) 
   
    vdiatm   = vdiatm + diatshift
    
    vco      = 48.d0*exp(alphax3*fv1(3-j))/fv1(3-j)
    dvco     = vco*(alphax3-(One/fv1(3-j)))
    
    vcolr    = -ac6b/(fv1(3-j)^6+damp6cox)
    dvcolr   = -Six*(fv1(3-j)^5)*vcolr/(fv1(3-j)^6+damp6cox)
    
    ham(3,3) = ham(3,3)+vco+vcolr
    
    vco      = 48.d0*exp(alphax2*fv1(3-j))/fv1(3-j)
    dvco     = vco*(alphax2-(One/fv1(3-j)))
    
    vcolr    = -ac6co/(fv1(3-j)^6+damp6cox)
    dvcolr   = -Six*(fv1(3-j)^5)*vcolr/(fv1(3-j)^6+damp6cox)
    
    vpair    = voo+vco
    
    rh       = fv1(j)*Half
    ro       = fv1(3-j)
    rbig2    = rh*rh+ro*ro-Two*rh*ro*cc
    rbig     = dsqrt(rbig2)
    rbigi    = One/rbig
    costheta = -(ro*ro-rh*rh-rbig2)*Half/(rbig*rh)
    drbigj   = (Half*rh-ro*(Half*cc+rh*dcc(j)))*rbigi
    drbig3j  = (ro-rh*(cc+ro*dcc(3-j)))*rbigi
    drbig3   = -rh*ro*dcc(3)*rbigi
    dcostj   = (Half/rbig)+(drbigj/rh)-costheta*(drbigj/rbig+Half/rh)
    dcost3j  = (((-ro/rbig)+drbig3j)/rh)-costheta*drbig3j/rbig
    dcost3   = drbig3/rh-costheta*drbig3/rbig
   
    fmat(1,1)  = One
    fmat(2,1)  = fv1(j)
    dfmat(1,1) = Zero
    dfmat(2,1) = One
    fmat(1,2)  = One
    fmat(2,2)  = costheta
    dfmat(1,2) = Zero
    dfmat(2,2) = One
       
    do k = 3:7
    
      km = k - 1
      fmat(k,1)  = fmat(2,1)  * fmat(km,1)
      dfmat(k,1) = dfmat(2,1) * fmat(km,1) + fmat(2,1) * dfmat(km,1)
      fmat(k,2)  = fmat(2,2)  * fmat(km,2)
      dfmat(k,2) = dfmat(2,2) * fmat(km,2) + fmat(2,2) * dfmat(km,2)
      
    end
       
    sum1   = Zero
    sum2   = Zero
    dsum11 = Zero
    dsum12 = Zero
    dsum21 = Zero
    dsum22 = Zero
       
    do k = 1:nparh
    
      kp     = k+nparh
      term   = fmat(ixpd(1,k),1)*fmat(ixpd(2,k),2)
      sum1   = sum1+par(k)*term
      dsum11 = dsum11+par(k)*dfmat(ixpd(1,k),1)*fmat(ixpd(2,k),2)
      dsum12 = dsum12+par(k)*fmat(ixpd(1,k),1)*dfmat(ixpd(2,k),2)
      sum2   = sum2+par(kp)*term
      dsum21 = dsum21+par(kp)*dfmat(ixpd(1,k),1)*fmat(ixpd(2,k),2)
      dsum22 = dsum22+par(kp)*fmat(ixpd(1,k),1)*dfmat(ixpd(2,k),2)
        
    end
       
    ex1      = exp(-gaus*((rbig-gaus1)^2))
    ex2      = exp(-gaus*((rbig-gaus2)^2))
    exr      = exp(-gdamp*fv1(j))
    calf     = exr*(ex1*sum1+ex2*sum2)
    dcalft   = Two*gaus*exr*((rbig-gaus1)*ex1*sum1+(rbig-gaus2)*ex2*sum2)
    dcalfj   = -drbigj*dcalft-gdamp*calf
    dcalf3j  = -drbig3j*dcalft
    dcalf3   = -drbig3*dcalft+exr*dcost3*(ex1*dsum12+ex2*dsum22)
    dcalfj   = dcalfj+exr*(ex1*(dsum11+dsum12*dcostj)+ex2*(dsum21+dsum22*dcostj))
    dcalf3j  = dcalf3j+exr*(ex1*dsum12*dcost3j+ex2*dsum22*dcost3j)
    exg      = exp(-gausv*calf)
    doff     = One/(One+exg)
    ddoff    = doff*doff*gausv*exg
    ddoffj   = ddoff*dcalfj
    ddoff3j  = ddoff*dcalf3j
    ddoff3   = ddoff*dcalf3
    realf    = One+calf*doff
    drealfj  = dcalfj*doff+calf*ddoffj
    drealf3j = dcalf3j*doff+calf*ddoff3j
    drealf3  = dcalf3*doff+calf*ddoff3
    ham(j,j) = vdiatm+vpair*realf+voolr+vcolr
       
  end
      
  cart(1,1) = Zero
  cart(2,1) = fv1(2)*fo
  cart(1,2) = Zero
  cart(2,2) = -fv1(2)*fc
  
  cthet = -(roo^2 - fv1(1)^2 - fv1(2)^2) * Half / (fv1(1)*fv1(2))
  
  dct1  = (One/fv1(2)) - (cthet/fv1(1))
  dct2  = (One/fv1(1)) - (cthet/fv1(2))
  dct3  = -roo / (fv1(1)*fv1(2))
  
  sthet = dsqrt(dabs(One-cthet^2))
  cott  = cthet/sthet
  
  cart(1,3) = fv1(1)*sthet
  cart(2,3) = fv1(1)*cthet+cart(2,2)
  
  d131 = dsqrt(dabs(One-cthet^2))-fv1(1)*dct1*cott
  d132 = -fv1(1)*dct2*cott
  d133 = -fv1(1)*dct3*cott
  d231 = cthet+fv1(1)*dct1
  d232 = fv1(1)*dct2-fc
  d233 = fv1(1)*dct3
  
  raw  = dsqrt(cart(1,3)^2+cart(2,3)^2)
  rawi = xms/raw
  rjac = xms*raw
  
  drj1 = (cart(1,3)*d131+cart(2,3)*d231)*rawi
  drj2 = (cart(1,3)*d132+cart(2,3)*d232)*rawi
  drj3 = (cart(1,3)*d133+cart(2,3)*d233)*rawi
  
  rho  = dsqrt(fv1(2)^2+rjac^2)
  rhoi = One/rho
  
  drh1 = drj1*rjac*rhoi
  drh2 = (fv1(2)+drj2*rjac)*rhoi
  drh3 = drj3*rjac*rhoi

  ex4 = exp(-xm(4)*rho) + 1.e-7
  ex5 = exp(-xm(5)*rho) + 1.e-7
  ex6 = exp(-xm(6)*rho) + 1.e-7
  ex7 = exp(-xm(7)*rho) + 1.e-7
  
  
  
  ham(1,2) = goff*ex4
  ham(2,1) = ham(1,2)
  ham(1,3) = goff3*ex6
  ham(2,3) = ham(1,3)
  ham(3,1) = ham(1,3)
  ham(3,2) = ham(1,3)
  
  xx = -xm(4) * (ham(1,2)-goff  * 1.e-7)
  xx = -xm(6) * (ham(1,3)-goff3 * 1.e-7)
  
  ndim = 3
  
  
      
  if (et0 ~= Zero)
      
    ndim  = 4
    ex1   = 48.d0 * exp(-alpc0*fv1(1)) / fv1(1)
    r15   = fv1(1)^5
    vlr1  = -One/(r15*fv1(1)+damp6co)
    part1 = ex1+c6c0*vlr1
    ex2   = 48.d0*exp(-alpc0*fv1(2))/fv1(2)
    r25   = fv1(2)^5
    vlr2  = -One/(fv1(2)*r25+damp6co)
    part2 = ex2+c6c0*vlr2
    ex3   = 64.d0*exp(-alpo0*roo)/roo
    r35   = roo^5
    vlr3  = -One/(roo*r35+damp6co)
    part3 = ex3+c6o0*vlr3
    v0o   = et0+part1+part2+part3
    xx    = Six*c6c0
    dv0o1 = -ex1*(alpc0+One/fv1(1))+xx*r15*vlr1*vlr1
    dv0o2 = -ex2*(alpc0+One/fv1(2))+xx*r25*vlr2*vlr2
    dv0o3 = -ex3*(alpo0+One/roo)+Six*c6o0*r35*vlr3*vlr3
       
    if (nfcnt > 0)
       
      scr(1)    = One
      scr(1+5)  = One
      scr(1+10) = One
      scr(2)    = fv1(1)-ret0
      scr(2+5)  = fv1(2)-ret0
      scr(2+10) = cc-cte0
      damp      = exp(-a1*(scr(2)^2+scr(7)^2)-a2*(scr(12)^2))
      damp2     = damp*Two
      xx        = -a2*scr(12)*damp2
      dm1       = -a1*scr(2)*damp2+xx*dcc(1)
      dm2       = -a1*scr(7)*damp2+xx*dcc(2)
      dm3       = xx*dcc(3)
        
      for k = 3:5
      
        k1  = k
        k1m = k1 - 1 
        scr(k1) = scr(k1m) * scr(2)
        k2  = k + 5
        k2m = k2 - 1
        scr(k2) = scr(k2m) * scr(2+5)
        k3  = k + 10
        k3m = k3 - 1 
        scr(k3) = scr(k3m) * scr(2+10)
         
      end
        
      voco=Zero
        
      for k = 1:nfcnt
        
        term0 = (scr(iqn(1,k)) * scr(iqn(2,k)+5) + scr(iqn(1,k)+5) * scr(iqn(2,k)))
        term  = scr(iqn(3,k)+10) * term0
        term3 = dfloat(iqn(3,k)-1) * scr(int(max(iqn(3,k)-One,One)+10.d0)) * term0
        term1 = scr(iqn(3,k)+10) * (dfloat(iqn(1,k)-1)*scr(int(max(iqn(1,k)-One,One))  ) * scr(iqn(2,k)+5) + dfloat(iqn(2,k)-1) * scr(iqn(1,k)+5) * scr(int(max(iqn(2,k)-One,One))  ))
        term2 = scr(iqn(3,k)+10) * (dfloat(iqn(1,k)-1)*scr(int(max(iqn(1,k)-One,One)+5)) * scr(iqn(2,k)  ) + dfloat(iqn(2,k)-1) * scr(iqn(1,k)  ) * scr(int(max(iqn(2,k)-One,One))+5))
       
        voco = voco + qff(k) * term
         
      end
      
      voco0 = voco
      voco  = v0o + damp * voco
        
    else
        
      voco  = v0o
      voco0 = Zero
         
    end
        
    ham(4,4) = voco
    ham(1,4) = vcupt0 * ex5
    ham(2,4) = vcupt0 * ex5
    ham(4,2) = vcupt0 * ex5
    ham(4,1) = vcupt0 * ex5
    xx = -xm(5) * (ham(1,4)-vcupt0 * 1.e-7)
  
    ham(3,4)    = vcupt1 * ex7
    
    xx = -xm(7) * (ham(3,4)-vcupt1 * 1.e-7)
    
    ham(4,3)    = ham(3,4)
        
  end
  
   
      
  for i = 1,ndim
    for j = 1,ndim
      ham2(j,i) = ham(j,i)
     end
  end
      
  call rs(4,ndim,ham,eig,1,ham,fv1,fv2,ierr)
      
  V = eig(1)
  
end
