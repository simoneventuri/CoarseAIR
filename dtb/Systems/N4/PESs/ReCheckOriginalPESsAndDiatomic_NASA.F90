cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine n4fitd(xyz1,xyz2,xyz3,xyz4,fit,dfit)
      implicit real*8 (a-h,o-z)
c
c     N4 PES including dissociation and N4 Td
c     DWS 8/4/09
c     xyzi: cartesian coordinates of atom i in a.u.
c
      parameter (lwork=50)
      dimension xyz1(3),xyz2(3),xyz3(3),xyz4(4)
      dimension vmat(4,4),eig(4),fv1(4),fv2(4),cpy(4,4),
     $     dfit(12),dpart(12,4),vec(4),work(lwork),
     $     iwork(20),ifail(4),dfit2(12),cart(3,4),vmatc(4,4)            7d22s14

      integer(kind=4)                  :: ic

c$$$      common/dwsdm/idoit
      save                                                              9d3s09
      do i=1,3                                                          9d8s09
         cart(i,1)=xyz1(i)                                                9d8s09
         cart(i,2)=xyz2(i)                                                9d8s09
         cart(i,3)=xyz3(i)                                                9d8s09
         cart(i,4)=xyz4(i)                                                9d8s09
      end do                                                            9d8s09
c$$$c$$$      do i=1,12                                                         9d8s09
c$$$c$$$       dfit2(i)=0d0                                                     9d8s09
c$$$c$$$      end do                                                            9d8s09
c$$$c$$$      wgt=-1d0/3d0                                                      9d8s09
c$$$c$$$      eps=1d-3                                                          9d8s09
c$$$c$$$      do ipass=1,2                                                      9d8s09
c$$$c$$$       do i=1,12                                                        9d8s09
c$$$c$$$        csav=cart(i,1)                                                  9d8s09
c$$$c$$$        cart(i,1)=csav+eps                                              9d8s09
c$$$c$$$        do ii=1,4
c$$$c$$$         do j=1,4
c$$$c$$$          vmat(j,ii)=1d-5
c$$$c$$$         end do
c$$$c$$$        end do
c$$$c$$$        call vtdd(cart,cart(1,2),cart(1,3),cart(1,4),vmat(1,1),dpart)   9d10s09
c$$$c$$$        call vn2n2d(cart,cart(1,2),cart(1,3),cart(1,4),vmat(2,2),dpart,       9d8s09
c$$$c$$$     $       shiftof0,0)                                                9d8s09
c$$$c$$$        call vn2n2d(cart,cart(1,3),cart(1,2),cart(1,4),vmat(3,3),dpart,       9d8s09
c$$$c$$$     $       shiftof0,0)                                                9d8s09
c$$$c$$$        call vn2n2d(cart,cart(1,4),cart(1,2),cart(1,3),vmat(4,4),dpart,       9d8s09
c$$$c$$$     $       shiftof0,0)                                                9d8s09
c$$$c$$$        vmat(1,1)=vmat(1,1)-shiftof0                                      9d3s09
c$$$c$$$        do ii=2,4
c$$$c$$$         vmat(1,ii)=6.0847866d-02
c$$$c$$$         vmat(ii,1)=6.0847866d-02
c$$$c$$$        end do
c$$$c$$$        zero=0d0
c$$$c$$$        call dsyevx('N','I','L',4,vmat,4,vl,ul,1,1,zero,neig,
c$$$c$$$     $     eig,vec,4,work,lwork,iwork,ifail,info)
c$$$c$$$        vpotp=eig(1)
c$$$c$$$        cart(i,1)=csav-eps                                              9d8s09
c$$$c$$$        do ii=1,4
c$$$c$$$         do j=1,4
c$$$c$$$          vmat(j,ii)=1d-5
c$$$c$$$         end do
c$$$c$$$        end do
c$$$c$$$        call vtdd(cart,cart(1,2),cart(1,3),cart(1,4),vmat(1,1),dpart)   9d10s09
c$$$c$$$        call vn2n2d(cart,cart(1,2),cart(1,3),cart(1,4),vmat(2,2),dpart,       9d8s09
c$$$c$$$     $       shiftof0,0)                                                9d8s09
c$$$c$$$        call vn2n2d(cart,cart(1,3),cart(1,2),cart(1,4),vmat(3,3),dpart,       9d8s09
c$$$c$$$     $       shiftof0,0)                                                9d8s09
c$$$c$$$        call vn2n2d(cart,cart(1,4),cart(1,2),cart(1,3),vmat(4,4),dpart,       9d8s09
c$$$c$$$     $       shiftof0,0)                                                9d8s09
c$$$c$$$        vmat(1,1)=vmat(1,1)-shiftof0                                      9d3s09
c$$$c$$$        do ii=2,4
c$$$c$$$         vmat(1,ii)=6.0847866d-02
c$$$c$$$         vmat(ii,1)=6.0847866d-02
c$$$c$$$        end do
c$$$c$$$        zero=0d0
c$$$c$$$        call dsyevx('N','I','L',4,vmat,4,vl,ul,1,1,zero,neig,
c$$$c$$$     $     eig,vec,4,work,lwork,iwork,ifail,info)
c$$$c$$$        vpotm=eig(1)
c$$$c$$$        cart(i,1)=csav                                                  9d8s09
c$$$c$$$        dfit2(i)=dfit2(i)+wgt*(vpotp-vpotm)*0.5d0/eps                   9d8s09
c$$$c$$$       end do                                                           9d8s09
c$$$c$$$       wgt=4d0/3d0                                                      9d8s09
c$$$c$$$       eps=eps*0.5d0                                                    9d8s09
c$$$c$$$      end do                                                            9d8s09
c$$$      write(6,*)('carts: ')
c$$$      call prntm2(xyz1,3,1,3)
c$$$      call prntm2(xyz2,3,1,3)
c$$$      call prntm2(xyz3,3,1,3)
c$$$      call prntm2(xyz4,3,1,3)
      do i=1,4
       do j=1,4
        vmat(j,i)=1d-5
       end do
      end do
      do i=1,12
       dpart(i,1)=0d0
       dpart(i,2)=0d0
       dpart(i,3)=0d0
       dpart(i,4)=0d0
      end do
      call vtdd(xyz1,xyz2,xyz3,xyz4,vmat(1,1),dpart)
      call vn2n2d(xyz1,xyz2,xyz3,xyz4,vmat(2,2),dpart(1,2),shiftof0,0)  9d3s09
      call vn2n2d(xyz1,xyz3,xyz2,xyz4,vmat(3,3),dpart(1,3),shiftof0,0)  9d3s09
      call vn2n2d(xyz1,xyz4,xyz2,xyz3,vmat(4,4),dpart(1,4),shiftof0,0)  9d3s09
      do i=1,3
       sv=dpart(i+3,3)
       dpart(i+3,3)=dpart(i+6,3)
       dpart(i+6,3)=sv
       sv=dpart(i+3,4)
       dpart(i+3,4)=dpart(i+6,4)
       dpart(i+6,4)=dpart(i+9,4)
       dpart(i+9,4)=sv
      end do
c$$$      write(6,3351)(vmat(i,i),i=1,4),shiftof0
 3351 format(5f10.5)
      vmat(1,1)=vmat(1,1)-shiftof0                                      9d3s09
      do i=2,4
       vmat(1,i)=6.0847866d-02
       vmat(i,1)=6.0847866d-02
      end do
c$$$      call rs(4,4,vmat,eig,0,vmat,fv1,fv2,ierr)
c$$$      fit=eig(1)
      zero=0d0
c$$$      if(idoit.le.4)then
c$$$       eig(1)=vmat(idoit,idoit)
c$$$       do i=1,4
c$$$        vec(i)=0d0
c$$$       end do
c$$$       vec(idoit)=1d0
c$$$      else
   
      ic = 0
      do i=1,4
       do j=1,4
          ic = ic + 1
          vmatc(j,i)=vmat(j,i)
       end do
      end do
     
      
      !     
      ! Call to MKL (Find the eigenvalues and eigenvectors)
      ! ---------------------------------------------------
      ! Compute only the 1st eigenvalue   ! MARCO

!     call dsyevx('V','I','L',4,vmat,4,vl,ul,1,1,zero,neig, eig,vec,4, 
!     &            work,lwork,iwork,ifail,info)

!     LAPACK  dsyevx (JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,    &
!                     ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, & 
!                     INFO)      
!     ---------------------------------------------------------------
        call dsyevx('V','I','L',4,vmat,4,vl,ul,1,1,zero,neig,
     $     eig,vec,4,work,lwork,iwork,ifail,info)


c$$$  call prntm2(vec,1,4,1)
      do i=1,4
         vec(i)=vec(i)**2
      end do
      if(info.ne.0)then
         write(6,*)('on return from dsyevx in n4fitd ')
         write(6,*)('info = '),info
         stop
      end if
c$$$      end if
      fit=eig(1)
c$$$      write(6,*)('vectors ')
c$$$      call prntm2(vec,1,4,1)
c$$$      if(vec(1)**2.gt.0.5d0)then
c$$$       write(6,*)('we are in Td part of pes '),fit
c$$$       write(6,*)('carts: ')
c$$$       call prntm2(xyz1,3,1,3)
c$$$       call prntm2(xyz2,3,1,3)
c$$$       call prntm2(xyz3,3,1,3)
c$$$       call prntm2(xyz4,3,1,3)
c$$$       call vtdd(xyz1,xyz2,xyz3,xyz4,vmat(1,1),dpart)
c$$$       call vn2n2d(xyz1,xyz2,xyz3,xyz4,vmat(2,2),dpart(1,2),shiftof0,0)
c$$$       call vn2n2d(xyz1,xyz3,xyz2,xyz4,vmat(3,3),dpart(1,3),shiftof0,0)
c$$$       call vn2n2d(xyz1,xyz4,xyz2,xyz3,vmat(4,4),dpart(1,4),shiftof0,0)
c$$$       vmat(1,1)=vmat(1,1)-shiftof0
c$$$       do i=1,4
c$$$        write(6,351)i,vmat(i,i),vec(i)
c$$$  351   format(i5,1p2e15.7)
c$$$       end do
c$$$c$$$       stop
c$$$      end if
c$$$      dpart(1,2)=sqrt(-32d0)
      do i=1,12
       dfit(i)=dpart(i,1)*vec(1)+dpart(i,2)*vec(2)
     $      +dpart(i,3)*vec(3)+dpart(i,4)*vec(4)
       if(dfit(i).ne.dfit(i).or.dfit(i).gt.1d20.or.dfit(i).lt.-1d20)then
        write(6,*)('got nan for dfit '),i
        write(6,*)('dparts '),dpart(i,1),dpart(i,2),
     $       dpart(i,3),dpart(i,4)
        write(6,*)('vecs '),vec(1),vec(2),vec(3),vec(4)
        write(6,*)('vmat: ')
        do k=1,4
         write(6,*)(vmat(k,j),j=1,4)
        end do
        if(dpart(i,2).ne.dpart(i,2).or.dpart(i,2).gt.1d20.or.
     $       dpart(i,2).le.-1d20)then
           call vn2n2d(xyz1,xyz2,xyz3,xyz4,vmat(2,2),dpart(1,2),
     $          shiftof0,1)
        else if(dpart(i,3).ne.dpart(i,3).or.dpart(i,3).gt.1d20.or.
     $          dpart(i,3).le.-1d20)then
           call vn2n2d(xyz1,xyz3,xyz2,xyz4,vmat(3,3),dpart(1,3),
     $          shiftof0,1)
        else if(dpart(i,4).ne.dpart(i,4).or.dpart(i,4).gt.1d20.or.
     $          dpart(i,4).lt.-1d20)then
           call vn2n2d(xyz1,xyz4,xyz2,xyz3,vmat(4,4),dpart(1,4),
     $          shiftof0,1)
        end if
        stop
       end if
c$$$       write(6,51)i,dfit(i),(dpart(i,j),vec(j),j=1,4)
c$$$   51  format(i5,1p9e15.7)
      end do
c$$$      rms=0d0                                                           9d8s09
c$$$      do i=1,12                                                         9d8s09
c$$$       rms=rms+(dfit(i)-dfit2(i))**2                                    9d10s09
c$$$      end do                                                            9d8s09
c$$$      rms=sqrt(rms/12d0)                                                9d8s09
c$$$c$$$      write(6,*)('der rms: '),rms
c$$$      if(.not.(rms.lt.1d-6))then                                        9d8s09
c$$$       write(6,*)('numerical and analytic ders do not agree: '),rms     9d8s09
c$$$       do i=1,12                                                        9d8s09
c$$$        write(6,3451)i,cart(i,1),dfit2(i),dfit(i)                       9d10s09
c$$$ 3451   format(i5,1p21e15.7)                                            9d8s09
c$$$       end do                                                           9d8s09
c$$$       stop                                                             9d8s09
c$$$      end if                                                            9d8s09
      if(dfit(1).ge.0d0)then
       iok=1
      else
       if(dfit(1).lt.0d0)then
        iok=1
       else
        write(6,*)('have not a number for first der: ')
        do i=1,4
         write(6,*)i,vec(i),dpart(1,i)
         if(dpart(1,i).le.0d0.or.dpart(1,i).gt.0d0)then
          iok=iok+1
         else
          if(i.eq.2)then
           call vn2n2d(xyz1,xyz2,xyz3,xyz4,vmat(2,2),dpart(1,2),
     $          shiftof0,1)
          else if(i.eq.3)then
           call vn2n2d(xyz1,xyz3,xyz2,xyz4,vmat(3,3),dpart(1,3),
     $          shiftof0,1)
          else if(i.eq.4)then
           call vn2n2d(xyz1,xyz4,xyz2,xyz3,vmat(4,4),dpart(1,4),
     $          shiftof0,1)
          end if
         end if
        end do
        stop
       end if
      end if
      return
      end
      subroutine vtdd(x1,x2,x3,x4,fit,dfit)
      implicit real*8 (a-h,o-z)
c     compute N4 Td potential and derivatives
      dimension x1(3),x2(3),x3(3),x4(3), dfit(12)
       dimension dvx12(12),dvx13(12),dvx14(12), dvx23(12),dvx24(12),
     b dvx34(12), dr12(12),dr13(12),dr14(12),dr23(12),dr24(12),dr34(12)
           save
      data dr12,dr13,dr14,dr23,dr24,dr34 / 72*0.d0 /
      data dvx12,dvx13,dvx14,dvx23,dvx24,dvx34 / 72*0.d0 /
      data retd,betatd,detd /2.76149524d0,1.075722195d0,0.0955135476d0/
           data ifirst / 0 /
           if(ifirst.eq.0) then
   db2td = -2d0*detd*betatd
           end if
c23456789012345678901234567890123456789012345678901234567890123456789012
      r12=0d0
      r13=0d0
      r14=0d0
      r23=0d0
      r24=0d0
      r34=0d0
      do i=1,3
       r12=r12+(x1(i)-x2(i))**2
       r13=r13+(x1(i)-x3(i))**2
       r14=r14+(x1(i)-x4(i))**2
       r23=r23+(x2(i)-x3(i))**2
       r24=r24+(x2(i)-x4(i))**2
       r34=r34+(x3(i)-x4(i))**2
      end do
      r12=sqrt(r12)
      r13=sqrt(r13)
      r14=sqrt(r14)
      r23=sqrt(r23)
      r24=sqrt(r24)
      r34=sqrt(r34)
c
          do i=1,3
      dr12(i) = (x1(i)-x2(i))/r12
      dr13(i) = (x1(i)-x3(i))/r13
      dr14(i) = (x1(i)-x4(i))/r14
      dr23(i+3) = (x2(i)-x3(i))/r23
      dr24(i+3) = (x2(i)-x4(i))/r24
      dr34(i+6) = (x3(i)-x4(i))/r34
      dr12(i+3) = -dr12(i)
      dr13(i+6) = -dr13(i)
      dr14(i+9) = -dr14(i)
      dr23(i+6) = -dr23(i+3)
      dr24(i+9) = -dr24(i+3)
      dr34(i+9) = -dr34(i+6)
      end do
c
cdbg           write(6,1) r12,r13,r14,r23,r24,r34
cdbg  1        format(5x,'in VTD', 3x,'6 r(a,b)',6f15.5/)
cdbg           do k=1,12
cdbg           write(6,2) k,dr12(k),dr13(k),dr14(k),dr23(k),dr24(k),dr34(k) 
cdbg  2        format(16x,'drdx',i2,6e15.5)
cdbg           end do
      
ccccc      retd=2.76149524d0
ccccc      betatd=1.075722195d0
ccccc      detd=0.0955135476d0
      ex12=exp(-betatd*(r12-retd))
      vx12=((ex12-1d0)**2)*detd
      ex13=exp(-betatd*(r13-retd))
      vx13=((ex13-1d0)**2)*detd
      ex14=exp(-betatd*(r14-retd))
      vx14=((ex14-1d0)**2)*detd
      ex23=exp(-betatd*(r23-retd))
      vx23=((ex23-1d0)**2)*detd
      ex24=exp(-betatd*(r24-retd))
      vx24=((ex24-1d0)**2)*detd
      ex34=exp(-betatd*(r34-retd))
      vx34=((ex34-1d0)**2)*detd
      fit=vx12+vx13+vx14+vx23+vx24+vx34
c     from my ccsd(t) caln
     $-218.47601249d0     
c     empirical shift to get X peak height about right
     $+0.10d0
c     compute derivatives wrt atomic cartesian coordinates
      dvx12dr = db2td*ex12*(ex12-1d0)
      dvx13dr = db2td*ex13*(ex13-1d0)      
    dvx14dr = db2td*ex14*(ex14-1d0)      
    dvx23dr = db2td*ex23*(ex23-1d0)      
    dvx24dr = db2td*ex24*(ex24-1d0)      
    dvx34dr = db2td*ex34*(ex34-1d0)   
cdbg          write(6,3) dvx12dr,dvx13dr,dvx14dr,dvx23dr,dvx24dr,dvx34dr
cdbg  3       format(5x,'dVTDdr',5x,6e15.5)
         do i=1,3
     dvx12(i) = dvx12dr*dr12(i)
     dvx13(i) = dvx13dr*dr13(i)
     dvx14(i) = dvx14dr*dr14(i)
     dvx23(i+3) = dvx23dr*dr23(i+3)
     dvx24(i+3) = dvx24dr*dr24(i+3)
     dvx34(i+6) = dvx34dr*dr34(i+6)   
     dvx12(i+3) = dvx12dr*dr12(i+3)
     dvx13(i+6) = dvx13dr*dr13(i+6)
     dvx14(i+9) = dvx14dr*dr14(i+9)
     dvx23(i+6) = dvx23dr*dr23(i+6)
     dvx24(i+9) = dvx24dr*dr24(i+9)
     dvx34(i+9) = dvx34dr*dr34(i+9)   
         end do
     do i=1,12
  dfit(i) = dvx12(i)+dvx13(i)+dvx14(i)+dvx23(i)+dvx24(i)+dvx34(i)
     end do     
cdbg          write(6,4) (dfit(k),k=1,12)
cdbg  4       format(5x,'dVTDdx',5x,6e15.5/15x,6e15.5)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine vn2n2d(x1,x2,x3,x4,fit,dfit,shiftof0,iprt)
      implicit real*8 (a-h,o-z)
      dimension x1(3),x2(3),x3(3),x4(3),xcm(3),rv(3,2),
     $     coef(20), rab(3),rcd(3)
      dimension dr1(12),dr2(12),dr3(12),dr4(12),dr5(12),dr6(12)
      dimension drv(3,2,12),drab(3,12),drcd(3,12),dxcm(3,12), 
     2     x12(3),x34(3), drbig(12), dfit(12), 
     3     dtheta1(12),dtheta2(12),dphi(12), dxa(12),dxb(12),ddot3(12)
      dimension dxa1(12),dxa2(12),dxa3(12), dxb1(12),dxb2(12),dxb3(12)
      save
      data ifirst/0/
      data x6,f2,f220,f221,f222/
     $   75.63d0,0.129334171818588d0,5.220578416136951d-002,
     $   -1.159697407223882d-002,1.450724604023408d-003/
      data x8,g2,g220,g221,g222/
     $     2489.d0,0.624578507475241d0,0.157105504988713d0,
     $     -2.188839577146523d-002,1.003716064096833d-003/
      data dr1,dr2,dr3,dr4,dr5,dr6 / 72*0.d0 / 
      data drab,drcd,dxcm,drv  / 72*0d0, 36*0.d0, 72*0.d0 /
      data dxa1,dxa2,dxa3, dxb1,dxb2,dxb3 / 72*0.d0 /
      data dxa,dxb /24*0.d0 /
c23456789012345678901234567890123456789012345678901234567890123456789012
      if(iprt.ne.0)then
       write(6,*)('input carts: ')
       do i=1,3
        write(6,5115)x1(i),x2(i),x3(i),x4(i)
       end do
      end if
      if(ifirst.eq.0)then
       pi=acos(-1d0)                                                    9d4s09
       bswth=4d0
c        eps added to rbig in case rbig=0   RLJ 8/09
       eps = 1.0d-6
       dpmn=-14.9681846800000d0
       damp=6d0
       damp6=damp**6
       damp4=damp**4
       damp5=damp**5
       beta=-2.18d0
       beta70 = 7d1*beta
       re = 2.1d0
       call diatccsdd(re,veccsd,dum,1)                                  9d3s09
       write(6,*)('ccsd potential at "re": '),veccsd                    9d3s09
       call vn2d(re,veuse,dum)                                          9d3s09
       write(6,*)('diatomic potential at "re": '),veuse                 9d3s09
       shiftof0=veccsd-veuse                                            9d3s09
       write(6,*)('shift of diatomic zero of energy: '),shiftof0        9d3s09
       shiftof0=2d0*shiftof0                                            9d3s09
       write(6,*)('shift of VTd '),shiftof0                             9d3s09
       re2 = 2d0*re
       coef(1)=1.0345018D+00
       coef(2)=9.7550153D+00
       coef(3)=2.3947270D-01
       coef(4)=-4.4270639D-01
       coef(5)=4.2206130D+00
       coef(6)=-9.6318443D+00
       coef(7)=-1.7070046D+01
       coef(8)=1.5595772D+00
       coef(9)=9.0354414D-01
       coef(10)=-3.7645959D+00
       coef(11)=6.2330259D-01
       coef(12)=1.5681477D+00
       coef(13)=8.1289901D-01
       coef(14)=-2.7320998D-01
       coef(15)=8.4521948D+00
       coef(16)=5.7363853D+00
       coef(17)=-1.3138353D+00
       coef(18)=-7.5160389D-01
       coef(19)=1.3283815D-01
       coef(20)=0.0038d0
c         for derivatives of r's and jacobi coordinates
       do i=1,3
       drab(i,i)     =  .5d0
       drab(i,i+3)   =  .5d0
       drcd(i,i+6)   =  .5d0
       drcd(i,i+9)   =  .5d0
       dxcm(i,i)     = -.5d0
       dxcm(i,i+3)   = -.5d0
       dxcm(i,i+6)   =  .5d0
       dxcm(i,i+9)   =  .5d0
       do j=1,6
     drv(i,1,j)   = -0.5d0
     drv(i,2,j+6) = -0.5d0
     end do
     drv(i,1,i)   = 0.5d0
     drv(i,2,i+6) = 0.5d0
       end do
c         for derivatives of the potential
       vibdispt = 2d0*sqrt(abs(dpmn))/(dpmn*dpmn)
       dp1 = 873.238d0*exp(-2.2d0)*2.1d0 - 42.6175d0
       dp2 = 873.238d0*exp(-2.2d0)*2.1d0 - 42.6175d0
       preqq = 1.5d0*sqrt(2d0)
       ifirst=1
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      r1=0d0
      r5=0d0
      r3=0d0
      r6=0d0
      r4=0d0
      r2=0d0
      do i=1,3
          x12(i) = x1(i)-x2(i)
          x34(i) = x3(i)-x4(i)
c
ccc       r1=r1+(x1(i)-x2(i))**2
       r1=r1+x12(i)**2
       r5=r5+(x1(i)-x3(i))**2
       r3=r3+(x1(i)-x4(i))**2
       r6=r6+(x2(i)-x3(i))**2
       r4=r4+(x2(i)-x4(i))**2
ccc       r2=r2+(x3(i)-x4(i))**2
       r2=r2+x34(i)**2
      end do
      r1=sqrt(r1)
      r5=sqrt(r5)
      r3=sqrt(r3)
      r6=sqrt(r6)
      r4=sqrt(r4)
      r2=sqrt(r2)
      dot1=0d0
      dot2=0d0
      rbig=0d0
      do j=1,3
       rab(j)=0.5d0*(x1(j)+x2(j))
       rcd(j)=0.5d0*(x3(j)+x4(j))
       xcm(j)=rcd(j)-rab(j)
       rv(j,1)=x1(j)-rab(j)
       rv(j,2)=x3(j)-rcd(j)
       rbig=rbig+(rcd(j)-rab(j))**2
       dot1=dot1+(x1(j)-x2(j))*(rcd(j)-rab(j))
       dot2=dot2+(x3(j)-x4(j))*(rcd(j)-rab(j))
      end do
      rbig=sqrt(rbig)
c$$$      theta1=acos(dot1/(r1*rbig))                                       6d20s14
c$$$      theta2=acos(dot2/(r2*rbig))
      theta1arg=dot1/(r1*rbig)                                          6d20s14
      if(abs(theta1arg).gt.1d0)then                                     6d20s14
       if(theta1arg.gt.0d0)then                                         6d20s14
        theta1arg=1d0                                                   6d20s14
       else                                                             6d20s14
        theta1arg=-1d0                                                  6d20s14
       end if                                                           6d20s14
      end if                                                            6d20s14
      theta1=acos(theta1arg)                                            6d20s14
      if(theta1.ne.theta1)then
       write(6,*)('theta1 nan '),dot1,r1,rbig,theta1arg
c$$$       stop
      end if
      theta2arg=dot2/(r2*rbig)                                          6d20s14
      if(abs(theta2arg).gt.1d0)then                                     6d20s14
       if(theta2arg.gt.0d0)then                                         6d20s14
        theta2arg=1d0                                                   6d20s14
       else                                                             6d20s14
        theta2arg=-1d0                                                  6d20s14
       end if                                                           6d20s14
      end if                                                            6d20s14
      theta2=acos(theta2arg)                                            6d20s14
      xa1=rv(2,1)*xcm(3)-rv(3,1)*xcm(2)
      xa2=rv(3,1)*xcm(1)-rv(1,1)*xcm(3)
      xa3=rv(1,1)*xcm(2)-rv(2,1)*xcm(1)
      xb1=rv(2,2)*xcm(3)-rv(3,2)*xcm(2)
      xb2=rv(3,2)*xcm(1)-rv(1,2)*xcm(3)
      xb3=rv(1,2)*xcm(2)-rv(2,2)*xcm(1)
      xa=sqrt(xa1**2+xa2**2+xa3**2)
      xb=sqrt(xb1**2+xb2**2+xb3**2)
      rxaxb = 1.d0/(xa*xb)
      if(abs(xa).le.1d-6.or.abs(xb).le.1d-6)then
       phi=0d0
       tphi = 0.d0
       dot3 = 0.d0
      else
       dot3=(xa1*xb1+xa2*xb2+xa3*xb3)*rxaxb
     tphi = 1.d0/(sqrt(abs(1.d0-dot3**2))+1d-14)                         9d4s09
       if(iprt.ne.0)then
        write(6,*)('dot3 = '),dot3,xa1,xb1,xa2,xb2,xa3,xb3,rxaxb
        write(6,*)('dot3-1d0 '),dot3-1d0,abs(abs(dot3)-1d0)
       end if
       ipst=0
       if(abs(abs(dot3)-1d0).lt.1d-10)then                              9d4s09
        if(dot3.gt.0d0)then                                             9d4s09
         phi=0d0                                                        9d4s09
         ipst=1
        else                                                            9d4s09
         phi=pi                                                         9d4s09
         ipst=2
        end if                                                          9d4s09
       else                                                             9d4s09
        phi=acos(dot3)
        ipst=3
       end if                                                           9d4s09
       if(iprt.gt.0)then
        write(6,*)('phi, ipst '),phi,ipst,pi
       end if
      end if
c$$$         write (6,1) r1,r2,r3,r4,r5,r6,rbig,theta1,theta2,phia,
c$$$     2     tphi,dot3
  1   format(5x,'r(i) & jacobi coord',6f12.5/4f15.8,2e15.5/)
c     derivatives of r's and Jacobi coordinates wrt atomic cartesians    (RLJ)
      do i=1,3
      dr1(i) = x12(i)/r1
      dr1(i+3) = -dr1(i)
c$$$      dr2(i+6) = x34(i)/r1
      dr2(i+6) = x34(i)/r2                                              8d31s09
      dr2(i+9) = -dr2(i+6)
      dr3(i) = (x1(i)-x4(i))/r3
      dr3(i+9) = -dr3(i)
      dr4(i+3) = (x2(i)-x4(i))/r4
      dr4(i+9) = -dr4(i+3)
      dr5(i) = (x1(i)-x3(i))/r5
      dr5(i+6) = -dr5(i)
      dr6(i+3) = (x2(i)-x3(i))/r6
      dr6(i+6) = -dr6(i+3)
      drbig(i) = -.5*xcm(i)/rbig
      drbig(i+3) =  drbig(i)
      drbig(i+6) = -drbig(i)
      drbig(i+9) = -drbig(i)
c      ddot1(i)   =  xcm(i) -.5d0*x12(i)
c      ddot1(i+3) = -xcm(i) -.5d0*x12(i)
c      ddot1(i+6) =          .5d0*x12(i)
c      ddot1(i+9) =          .5d0*x12(i)
c      ddot2(i)   =         -.5d0*x34(i)
c      ddot2(i+3) =         -.5d0*x34(i)
c      ddot2(i+3) =  xcm(i) +.5d0*x34(i)
c      ddot2(i+3) = -xcm(i) +.5d0*x34(i)
      end do
         trat1 = dot1/(r1*rbig)
         trat2 = dot2/(r2*rbig)
         tpre1 = 1.d0/sqrt(1.d0 - trat1**2)
         xpre1 = 1.d0/(r1*rbig)
         ypre1 = .5d0*dot1/(r1*rbig**2)
         zpre1 = dot1/(rbig*r1**2)
         tpre2 = 1.d0/sqrt(1.d0 - trat2**2)
         xpre2 = 1.d0/(r2*rbig)
         ypre2 = .5d0*dot2/(r2*rbig**2)
         zpre2 = dot2/(rbig*r2**2)
      dws1=-1d0/sqrt(abs(r1*r1*rbig*rbig-dot1*dot1)+1d-14)              8d11s14
      dws2=-1d0/sqrt(abs(r2*r2*rbig*rbig-dot2*dot2)+1d-14)              8d11s14
      if(iprt.ne.0)then
      write(6,*)('dws1, dws2 '),dws1,dws2
      write(6,*)('sqrt arg '),r1*r1*rbig*rbig,dot1*dot1,
     $     r1*r1*rbig*rbig-dot1*dot1
      write(6,*)('r1,rbig '),r1,r2,rbig
      write(6,*)('dot1 '),dot1
      write(6,*)('rab '),rab
      write(6,*)('rcd '),rcd
      write(6,*)('x12 '),x12
      write(6,*)('dr1 '),dr1
      write(6,*)('drbig '),drbig
      end if
      do i=1,3
c23456789012345678901234567890123456789012345678901234567890123456789012
c$$$      dtheta1(i) = trat1*(xpre1*(xcm(i) -.5d0*x12(i)) + 
c$$$     b    ypre1*xcm(i)/rbig -zpre1*x12(i)/r1)
c$$$      dtheta1(i+3) = trat1*(xpre1*(xcm(i) -.5d0*x12(i)) + 
c$$$     b    ypre1*xcm(i)/rbig +zpre1*x12(i)/r1)
c$$$      dtheta1(i+6) = trat1*(xpre1*.5d0*x12(i) - ypre1*xcm(i)/rbig) 
c$$$      dtheta1(i+9) = dtheta1(i+6)
      dtheta1(i)=dws1*(rcd(i)-rab(i)-0.5d0*x12(i)-dot1*((dr1(i)/r1)     8d31s09
     $     +(drbig(i)/rbig)))                                           8d31s09
      dtheta1(i+3)=dws1*(rab(i)-rcd(i)-0.5d0*x12(i)-dot1*               8d31s09
     $     ((dr1(i+3)/r1)+(drbig(i+3)/rbig)))                           8d31s09
      dtheta1(i+6)=dws1*(0.5d0*x12(i)-dot1*((dr1(i+6)/r1)               8d31s09
     $     +(drbig(i+6)/rbig)))                                         8d31s09
      dtheta1(i+9)=dws1*(0.5d0*x12(i)-dot1*((dr1(i+9)/r1)               8d31s09
     $     +(drbig(i+9)/rbig)))                                         8d31s09
c
c$$$      dtheta2(i+6) = trat2*(xpre2*(xcm(i) -.5d0*x34(i)) + 
c$$$     b    ypre2*xcm(i)/rbig -zpre2*x34(i)/r2)
c$$$      dtheta2(i+9) = trat2*(xpre2*(xcm(i) -.5d0*x34(i)) + 
c$$$     b    ypre2*xcm(i)/rbig +zpre2*x34(i)/r2)
c$$$      dtheta2(i)   = trat2*(xpre2*.5d0*x34(i) - ypre2*xcm(i)/rbig) 
c$$$      dtheta2(i+3) = dtheta2(i)
      dtheta2(i)=dws2*(-0.5d0*x34(i)-dot2*((dr2(i)/r2)                  8d31s09
     $     +(drbig(i)/rbig)))                                           8d31s09
      dtheta2(i+3)=dws2*(-0.5d0*x34(i)-dot2*                            8d31s09
     $     ((dr2(i+3)/r2)+(drbig(i+3)/rbig)))                           8d31s09
      dtheta2(i+6)=dws2*(rcd(i)-rab(i)+0.5d0*x34(i)-dot2*((dr2(i+6)/r2) 8d31s09
     $     +(drbig(i+6)/rbig)))                                         8d31s09
      dtheta2(i+9)=dws2*(rab(i)-rcd(i)+0.5d0*x34(i)-dot2*((dr2(i+9)/r2) 8d31s09
     $     +(drbig(i+9)/rbig)))                                         8d31s09
      end do
c     derivative of phi
      dxa1(2) =  .5d0*(rv(3,1) + xcm(3))
      dxa1(3) = -.5d0*(rv(2,1) + xcm(2)) 
      dxa1(5) =  .5d0*(rv(3,1) - xcm(3)) 
      dxa1(6) = -.5d0*(rv(2,1) - xcm(2))    
      dxa1(8) = -.5d0*rv(3,1)
      dxa1(9) =  .5d0*rv(2,1)
      dxa1(11) = -.5d0*rv(3,1)
      dxa1(12) =  .5d0*rv(2,1)
      dxb1(2) =  .5d0*rv(3,2)
      dxb1(3) = -.5d0*rv(2,2)
      dxb1(5) =  .5d0*rv(3,2)
      dxb1(6) = -.5d0*rv(2,2)   
      dxb1(8) = -.5d0*(rv(3,2) - xcm(3))
      dxb1(9) =  .5d0*(rv(2,2) - xcm(2)) 
      dxb1(11) = -.5d0*(rv(3,2) + xcm(3)) 
      dxb1(12) =  .5d0*(rv(2,2) + xcm(2))    
c  
      dxa2(1) = -.5d0*(rv(3,1) + xcm(3))
      dxa2(3) =  .5d0*(rv(1,1) + xcm(1)) 
      dxa2(4) = -.5d0*(rv(3,1) - xcm(3)) 
      dxa2(6) =  .5d0*(rv(1,1) - xcm(1))    
      dxa2(7) =  .5d0*rv(3,1)
      dxa2(9) = -.5d0*rv(1,1)
      dxa2(10) =  .5d0*rv(3,1)
      dxa2(12) = -.5d0*rv(1,1)
      dxb2(1) = -.5d0*rv(3,2)
      dxb2(3) =  .5d0*rv(1,2)
      dxb2(4) = -.5d0*rv(3,2)
      dxb2(6) =  .5d0*rv(1,2)   
      dxb2(7) =  .5d0*(rv(3,2) - xcm(3))
      dxb2(9) = -.5d0*(rv(1,2) - xcm(1)) 
      dxb2(10) =  .5d0*(rv(3,2) + xcm(3)) 
cccc      dxb2(12) = -.5d0*(rv(1,2) - xcm(1))   fixed 8/31 RLJ 
      dxb2(12) = -.5d0*(rv(1,2) + xcm(1))    
c
      dxa3(1) =  .5d0*(rv(2,1) + xcm(2))
      dxa3(2) = -.5d0*(rv(1,1) + xcm(1)) 
      dxa3(4) =  .5d0*(rv(2,1) - xcm(2)) 
      dxa3(5) = -.5d0*(rv(1,1) - xcm(1))    
      dxa3(7) = -.5d0*rv(2,1)
      dxa3(8) =  .5d0*rv(1,1)
      dxa3(10) = -.5d0*rv(2,1)
      dxa3(11) =  .5d0*rv(1,1)
      dxb3(1) =  .5d0*rv(2,2)
      dxb3(2) = -.5d0*rv(1,2)
      dxb3(4) =  .5d0*rv(2,2)
      dxb3(5) = -.5d0*rv(1,2)   
      dxb3(7) = -.5d0*(rv(2,2) - xcm(2))
      dxb3(8) =  .5d0*(rv(1,2) - xcm(1)) 
      dxb3(10) = -.5d0*(rv(2,2) + xcm(2)) 
      dxb3(11) =  .5d0*(rv(1,2) + xcm(1))    
      rxa1 = xa1/xa
      rxa2 = xa2/xa
      rxa3 = xa3/xa
      rxb1 = xb1/xb
      rxb2 = xb2/xb
      rxb3 = xb3/xb
      dxa(1) = rxa2*dxa2(1) + rxa3*dxa3(1)
      dxa(4) = rxa2*dxa2(4) + rxa3*dxa3(4)
      dxa(7) = rxa2*dxa2(7) + rxa3*dxa3(7)
      dxa(10) = rxa2*dxa2(10) + rxa3*dxa3(10)
      dxa(2) = rxa1*dxa1(2) + rxa3*dxa3(2)
      dxa(5) = rxa1*dxa1(5) + rxa3*dxa3(5)
      dxa(8) = rxa1*dxa1(8) + rxa3*dxa3(8)
      dxa(11) = rxa1*dxa1(11) + rxa3*dxa3(11)
      dxa(3) = rxa1*dxa1(3) + rxa2*dxa2(3)
      dxa(6) = rxa1*dxa1(6) + rxa2*dxa2(6)
      dxa(9) = rxa1*dxa1(9) + rxa2*dxa2(9)
      dxa(12) = rxa1*dxa1(12) + rxa2*dxa2(12)
      dxb(1) = rxb2*dxb2(1) + rxb3*dxb3(1)
      dxb(4) = rxb2*dxb2(4) + rxb3*dxb3(4)
      dxb(7) = rxb2*dxb2(7) + rxb3*dxb3(7)
      dxb(10) = rxb2*dxb2(10) + rxb3*dxb3(10)
      dxb(2) = rxb1*dxb1(2) + rxb3*dxb3(2)
      dxb(5) = rxb1*dxb1(5) + rxb3*dxb3(5)
      dxb(8) = rxb1*dxb1(8) + rxb3*dxb3(8)
      dxb(11) = rxb1*dxb1(11) + rxb3*dxb3(11)
      dxb(3) = rxb1*dxb1(3) + rxb2*dxb2(3)
      dxb(6) = rxb1*dxb1(6) + rxb2*dxb2(6)
      dxb(9) = rxb1*dxb1(9) + rxb2*dxb2(9)
      dxb(12) = rxb1*dxb1(12) + rxb2*dxb2(12)
      do i=1,12
      ddot3(i) = (xa1*dxb1(i)+xb1*dxa1(i) + xa2*dxb2(i)+xb2*dxa2(i) + 
     b xa3*dxb3(i)+xb3*dxa3(i))*rxaxb -dot3*rxaxb*(xa*dxb(i)+xb*dxa(i))
c$$$      dphi(i) = tphi*ddot3(i)
      dphi(i) = -tphi*ddot3(i)                                          8d31s09
      end do
      if(iprt.ne.0)then
       write(6,*)('tphi '),tphi
       do i=1,12
        write(6,5115)ddot3(i)
       end do
      end if
cdbg       write(6,19) (k,dxa1(k),dxa2(k),dxa3(k),dxa1(k),dxa2(k),dxa3(k), 
cdbg     2  ddot3(k),k=1,12)
 19    format('  debug ddot3',i5,7e14.4)
c$$$         do k=1,12
c$$$         write (6,2) dr1(k),dr2(k),dr3(k),dr4(k),dr5(k),dr6(k),
c$$$     2      drbig(k),dtheta1(k),dtheta2(k),dphi(k)
c$$$         end do
  2   format(5x,'dr(i).dx(j) & deriv. of jacobi coord',10e12.4)
c     get diatomic potential curves
      if(iprt.ne.0)then                                                 9d3s09
      write(6,515)r1,r2,rbig,theta1,theta2,phi                          9d3s09
  515 format('jacobis: ',6f10.5)                                        9d3s09
      write(6,*)('jacobi ders ')
      do i=1,12
       write(6,5115)dr1(i),dr2(i),drbig(i),dtheta1(i),dtheta2(i),dphi(i)
 5115  format(1p6e26.18)
      end do
      write(6,*)('other rij ders ')
      do i=1,12
       write(6,5115)dr3(i),dr4(i),dr5(i),dr6(i)
      end do
      end if                                                            9d3s09
c$$$      call diatccsdd(r1,v12,dv12dr1,1)                                           7d27s09
c$$$      call diatccsdd(r2,v34,dv34dr2,1)                                           7d27s09
      call vn2d(r1,v12,dv12dr1)
      call vn2d(r2,v34,dv34dr2)
      call diatccsdd(r1,dp1,ddp1dr1,3)
      call diatccsdd(r2,dp2,ddp2dr2,3)
      call diatccsdd(r1,q1,dq1dr1,2)
      call diatccsdd(r2,q2,dq2dr2,2)
cdbg         write(6,3) v12,v34,q1,q2,dv12dr1,dv34dr2,dq1dr1,dq2dr2
  3   format(5x,'v12,v34,q1,q2',4f15.7/13x,'deriv',4e15.5)
      vibdispa=dp1*dp2/(dpmn*dpmn)
      vibdispb=vibdispa*sqrt(abs(dpmn))*2d0
     $     /(sqrt(abs(dp1))+sqrt(abs(dp2)))
c$$$      vibdispb=vibdispt*(dp1*dp2)/(sqrt(abs(dp1))+sqrt(abs(dp2)))
c$$$      dvibdisp1 =  0.d0                                                 (d vibdisp/d r1)
c$$$      dvibdisp2 =  0.d0                                                 (d vibdisp/d r2)
      dvibdisp1 =ddp1dr1*vibdispb*(1d0-0.5d0*sqrt(abs(dp1))*bot)/dp1    2d9s10
      dvibdisp2 =ddp2dr2*vibdispb*(1d0-0.5d0*sqrt(abs(dp2))*bot)/dp2    2d9s10
      c6000p=x6*vibdispb
      c8000p=x8*vibdispb
      q12 = q1*q2
      ct1=cos(theta1)
      ct2=cos(theta2)
      st1=sin(theta1)
      st2=sin(theta2)
      cph = cos(phi)
      sph = sin(phi)
cdbg       write(6,31) vibdispb,c6000p,c8000p, ct1,ct2,st1,st2,cph,sph
 31    format(5x,'vibdispb,c6000(pre),c8000(pre)',3e15.5/
     2   5x,'ct1,ct2,st1,st2,cph,sph',6f15.8/)
      cp=cos(2d0*phi)
      sp=sin(2d0*phi)
      cst1 = ct1*st1
      cst2 = ct2*st2
      y2a=ct1*ct1
      y2b=ct2*ct2
      z2a=st1*st1
      z2b=st2*st2
      y2=y2a+y2b
cccc      y3=ct1*ct1*ct2*ct2
      y3=y2a*y2b
cccc      y4=st1*st1*st2*st2*cp
      y4=z2a*z2b*cp 
      y5a=y2a*y2a
      y5b=y2b*y2b
      y5=y5a+y5b
cccc      y6a=st1*st1*cp                                                    6d9s09
cccc      y6b=st2*st2*cp                                                    6d9s09
      y6a=z2a*cp                                                            6d9s09
      y6b=z2b*cp                                                            6d9s09
      y6=y6a+y6b
cdbg         write(6,34) cst1,cst2,cp,sp,y2a,y2b,z2a,z2b, y2,y3,y4,y5,y6
  34     format(5x,'cst1,cst2,cp,sp,y2a,etc',8f10.6/5x,'y2..y6',5e15.5/)
cccc      p2a=0.5d0*(3d0*(ct1**2)-1d0)
cccc      p2b=0.5d0*(3d0*(ct2**2)-1d0)
      p2a=0.5d0*(3d0*y2a-1d0)
      p2b=0.5d0*(3d0*y2b-1d0)
cccc      p21a=-3d0*st1*ct1
cccc      p21b=-3d0*st2*ct2
      p21a=-3d0*cst1 
      p21b=-3d0*cst2 
cccc      p22a=3d0*(st1**2)
cccc      p22b=3d0*(st2**2)
      p22a=3d0*z2a
      p22b=3d0*z2b
c
      c6000=-c6000p*(1d0+f2*(p2a+p2b)+f220*p2a*p2b
     $     +f221*p21a*p21b+f222*p22a*p22b)
      dc6000r1=-dvibdisp1*x6*(1d0+f2*(p2a+p2b)+f220*p2a*p2b             2d9s10
     $     +f221*p21a*p21b+f222*p22a*p22b)                              2d9s10
      dc6000r2=-dvibdisp2*x6*(1d0+f2*(p2a+p2b)+f220*p2a*p2b             2d9s10
     $     +f221*p21a*p21b+f222*p22a*p22b)                              2d9s10
c23456789012345678901234567890123456789012345678901234567890123456789012
      dc6000t1 = -c6000p*(cst1*(-3d0*f2-3d0*f220*p2b+
     2    +18d0*f222*z2b) + 9d0*f221*cst2*(y2a-z2a))
      dc6000t2 = -c6000p*(cst2*(-3d0*f2-3d0*f220*p2a+
     2    +18d0*f222*z2a) + 9d0*f221*cst1*(y2b-z2b))
      c8000=-c8000p*(1d0+g2*(p2a+p2b)+g220*p2a*p2b
     $     +g221*p21a*p21b+g222*p22a*p22b)
      dc8000r1=-dvibdisp1*x8*(1d0+g2*(p2a+p2b)+g220*p2a*p2b             2d9s10
     $     +g221*p21a*p21b+g222*p22a*p22b)                              2d9s10
      dc8000r2=-dvibdisp2*x8*(1d0+g2*(p2a+p2b)+g220*p2a*p2b             2d9s10
     $     +g221*p21a*p21b+g222*p22a*p22b)                              2d9s10
      dc8000t1 = -c8000p*(cst1*(-3d0*g2-3d0*g220*p2b+
     2    +18d0*g222*z2b) + 9d0*g221*cst2*(y2a-z2a))
      dc8000t2 = -c8000p*(cst2*(-3d0*g2-3d0*g220*p2a+
     2    +18d0*g222*z2a) + 9d0*g221*cst1*(y2b-z2b))
c
      qq=preqq*(0.5d0*p2a*p2b -(p21a*p21b/9d0)*cph
     $     +(p22a*p22b/72d0)*cp)*q12
      if(iprt.ne.0)then
       write(6,*)('q1,q2 '),qq,dq1dr1,dq2dr2,q1,q2
      end if
      if(qq.ne.0d0)then                                                 6d11s15
       dqqdr1 = qq*(dq1dr1/q1)                                          6d11s15
       dqqdr2 = qq*(dq2dr2/q2)                                          6d11s15
      else                                                              6d11s15
       dqqdr1=0d0                                                       6d11s15
       dqqdr2=0d0                                                       6d11s15
      end if                                                            6d11s15
      dqqdph = preqq*((p21a*p21b/9d0)*sph -(p22a*p22b/36.d0)*sp)*q12
      dqqdt1 = preqq*(-1.5d0*cst1*p2b -cst2*(y2a-z2a)*cph
     2  +0.25d0*cst1*z2b*cp)*q12
c$$$      dqqdt2 = preqq*(-1.5d0*cst2*p2b -cst1*(y2b-z2b)*cph
      dqqdt2 = preqq*(-1.5d0*cst2*p2a -cst1*(y2b-z2b)*cph               8d31s09
     2  +0.25d0*cst2*z2a*cp)*q12
cdbg       write(6,4) c6000,c8000,qq, dc6000t1,dc6000t2, dc8000t1,dc8000t2,
cdbg     2  dqqdr1,dqqdr2,dqqdt1,dqqdt2,dqqdph
  4    format(5x,'c6,c8,qq', 3e15.5/7x,'dc6,dc8', 4e15.5/
     2  11x,'dqq',5e15.5/)
       denom1=1.d0/((rbig**6)+damp6)
       denom2a=1.d0/((rbig**4)+damp4)
       denom2 = denom2a*denom2a
       denom3=1.d0/((rbig**5)+damp5)
cccc      vlr=(c6000/((rbig**6)+damp6))
cccc     $     +(c8000/(((rbig**4)+damp4)**2))
cccc     $     +(qq/((rbig**5)+damp5))
      vlr=(c6000*denom1) + (c8000*denom2) + (qq*denom3)
      ex1=exp(beta*r3)+exp(beta*r4)+exp(beta*r5)+exp(beta*r6)
      pairwise=ex1*7d1
      sargtyp=(r1+r2)/(rbig+eps)                                        8d3s09
      swth=1d0/(1d0+exp(bswth*(sargtyp-1.2d0)))                         8d3s09
      if(swth.gt.1d-12)then                                             9d3s09
       dswth=-swth*swth*bswth*exp(bswth*(sargtyp-1.2d0))                 8d31s09
      else                                                              9d3s09
       dswth=0d0                                                        9d3s09
      end if                                                            9d3s09
      if(iprt.ne.0)then
       write(6,*)('dswth '),swth,bswth,sargtyp
      end if
      dswthab=dswth/(rbig+eps)                                          8d31s09
      if(iprt.ne.0)then
       write(6,*)('dswthab '),dswthab,dswth,rbig,eps
      end if
      dswthc=-(r1+r2)*dswthab/(rbig+eps)                                8d31s09
      swth1 = bswth*swth*(1d0-swth)
      dswthdr1 = -swth1/(rbig+eps)
      dswthdr2 = -swth1/(rbig+eps)
      dswthdrb = swth1*(r1+r2)/(rbig+eps)**2
      vlro=vlr
      vlr=vlr*swth                                                      7d27s09
c
c$$$      dvlrdr1=vlro*dswthr1 
c$$$      dvlrdr2=vlro*dswthr2 
      dvlrdr1=vlro*dswthab                                              8d31s09
      dvlrdr2=vlro*dswthab                                              8d31s09
c$$$      dvlrdt1= swth*(dc6000dt1*denom1 + dc8000dt1*denom2
c$$$     2           + dqqdt1*denom3)
c$$$      dvlrdt2= swth*(dc6000dt2*denom1 + dc8000dt2*denom2
c$$$     2           + dqqdt2*denom3)
      dvlrdt1= swth*(dc6000t1*denom1 + dc8000t1*denom2                  8d31s09
     2           + dqqdt1*denom3)
      dvlrdt2= swth*(dc6000t2*denom1 + dc8000t2*denom2                  8d31s09
     2           + dqqdt2*denom3)
c$$$      dvlrdrb=vlro*dswthrb - swth*(6d0*c6000*(denom1**2)*(rbig**5)
c$$$     2 +8d0*c8000*denom2*denom2a*(rbig**3)+5d0*qq*(denom3**2)*(rbig**4))
      dvlrdrb=vlro*dswthc - swth*(6d0*c6000*(denom1**2)*(rbig**5)       8d31s09
     2 +8d0*c8000*denom2*denom2a*(rbig**3)+5d0*qq*(denom3**2)*(rbig**4))
c$$$      write(6,*)('dvlrdr1 = '),dvlrdr1,vlro,dswthr1
      if(iprt.ne.0)then
       write(6,*)('dvlrdr1 '),dvlrdr1,swth,dqqdr1,denom3
       write(6,*)('dvlrdr2 '),dvlrdr2,dqqdr2,dc6000r2,denom1,
     $      dc8000r2,denom2
      end if
      dvlrdr1=dvlrdr1+swth*(dqqdr1*denom3+dc6000r1*denom1               2d9s10
     $     +dc8000r1*denom2)                                            2d9s10
      dvlrdr2=dvlrdr2+swth*(dqqdr2*denom3+dc6000r2*denom1               2d9s10
     $     +dc8000r2*denom2)                                            2d9s10
      dvlrph=denom3*dqqdph*swth                                         8d31s09
c$$$      do i=1,12
c     ok
c$$$       dfit(i)=(dvlrdr1+swth*dqqdr1*denom3)*dr1(i)
c$$$     $        +(dvlrdr2+swth*dqqdr2*denom3)*dr2(i)
c$$$     $      +(denom3*dqqdph*dphi(i))*swth
c$$$     $      +dvlrdt1*dtheta1(i)+dvlrdt2*dtheta2(i)
c$$$     $      +dvlrdrb*drbig(i)
c$$$       dfit(i)=(denom3*(dqqdr1*dr1(i)+dqqdr2*dr2(i)+dqqdph*dphi(i))
c$$$     $      +(denom1*dc6000t1+denom2*dc8000t1+denom3*dqqdt1)*dtheta1(i)
c$$$     $      +(denom1*dc6000t2+denom2*dc8000t2+denom3*dqqdt2)*dtheta2(i)
c$$$     $      +drbig(i)*(-6d0*c6000*denom1*denom1*(rbig**5)
c$$$     $      -8d0*c8000*denom2*denom2a*(rbig**3)
c$$$     $      -5d0*qq*denom3*denom3*(rbig**4)))*swth
c$$$     $      +vlro*(dswthab*(dr1(i)+dr2(i))+dswthc*drbig(i))
c$$$       dfit(i)=dswthab*(dr1(i)+dr2(i))+dswthc*drbig(i)
c$$$       dfit(i)=denom3*(dqqdr1*dr1(i)+dqqdr2*dr2(i)+dqqdph*dphi(i))
c$$$     $      +(denom1*dc6000t1+denom2*dc8000t1+denom3*dqqdt1)*dtheta1(i)
c$$$     $      +(denom1*dc6000t2+denom2*dc8000t2+denom3*dqqdt2)*dtheta2(i)
c$$$     $      +drbig(i)*(-6d0*c6000*denom1*denom1*(rbig**5)
c$$$     $      -8d0*c8000*denom2*denom2a*(rbig**3)
c$$$     $      -5d0*qq*denom3*denom3*(rbig**4))
c$$$       dfit(i)=-8d0*drbig(i)*denom2*denom2a*(rbig**3)
c$$$       dfit(i)=-6d0*drbig(i)*denom1*denom1*(rbig**5)
c$$$       dfit(i)=dc8000t1*dtheta1(i)+dc8000t2*dtheta2(i)
c$$$       dfit(i)=dc6000t1*dtheta1(i)+dc6000t2*dtheta2(i)
c$$$       dfit(i)=(dqqdr1*dr1(i)+dqqdr2*dr2(i)
c$$$     $      +dqqdt1*dtheta1(i)+dqqdt2*dtheta2(i)+dqqdph*dphi(i))
c$$$       dfit(i)=-1.5d0*cst1*dtheta1(i)
c$$$      end do
c
cccc      ddr1=(r1-2.1d0)/(r1+2.1d0)
cccc      ddr2=(r2-2.1d0)/(r2+2.1d0)
         r1re = 1d0/(r1+re)
         r2re = 1d0/(r2+re)
         r1re2 = r1re*r1re
         r2re2 = r2re*r2re
      ddr1=(r1-2.1d0)*r1re
      ddr2=(r2-2.1d0)*r2re
      dr12=ddr1*ddr2
      dr1p2=ddr1+ddr2
c$$$      do i=1,12
c$$$       dfit(i)=beta*7d1*(exp(beta*r3)*dr3(i)
c$$$     $      +exp(beta*r4)*dr4(i)
c$$$     $      +exp(beta*r5)*dr5(i)
c$$$     $      +exp(beta*r6)*dr6(i))
c$$$      end do
c$$$      fit=pairwise
c$$$      if(dr1(12).ne.-132d0)return
      bra=pairwise*swth
      brb=bra*rbig
c
        dpairdr3 = beta70*exp(beta*r3)
        dpairdr4 = beta70*exp(beta*r4)
        dpairdr5 = beta70*exp(beta*r5)
        dpairdr6 = beta70*exp(beta*r6)
c$$$        dbradr1 = pairwise*dswthdr1
c$$$        dbradr2 = pairwise*dswthdr2
        dbradr1 = pairwise*dswthab                                      8d31s09
        if(iprt.ne.0)then
         write(6,*)('pairwise,dswthab '),pairwise,dswthab,dbradr1
        end if
        dbradr2 = pairwise*dswthab                                      8d31s09
        dbradr3 = swth*dpairdr3
        dbradr4 = swth*dpairdr4
        dbradr5 = swth*dpairdr5
        dbradr6 = swth*dpairdr6
c$$$        do i=1,12
c$$$         dfit(i)=dbradr1*dr1(i)+dbradr2*dr2(i)
c$$$     $        +dbradr3*dr3(i)+dbradr4*dr4(i)
c$$$     $        +dbradr5*dr5(i)+dbradr6*dr6(i)
c$$$     $        +pairwise*dswthc*drbig(i)
c$$$        end do
c$$$        fit=bra
c$$$        if(dr1(12).ne.-132d0)return
        dbrbdrb = bra
        dbrbdr1 = rbig*dbradr1
        dbrbdr2 = rbig*dbradr2
        dbrbdr3 = rbig*dbradr3
        dbrbdr4 = rbig*dbradr4
        dbrbdr5 = rbig*dbradr5
        dbrbdr6 = rbig*dbradr6
c$$$        do i=1,12
c$$$         dfit(i)=dbrbdr1*dr1(i)+dbrbdr2*dr2(i)
c$$$     $        +dbrbdr3*dr3(i)+dbrbdr4*dr4(i)
c$$$     $        +dbrbdr5*dr5(i)+dbrbdr6*dr6(i)
c$$$     $        +(pairwise*dswthc*rbig+bra)*drbig(i)
c$$$        end do
c$$$        fit=brb
c$$$        if(dr1(12).ne.-132d0)return
c$$$        write(6,5) vlr,bra,brb,pairwise,dvlrdr1,
c$$$     1    dvlrdr2,dvlrdt1,dvlrdt2,dvlrdrb,
c$$$     2    dpairdr3,dpairdr4,dpairdr5,dpairdr6, dbradr1,dbradr2,dbradr3,
c$$$     3    dbradr4,dbradr5,dbradr6, dbrbdr1,dbrbdr2,dbrbdr3,dbrbdr4,
c$$$     4    dbrbdr5,dbrbdr6
c$$$c23456789012345678901234567890123456789012345678901234567890123456789012
c$$$  5     format(5x,'vlr,bra,brb,pair',4e15.5/12x,
c$$$     1  'dvlr',5e15.5/11x,'dpair',
c$$$     2  4e15.5/12x,'dbra',6e15.5/12x,'dbrb',7e15.5/)
c
cccc      fit=bra*(coef(1)+coef(6)*y2+coef(10)*y3+coef(12)*y4
cccc     $     +coef(15)*y5+coef(18)*y6
cccc     $       +dr1p2*(coef(2)+coef(7)*y2+coef(13)*y4+coef(16)*y5))
cccc     $     +brb*(coef(3)+coef(8)*y2+coef(11)*y3+coef(14)*y4
cccc     $     +coef(17)*y5+coef(19)*y6
cccc     $     +dr1p2*(coef(4)+coef(9)*y2)+dr12*coef(5))
      func1a = (coef(2) +coef(7)*y2 +coef(13)*y4 +coef(16)*y5)
      func2a = (coef(4) +coef(9)*y2)  
      func1= coef(1) +coef(6)*y2 +coef(10)*y3 +coef(12)*y4
     $     +coef(15)*y5 +coef(18)*y6 +dr1p2*func1a
      func2 = coef(3) +coef(8)*y2 +coef(11)*y3 +coef(14)*y4
     $     +coef(17)*y5 +coef(19)*y6 +dr1p2*func2a +dr12*coef(5)
c
      dfunc1dr1 = func1a*re2*r1re2
      dfunc1dr2 = func1a*re2*r2re2
      dfunc1dt1 = -2d0*cst1*(coef(6)+coef(10)*y2b-coef(12)*z2b*cp
c$$$     2   +2d0*coef(15)*y2a-coef(18)*cp -dr1p2*(coef(7)-coef(13)*z2b*cp
     2   +2d0*coef(15)*y2a-coef(18)*cp +dr1p2*(coef(7)-coef(13)*z2b*cp  8d31s09
     3   +2d0*coef(16)*y2a))
      dfunc1dt2 = -2d0*cst2*(coef(6)+coef(10)*y2a-coef(12)*z2a*cp
c$$$     2   +2d0*coef(15)*y2b-coef(18)*cp -dr1p2*(coef(7)-coef(13)*z2a*cp
     2   +2d0*coef(15)*y2b-coef(18)*cp +dr1p2*(coef(7)-coef(13)*z2a*cp  8d31s09
     3   +2d0*coef(16)*y2b))
      dfunc1dph = -2d0*sp*((coef(12)+coef(13)*dr1p2)*z2a*z2b  
     2   +coef(18)*(z2a+z2b))
      dfunc2dr1 = (re2*r1re2)*(func2a+coef(5)*ddr2)
      dfunc2dr2 = (re2*r2re2)*(func2a+coef(5)*ddr1)
      dfunc2dt1 = -2d0*cst1*(coef(8)+coef(11)*y2b-coef(14)*z2b*cp
     2   +2d0*coef(17)*y2a-coef(19)*cp +dr1p2*coef(9))
      dfunc2dt2 = -2d0*cst2*(coef(8)+coef(11)*y2a-coef(14)*z2a*cp
     2   +2d0*coef(17)*y2b-coef(19)*cp +dr1p2*coef(9))
      dfunc2dph = -2d0*sp*(coef(14)*z2a*z2b +coef(19)*(z2a+z2b))
c$$$      do i=1,12
c     ok func2
c$$$       dfit(i)=dfunc2dr1*dr1(i)+dfunc2dr2*dr2(i)
c$$$     $      +dfunc2dt1*dtheta1(i)+dfunc2dt2*dtheta2(i)
c$$$     $      +dfunc2dph*dphi(i)
c     ok y2 (6)
c$$$       dfit(i)=-2d0*(cst1*dtheta1(i)+cst2*dtheta2(i))
c     ok y3 (10)
c$$$       dfit(i)=-2d0*(cst1*y2b*dtheta1(i)+cst2*y2a*dtheta2(i))
c     ok y4 (12)
c$$$       dfit(i)=2d0*(cst1*z2b*cp*dtheta1(i)+cst2*z2a*cp*dtheta2(i))
c$$$     $      -2d0*sp*z2a*z2b*dphi(i)
c     ok y5 (15)
c$$$       dfit(i)=-4d0*(cst1*y2a*dtheta1(i)+cst2*y2b*dtheta2(i))
c     ok y6 (18)
c$$$       dfit(i)=2d0*cp*(cst1*dtheta1(i)+cst2*dtheta2(i))
c$$$     $      -2d0*sp*(z2a+z2b)*dphi(i)
c     ok fun1a
c$$$       dfit(i)=2d0*cst1*dtheta1(i)*(-coef(7)+coef(13)*z2b*cp
c$$$     $      -2d0*coef(16)*y2a)
c$$$     $      +2d0*cst2*dtheta2(i)*(-coef(7)+coef(13)*z2a*cp
c$$$     $      -2d0*coef(16)*y2b)
c$$$     $      -2d0*sp*coef(13)*z2a*z2b*dphi(i)
c     ok func1
c$$$       dfit(i)=dfunc1dr1*dr1(i)+dfunc1dr2*dr2(i)
c$$$     $      +dfunc1dt1*dtheta1(i)+dfunc1dt2*dtheta2(i)
c$$$     $      +dfunc1dph*dphi(i)
c$$$      end do
c$$$      fit=func2
c$$$      if(r1.ne.-132d0)return
c
      fit=bra*func1 + brb*func2 +pairwise+vlr+v12+v34+coef(20)
      dfitr1 = dbradr1*func1+bra*dfunc1dr1+dbrbdr1*func2
     2    +brb*dfunc2dr1 + dvlrdr1+dv12dr1
      if(iprt.gt.0)then
       write(6,*)('dfitr1: '),dbradr1,func1,bra,dfunc1dr1,
     $      dbrbdr1,func2,brb,dfunc2dr1,dvlrdr1,dv12dr1
      end if
      dfitr2 = dbradr2*func1+bra*dfunc1dr2+dbrbdr2*func2
     2    +brb*dfunc2dr2 + dvlrdr2+dv34dr2
      dfitt1 = bra*dfunc1dt1 + brb*dfunc2dt1 + dvlrdt1 
      dfitt2 = bra*dfunc1dt2 + brb*dfunc2dt2 + dvlrdt2 
      dfitph = bra*dfunc1dph + brb*dfunc2dph
     $    +dvlrph                                                       8d31s09
c$$$      dfitrb = dbrbdrb*func2 + dvlrdrb
      dfitrb = dvlrdrb                                                  8d31s09
     $    +pairwise*dswthc*(func1+rbig*func2)+bra*func2                 8d31s09
      dfitr3 = dbradr3*func1 + dbrbdr3*func2 + dpairdr3
      dfitr4 = dbradr4*func1 + dbrbdr4*func2 + dpairdr4
      dfitr5 = dbradr5*func1 + dbrbdr5*func2 + dpairdr5
      dfitr6 = dbradr6*func1 + dbrbdr6*func2 + dpairdr6
      if(iprt.gt.0)then
       write(6,*)('dfitr1: '),dfitr1
       write(6,*)('dfitr2: '),dfitr2,dbradr2,func1,bra,dfunc1dr2,
     $      dbrbdr2,func2,brb,dfunc2dr2,dvlrdr2,dv34dr2

       write(6,*)('dfitt1: '),dfitt1
       write(6,*)('dfitt2: '),dfitt2
       write(6,*)('dfitph: '),dfitph
       write(6,*)('dfitrb: '),dfitrb
       write(6,*)('dfitr3: '),dfitr3
       write(6,*)('dfitr4: '),dfitr4
       write(6,*)('dfitr5: '),dfitr5
       write(6,*)('dfitr6: '),dfitr6
      end if
      do i=1,12
      dfit(i) = dfitr1*dr1(i)+dfitr2*dr2(i)+dfitt1*dtheta1(i)
     2   + dfitt2*dtheta2(i)+dfitph*dphi(i)+dfitrb*drbig(i)
     3   + dfitr3*dr3(i)+dfitr4*dr4(i)+dfitr5*dr5(i)+dfitr6*dr6(i)
      end do
c$$$      do i=1,12
c     ok
c$$$      fit=brb*func2
c$$$       dfit(i)=func2*(dbrbdr1*dr1(i)+dbrbdr2*dr2(i)
c$$$     $        +dbrbdr3*dr3(i)+dbrbdr4*dr4(i)
c$$$     $        +dbrbdr5*dr5(i)+dbrbdr6*dr6(i)
c$$$     $        +(pairwise*dswthc*rbig+bra)*drbig(i))
c$$$     $      +brb*(dfunc2dr1*dr1(i)+dfunc2dr2*dr2(i)
c$$$     $      +dfunc2dt1*dtheta1(i)+dfunc2dt2*dtheta2(i)
c$$$     $      +dfunc2dph*dphi(i))
c$$$      fit=bra*func1
c$$$       dfit(i)=func1*(dbradr1*dr1(i)+dbradr2*dr2(i)
c$$$     $        +dbradr3*dr3(i)+dbradr4*dr4(i)
c$$$     $        +dbradr5*dr5(i)+dbradr6*dr6(i)
c$$$     $        +pairwise*dswthc*drbig(i))
c$$$     $      +bra*(dfunc1dr1*dr1(i)+dfunc1dr2*dr2(i)
c$$$     $      +dfunc1dt1*dtheta1(i)+dfunc1dt2*dtheta2(i)
c$$$     $      +dfunc1dph*dphi(i))
c$$$      end do
c$$$      write(6,6) func1,func2,fit, dfunc1dr1,dfunc1dr2,dfunc1dt1,
c$$$     2   dfunc1dt2,dfunc1dph,dfunc2dr1,dfunc2dr2,dfunc2dt1,dfunc2dt2,
c$$$     3   dfunc1dph 
c$$$  6   format(5x,'func1,func2,fit',3e15.5/5x,'dfunc1',5e15.5/
c$$$     2  5x,'dfunc2',5e15.5/)
c$$$      write(6,7) dfitr1,dfitr2,dfitt1,dfitt2,dfitph,dfitrb,dfitr3,
c$$$     2  dfitr4,dfitr5,dfitr6
c$$$  7   format(5x,'dfit',6e15.5/10x,4e15.5/)
c$$$      write(6,8) (k, dfit(k),k=1,12)
c$$$  8   format(5x,'dfitdx',i2,e15.5)
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diatccsdd(r,ans,dansdr,iwant)
      implicit real*8 (a-h,o-z)
c
c     iwant=1 for pes
c     iwant=2 for qm
c     iwant=3 for dipp
c
      parameter (idp=10,idq=5,idd=5)
      dimension pes(idp),qm(idq),dip(idd)
      data ao,einf,c6,d,co/2.d0,-109.036558873442d0,
     $      -58.8938583622750d0,2.5d0,49d0/
      data npow,a1,re/6,6d0,2.1d0/
      data pes/   146.736533087545d0,   101.160591312956d0,
     $       -251.323812509419d0,  -499.020252425371d0,
     $     -210.074594675137d0,   135.174682961986d0,
     $     214.769359359500d0,   425.317131479472d0,
     $        570.993394092895d0,   266.199319038817d0/
      data a2/2d0/
      data qm/   17.3398942423333d0,   10.2017829193229d0,
     $     -4.71942968995740d0,  -9.61643152091987d0,
     $  -2.12839249419926d0/
      data a3/2d0/
      data dip/  -42.6175278814367d0,   873.328188249918d0,
     $        984.492872053711d0,   518.057930455541d0,
     $        124.185667297002d0/
      if(iwant.eq.1)then
c               ccsd diatomic potential
       r2=r*r
       r4=r2*r2
       r6=r2*r4
       d2=d*d
       d4=d2*d2
       d6=d2*d4
       c8=c6*20d0
       c10=c6*500d0
       vlr=einf+c6/(r6+d6)
     $      +c8/((r4+d4)**2)
     $      +c10/((r2+d2)**5)
       vrep=co*exp(-ao*r)/r
       dvlrdr= -6d0*c6*(r4*r)/(r6+d6)**2
     $      -8d0*c8*(r2*r)/((r4+d4)**3)
     $      -10d0*c10*r/((r2+d2)**6)
       dvrepdr= -vrep*ao - vrep/r
       dr=r-re
       ex=exp(-a1*r)*(r**npow)
       dexdr = -a1*ex + ex*npow/r
       poly=pes(idp)
       do i=idp-1,1,-1
        poly=pes(i)+poly*dr
       end do
c23456789012345678901234567890123456789012345678901234567890123456789012
        dpoly = ((((((((9d0*pes(10)*dr+8d0*pes(9))*dr+7d0*pes(8))*dr
     2   +6d0*pes(7))*dr+5d0*pes(6))*dr+4d0*pes(5))*dr+3d0*pes(4))*dr
     3   +2d0*pes(3))*dr+pes(2))
       ans=vrep+poly*ex+vlr
       dansdr = dvrepdr+dvlrdr+poly*dexdr + ex*dpoly
c23456789012345678901234567890123456789012345678901234567890123456789012
c$$$       write(6,2)r,vlr,vrep,poly,ex,poly*ex
    2  format(1p6e15.7)
      else if(iwant.eq.2)then
c       quadrupole
       ex=exp(-a2*r)
       dr=r-re
       sum=qm(idq)
       do i=idq-1,1,-1
        sum=qm(i)+sum*dr
       end do
       ans=sum*ex*r
       dansdr = sum*ex -a2*ans +r*ex*(((4d0*qm(5)*dr+3d0*qm(4))*dr
     2    +2d0*qm(3))*dr+qm(2))
c$$$       write(6,2)r,ex,sum,ans,ex*r
      else if(iwant.eq.3)then
c         dipole polarizability
       ex=exp(-a3*r)
       dr=r-re
       sum=dip(idd)
       dsum=0d0                                                         2d9s10
       do i=idd-1,2,-1
        dsum=sum+dsum*dr                                                2d9s10
        sum=dip(i)+sum*dr
       end do
       ans=sum*ex*r+dip(1)
c$$$       dansdr = sum*ex -a3*ans +r*ex*(((4d0*dip(5)*dr+3d0*dip(4))*dr
c$$$     2    +2d0*dip(3))*dr+dip(2))
       dansdr=ex*r*dsum+sum*ex-a3*r*ex*sum                              2d9s10
      else
       write(6,*)('unknown iwant: '),iwant
       stop
      end if
      return
      end



      subroutine vdiat(n,v,r,dvdr)
      implicit real*8 (a-h,o-z)
      dimension v(n),r(n),dvdr(n)
      do i=1,n
       ruse=max(1d-4,r(i))
       call vn2d(ruse,v(i),dvdr(i))
       if(dvdr(i).ne.dvdr(i))then
        write(6,*)('got nan in vdiat!!! '),r(i)
        stop
       end if
      end do
      return
      end
      subroutine vn2d(r,vnn,dr)
c     determine two body (N2) potential energy vnn
c     fit to acpf 2s core + LeRoy
      implicit real*8(a-h,o-z)
      dimension csr(12)
      data csr /-3.89327515161896D+02, -6.19410598346194D+02,
     $     -5.51461129947346D+02,
     $      -3.54896837660797D+02, -1.08347448451266D+02,
     $     -6.59348244094835D+01,  1.30457802135760D+01,
     $      7.20557758084161D+01, -1.81062803146583D+01,
     $     -2.84137057950277D+01,  1.40509401686544D+01,
     $     -1.84651681798865D+00/
      data co/  49D0/
      data ao/  2.4D+00/
      data a1/  5.38484562702061D+00/
      data clr/ 23.9717901220746d0/
      data d /3d0/
      data nrep /12/
      data nlr /10/
      data re /2.1d0/
      vdlr=0.d0
      vrep=co*dexp(-ao*r)/r
      drep=-vrep*(ao+(1d0/r))
c$$$      vrep=co*dexp(-ao*r)
c$$$      drep=-vrep*ao
      suma=csr(nrep)
      dsuma=0d0
      dr=r-re
      do n=nrep-1,1,-1
       dsuma=suma+dr*dsuma
       suma=csr(n)+suma*dr
      end do
c$$$      dsuma=suma+r*dsuma
c$$$      suma=suma*r
c$$$       write(6,*)r,suma
c$$$      suma=0.d0
c$$$      xs=1.0d0
c$$$       do n=1,nrep
c$$$       xs=r*xs
c$$$       suma=suma+xs*csr(n)
c$$$       enddo
      ex=exp(-a1*r)*(r**6)
      vsr=suma*ex
c$$$      vrep=vrep+suma*ex
      drep=drep+((-a1+(6d0/r))*suma+dsuma)*ex
 30   r2=r*r
      c6=clr
      c8=clr*20d0
      c10=clr*500d0
      r3=r2*r
      r4=r2*r2
      r5=r4*r
      r6=r4*r2
      d2=d*d
      d4=d2*d2
      d6=d4*d2
      vdlr=vdlr-c6/(r6+d6)
      dr=drep+c6*6d0*r5/((r6+d6)**2)
      vdlr=vdlr-c8/((r4+d4)**2)
      dr=dr+8d0*c8*r3/((r4+d4)**3)
      vdlr=vdlr-c10/((r2+d2)**5)
      dr=dr+10d0*c10*r/((r2+d2)**6)
 50   vnn=vrep+vdlr+vsr
c$$$      write(6,*)vdlr,vsr,vrep,suma
c$$$      vnn=vrep
      return
      end


program ComputeDiat

  use PES_N4_singlet_umn_v3 ,only: pot, pot_Mod

  implicit none 
      
  double precision :: X(4),   Y(4),  Z(4)
  double precision :: E
  double precision :: dEdX(4), dEdY(4), dEdZ(4)
  double precision :: R(6)
  double precision :: dVdR(6)
  double precision ,parameter     ::     Pi  = acos( - 1.d0 )
  
  integer                    :: Unit, Status
  character(200)             :: FileName
  double precision,parameter :: Cconv = 0.52917721092d0
  double precision,parameter :: Econv = 0.159360144d-2
  integer                    :: i, j, iA
  double precision           :: RStart, REnd, hGrid
  double precision           :: RStart1, REnd1, hGrid1
  double precision           :: RStart2, REnd2, hGrid2
  integer                    :: NPoints, NPoints1, NPoints2
  double precision           :: alpha
  character(3)               :: alphaChar
  double precision           :: alphaVec(7)
!  
!  X    = [0.0d0,     0.0d0,   1.d10]
!  Y    = [0.0d0,     0.0d0,   1.d10]    
!  Z    = [0.0d0,     0.0d0,   1.d10]  
!  
!  NPoints = 10000
!  RStart  = 1.5d0
!  REnd    = 10.d0
!  hGrid = (REnd - RStart) / (NPoints-1)
!  X(2)  = RStart
!  
!  FileName = './DiatPot.dat'
!  open( File=trim(adjustl(FileName)), NewUnit=Unit, status='REPLACE', iostat=Status )
!  
!    do i=1,NPoints
!    
!      call pot(X,Y,Z,E,dEdX,dEdY,dEdZ)
!      
!      write(Unit,'(3d20.10)') X(2), E, dEdX(2)
!      
!      X(2) = X(2) + hGrid
!      
!    end do
!    
!  close(Unit)
!  
  
  NPoints1 = 300
  NPoints2 = 300
  RStart1  = 1.5d0
  REnd1    = 10.d0
  RStart2  = 1.5d0
  REnd2    = 10.d0
  hGrid1   = (REnd1 - RStart1) / (NPoints1-1)
  hGrid2   = (REnd2 - RStart2) / (NPoints2-1)

  R        = 1.e3
  
  alphaVec = [60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0]
  
  do iA=1,size(alphaVec,1)
    alpha = alphaVec(iA)
    write(alphaChar,'(I3)') int(alpha)
  
    FileName = './PESFromGrid.csv.' // trim(adjustl(alphaChar))
    open( File=trim(adjustl(FileName)), NewUnit=Unit, status='REPLACE', iostat=Status )
      write(Unit,'(A)') 'R1,R2,R3,V,dV1,dV2,dV3'
      
      R(1) = RStart1
      do i=1,NPoints1
        R(4) = RStart2
        do j=1,NPoints2
        
          R(2) = sqrt( R(1)**2 + R(4)**2 - 2.d0 * R(1) * R(4) * cos(alpha/180.d0*Pi) ) 
          
          call n4fitd(R*0.52917721092, E, dVdR)
          write(Unit,'(e20.10,6(A,e20.10))') R(1), ',', R(2), ',', R(4), ',', E, ',', dVdR(1), ',', dVdR(2), ',', dVdR(4)
          
          R(4) = R(4) + hGrid2
        end do
        R(1) = R(1) + hGrid1
      end do
      
    close(Unit)
    
  end do
  
  
end program ComputeDiat
