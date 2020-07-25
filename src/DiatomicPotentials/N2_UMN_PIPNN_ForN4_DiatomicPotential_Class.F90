! -*-F90-*-
!===============================================================================================================
! 
! Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
! 
! Copyright (C) 2018 Berkan Bolkan and Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
!
! Based on "VVTC" (Vectorized Variable stepsize Trajectory Code) by David Schwenke (NASA Ames Research Center). 
! 
! This program is free software; you can redistribute it and/or modify it under the terms of the 
! Version 2.1 GNU Lesser General Public License as published by the Free Software Foundation. 
! 
! This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU Lesser General Public License for more details. 
! 
! You should have received a copy of the GNU Lesser General Public License along with this library; 
! if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
! 
!---------------------------------------------------------------------------------------------------------------
!===============================================================================================================

Module N2_UMN_PIPNN_ForN4_DiatomicPotential_Class

  use Parameters_Module         ,only:  rkp, Zero, Half, One, Two, Three, Four, Five, Six, Seven, Eight, B_To_Ang, Kcm_To_Hartree, KcmAng_To_HartB, Hartree_To_Kcm
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  Private
  public  ::    N2_UMN_PIPNN_ForN4_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    N2_UMN_PIPNN_ForN4_DiatomicPotential_Type
  contains
    procedure         ::    Initialize        =>    Initialize_N2_UMN_PIPNN_ForN4_DiatomicPotential
    procedure         ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_N2
    procedure         ::    DiatomicPotential =>    DiatomicPotential_N2
  End Type

  character(*)  ,parameter  ::    Name_DiatPot    = 'UMN_PIPNN_ForN4'
  logical       ,parameter  ::    i_Debug_Global = .False.
  
  real(rkp)                 ,parameter :: re   = 1.098_rkp
  real(rkp)                 ,parameter :: de   = 225.213_rkp
  real(rkp) ,dimension(0:6) ,parameter :: as   = [2.7475450369759_rkp, 0.218868498415108_rkp, 0.248885765371433_rkp, -0.229295466336412_rkp, -0.653389048592838_rkp, 1.03611964035396_rkp, 1.71287482791961_rkp]
  real(rkp)                 ,parameter :: VRef = Zero
  contains

!________________________________________________________________________________________________________________________________!
Subroutine Initialize_N2_UMN_PIPNN_ForN4_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(N2_UMN_PIPNN_ForN4_DiatomicPotential_Type)  ,intent(out)  ::    This
  type(Input_Type)                      ,intent(in)   ::    Input
  character(:) ,allocatable             ,intent(in)   ::    SpeciesName
  integer                               ,intent(in)   ::    iMol
  real(rkp)                             ,intent(in)   ::    Mass1
  real(rkp)                             ,intent(in)   ::    Mass2
  logical ,optional                     ,intent(in)   ::    i_Debug
  
  logical                                             ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_N2_UMN_PIPNN_ForN4_DiatomicPotential" )
  !i_Debug_Loc   =     Logger%On()
  
  allocate( This%Name        ,source = trim(Name_DiatPot) )
  allocate( This%SpeciesName ,source = trim(adjustl(SpeciesName)) )
  This%iMol         =    iMol
  This%Initialized  =   .True.

  This%RedMass = Mass1 * Mass2 / (Mass1 + Mass2)
  This%xmui    = One  / This%RedMass            ! Computing the inverse of the target reduced mass [1/a.u.]
  This%xmui2   = Half * This%xmui
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Function DiatomicPotential_N2( This, R ) result( V )

  class(N2_UMN_PIPNN_ForN4_DiatomicPotential_Type)  ,intent(in)  :: This
  real(rkp)                                   ,intent(in)  :: R
  real(rkp)                                                :: V
    
  real(rkp)                                                :: RAng
  real(rkp)                                                :: VDiat
  
  RAng = R * B_To_Ang                                                                                                                   

  call Ev2gm2(RAng, VDiat)
  
  V = VDiat * Kcm_To_Hartree + VRef
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!  

!________________________________________________________________________________________________________________________________!
Elemental Subroutine Compute_Vd_dVd_N2( This, R, V, dV )

  class(N2_UMN_PIPNN_ForN4_DiatomicPotential_Type)  ,intent(in)  :: This
  real(rkp)                                   ,intent(in)  :: R
  real(rkp)                                   ,intent(out) :: V
  real(rkp)                                   ,intent(out) :: dV
  
  real(rkp)                                    :: RAng
  real(rkp)                                    :: VDiat
  real(rkp)                                    :: dVDiat
  
  RAng = R * B_To_Ang
  
  call Ev2gm2_Grad(RAng, VDiat, dVDiat)
  
  V  =  VDiat * Kcm_To_Hartree + VRef                                                                       
  dV = dVDiat * Kcm_To_Hartree * B_To_Ang
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine Ev2gm2( R, V )

!**********************************************************************
!
! This subroutine evaulates the 2-body potential energy and gradient
! for given r with two kinds of potentials:
!
! A polinomial based on MEG variables and a Lennard-Jones term
! V(r) = 0 for r -> inifity
! u=exp(-(r-red)/alpha-(r-red)^2/beta))
! V(r)=-De*(cs1*u^1+cs2*u^2+cs3*u^3+cs4*u^4+cs5*u^5+cs6*u^6)
!
!**********************************************************************

  use Parameters_Module     ,only:  Zero, Two

    real(rkp)               ,intent(in)  :: R
    real(rkp)               ,intent(out) :: V
 
    real(rkp)                            :: y, y2, y3, y4, y5, y6
    real(rkp)                            :: r3, r4 
    real(rkp)                            :: re4
    real(rkp)                            :: fy
    real(rkp)                            :: u
    real(rkp)                            :: dfdy
    real(rkp)                            :: rr
    real(rkp)                            :: dydr
    real(rkp)                            :: dfdr
    integer                              :: k
    ! Dispersion variables
    real(rkp)                            :: dist
    real(rkp)                            :: disp
    real(rkp)                            :: dispdr

    v   = Zero
    r3  = r**3
    r4  = r**4
    re4 = re**4
    y   = (r4 - re4)/(r4 + re4)
    y2  = y*y
    y3  = y*y2
    y4  = y2*y2
    y5  = y2*y3
    y6  = y3*y3

    fy = as(0) + as(1)*y + as(2)*y2 + as(3)*y3 + as(4)*y4 + as(5)*y5 + as(6)*y6

    u  = exp(-fy*(r-re))
    v  = de*(One-u)*(One-u)-de

    ! Add D3 dispersion correction
    dist=r
    call d3disp(dist,disp,dispdr,0)
    v=v+disp
  
End Subroutine Ev2gm2
!--------------------------------------------------------------------------------------------------------------------------------! 


!________________________________________________________________________________________________________________________________!
Pure Subroutine Ev2gm2_Grad( R, V, dV ) 
  
  real(rkp)               ,intent(in)  :: R
  real(rkp)               ,intent(out) :: v
  real(rkp)               ,intent(out) :: dV
  
  real(rkp)                            :: y, y2, y3, y4, y5, y6
  real(rkp)                            :: r3, r4 
  real(rkp)                            :: re4
  real(rkp)                            :: fy
  real(rkp)                            :: u
  real(rkp)                            :: dfdy
  real(rkp)                            :: rr
  real(rkp)                            :: dydr
  real(rkp)                            :: dfdr
  integer                              :: k
  ! Dispersion variables
  real(rkp)                            :: dist
  real(rkp)                            :: disp
  real(rkp)                            :: dispdr

  v   = Zero
  r3  = r**3
  r4  = r**4
  re4 = re**4
  y   = (r4 - re4)/(r4 + re4)
  y2  = y*y
  y3  = y*y2
  y4  = y2*y2
  y5  = y2*y3
  y6  = y3*y3

  fy = as(0) + as(1)*y + as(2)*y2 + as(3)*y3 + as(4)*y4 + as(5)*y5 + as(6)*y6

  u  = exp(-fy*(r-re))
  v  = de*(One-u)*(One-u)-de


  ! Compute the gradient if needed
  dV   = Zero
  dfdy = as(1) + Two*as(2)*y + Three*as(3)*y2 + Four*as(4)*y3 + Five*as(5)*y4 + Six*as(6)*y5
  rr   = r4 + re4
  dydr = Eight*r3*re4/(rr*rr)
  dfdr = dfdy*dydr
  dV   = Two*de*(One-u)*u*(dfdr*(r-re)+fy)

  ! Add D3 dispersion correction
  ! Add analytical gradient of D3 dispersion correction
  dist = r 
  call d3disp_Grad(dist,disp,dispdr,0)
  v    = v  + disp
  dV   = dV + dispdr
   
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine d3disp(distAng, disp, dispdr, igrad)
!**********************************************************************
! Dispersion correction based on Grimme's D3(BJ) calculation for
! diatomic pairs
!
! Several subroutines of DFTD3 V3.1 Rev 1 by Grimme were merged into
! subroutine edisp and they have been heavily modified to calculate only
! that dispersion energy correction that is needed.
!
! S. Grimme, J. Antony, S. Ehrlich and H. Krieg
! J. Chem. Phys, 132 (2010), 154104
! and
! S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011),
! 1456-1465
!
! The C6 values are fixed.
!
!**********************************************************************
  
  real(rkp)               ,intent(in)  :: distAng
  real(rkp)               ,intent(out) :: disp
  real(rkp)               ,intent(out) :: dispdr
  integer                 ,intent(in)  :: igrad

  real(rkp)                            :: dist
  real(rkp)                            :: s6, s8, rs6, rs8
  real(rkp)                            :: e6, e8
  real(rkp)                            :: c6
  real(rkp)                            :: e6dr
  real(rkp)                            :: e8dr
  integer                              :: iz
  integer                              :: i, j
  real(rkp) ,dimension(94)             :: r2r4
 

  ! Generalized parameters for BJ damping from P. Verma, B. Wang,
  ! L. E. Fernandez, and D. G. Truhlar, J. Phys. Chem. A 121, 2855 (2017)
  s6  = One
  s8  = Two
  rs6 = 0.5299d0
  rs8 = 2.20d0

  dist = distAng / B_To_Ang

  ! iz for N2 system
  iz = 7
  ! C6 for N2 system
  c6 = 19.7d0

  ! Calculate dispersion correction
  call edisp(94,5,2,dist,iz,rs6,rs8,e6,e8,e6dr,e8dr,c6,0)

  disp = Zero
  disp = disp + (-s6*e6-s8*e8)*Hartree_To_Kcm

End Subroutine d3disp
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine d3disp_Grad(distAng, disp, dispdr, igrad)
!**********************************************************************
! Dispersion correction based on Grimme's D3(BJ) calculation for
! diatomic pairs
!
! Several subroutines of DFTD3 V3.1 Rev 1 by Grimme were merged into
! subroutine edisp and they have been heavily modified to calculate only
! that dispersion energy correction that is needed.
!
! S. Grimme, J. Antony, S. Ehrlich and H. Krieg
! J. Chem. Phys, 132 (2010), 154104
! and
! S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011),
! 1456-1465
!
! The C6 values are fixed.
!
!**********************************************************************
  
  real(rkp)               ,intent(in)  :: distAng
  real(rkp)               ,intent(out) :: disp
  real(rkp)               ,intent(out) :: dispdr
  integer                 ,intent(in)  :: igrad

  real(rkp)                            :: dist
  real(rkp)                            :: s6, s8, rs6, rs8
  real(rkp)                            :: e6, e8
  real(rkp)                            :: c6
  real(rkp)                            :: e6dr
  real(rkp)                            :: e8dr
  integer                              :: iz
  integer                              :: i, j
  real(rkp) ,dimension(94)             :: r2r4
 

  ! Generalized parameters for BJ damping from P. Verma, B. Wang,
  ! L. E. Fernandez, and D. G. Truhlar, J. Phys. Chem. A 121, 2855 (2017)
  s6  = One
  s8  = Two
  rs6 = 0.5299d0
  rs8 = 2.20d0

  dist = distAng / B_To_Ang

  ! iz for N2 system
  iz = 7
  ! C6 for N2 system
  c6 = 19.7d0

  ! Calculate dispersion correction
  call edisp_Grad(94,5,2,dist,iz,rs6,rs8,e6,e8,e6dr,e8dr,c6,0)

  disp   = (-s6*e6-s8*e8)     * Hartree_To_Kcm

  dispdr = (-s6*e6dr-s8*e8dr) * Hartree_To_Kcm/B_To_Ang

End Subroutine d3disp_Grad
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine edisp(max_elem, maxc, n, dist, iz, rs6, rs8, e6, e8, e6dr, e8dr, c6a, igrad)
!**********************************************************************
! compute energy
!**********************************************************************

  integer                                             ,intent(in)  :: max_elem
  integer                                             ,intent(in)  :: maxc
  integer                                             ,intent(in)  :: n
  real(rkp)                                           ,intent(in)  :: dist
  integer                                             ,intent(in)  :: iz
  real(rkp)                                           ,intent(in)  :: rs6
  real(rkp)                                           ,intent(in)  :: rs8
  real(rkp)                                           ,intent(out) :: e6
  real(rkp)                                           ,intent(out) :: e8
  real(rkp)                                           ,intent(out) :: e6dr
  real(rkp)                                           ,intent(out) :: e8dr
  real(rkp)                                           ,intent(in)  :: c6a
  integer                                             ,intent(in)  :: igrad

  integer                                                          :: iat,jat
  real(rkp)                                                        :: r,tmp,c6,c8,a1,a2
  real(rkp)                                                        :: damp6,damp8
  integer                                                          :: step

  !  r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
  !  the large number of digits is just to keep the results consistent
  !  with older versions. They should not imply any higher accuracy than
  !  the old values
  real(rkp) ,dimension(94)                              ,parameter :: r2r4 = [2.00734898d0,  1.56637132d0,  5.01986934d0,  3.85379032d0, &
                                                                              3.64446594d0,  3.10492822d0,  2.71175247d0,  2.59361680d0, &
                                                                              2.38825250d0,  2.21522516d0,  6.58585536d0,  5.46295967d0, &
                                                                              5.65216669d0,  4.88284902d0,  4.29727576d0,  4.04108902d0, &
                                                                              3.72932356d0,  3.44677275d0,  7.97762753d0,  7.07623947d0, &
                                                                              6.60844053d0,  6.28791364d0,  6.07728703d0,  5.54643096d0, &
                                                                              5.80491167d0,  5.58415602d0,  5.41374528d0,  5.28497229d0, &
                                                                              5.22592821d0,  5.09817141d0,  6.12149689d0,  5.54083734d0, &
                                                                              5.06696878d0,  4.87005108d0,  4.59089647d0,  4.31176304d0, &
                                                                              9.55461698d0,  8.67396077d0,  7.97210197d0,  7.43439917d0, &
                                                                              6.58711862d0,  6.19536215d0,  6.01517290d0,  5.81623410d0, &
                                                                              5.65710424d0,  5.52640661d0,  5.44263305d0,  5.58285373d0, &
                                                                              7.02081898d0,  6.46815523d0,  5.98089120d0,  5.81686657d0, &
                                                                              5.53321815d0,  5.25477007d0, 11.02204549d0,  0.15679528d0, &
                                                                              9.35167836d0,  9.06926079d0,  8.97241155d0,  8.90092807d0, &
                                                                              8.85984840d0,  8.81736827d0,  8.79317710d0,  7.89969626d0, &
                                                                              8.80588454d0,  8.42439218d0,  8.54289262d0,  8.47583370d0, &
                                                                              8.45090888d0,  8.47339339d0,  7.83525634d0,  8.20702843d0, &
                                                                              7.70559063d0,  7.32755997d0,  7.03887381d0,  6.68978720d0, &
                                                                              6.05450052d0,  5.88752022d0,  5.70661499d0,  5.78450695d0, &
                                                                              7.79780729d0,  7.26443867d0,  6.78151984d0,  6.67883169d0, &
                                                                              6.39024318d0,  6.09527958d0, 11.79156076d0, 11.10997644d0, &
                                                                              9.51377795d0,  8.67197068d0,  8.77140725d0,  8.65402716d0, &
                                                                              8.53923501d0,  8.85024712d0 ]

  ! these new data are scaled with k2=4./3.  and converted a_0 via
  ! B_To_Ang=0.52917726d0
  real(rkp) ,dimension(94)                              ,parameter :: rcov = [0.80628308d0, 1.15903197d0, 3.02356173d0, 2.36845659d0, &
                                                                              1.94011865d0, 1.88972601d0, 1.78894056d0, 1.58736983d0, &
                                                                              1.61256616d0, 1.68815527d0, 3.52748848d0, 3.14954334d0, &
                                                                              2.84718717d0, 2.62041997d0, 2.77159820d0, 2.57002732d0, &
                                                                              2.49443835d0, 2.41884923d0, 4.43455700d0, 3.88023730d0, &
                                                                              3.35111422d0, 3.07395437d0, 3.04875805d0, 2.77159820d0, &
                                                                              2.69600923d0, 2.62041997d0, 2.51963467d0, 2.49443835d0, &
                                                                              2.54483100d0, 2.74640188d0, 2.82199085d0, 2.74640188d0, &
                                                                              2.89757982d0, 2.77159820d0, 2.87238349d0, 2.94797246d0, &
                                                                              4.76210950d0, 4.20778980d0, 3.70386304d0, 3.50229216d0, &
                                                                              3.32591790d0, 3.12434702d0, 2.89757982d0, 2.84718717d0, &
                                                                              2.84718717d0, 2.72120556d0, 2.89757982d0, 3.09915070d0, &
                                                                              3.22513231d0, 3.17473967d0, 3.17473967d0, 3.09915070d0, &
                                                                              3.32591790d0, 3.30072128d0, 5.26603625d0, 4.43455700d0, &
                                                                              4.08180818d0, 3.70386304d0, 3.98102289d0, 3.95582657d0, &
                                                                              3.93062995d0, 3.90543362d0, 3.80464833d0, 3.82984466d0, &
                                                                              3.80464833d0, 3.77945201d0, 3.75425569d0, 3.75425569d0, &
                                                                              3.72905937d0, 3.85504098d0, 3.67866672d0, 3.45189952d0, &
                                                                              3.30072128d0, 3.09915070d0, 2.97316878d0, 2.92277614d0, &
                                                                              2.79679452d0, 2.82199085d0, 2.84718717d0, 3.32591790d0, &
                                                                              3.27552496d0, 3.27552496d0, 3.42670319d0, 3.30072128d0, &
                                                                              3.47709584d0, 3.57788113d0, 5.06446567d0, 4.56053862d0, &
                                                                              4.20778980d0, 3.98102289d0, 3.82984466d0, 3.85504098d0, &
                                                                              3.88023730d0, 3.90543362d0 ]

  e6   = Zero
  e8   = Zero

  e6dr = Zero
  e8dr = Zero

  a1   = rs6
  a2   = rs8

  ! DFT-D3
  step=0
  do iat=1,n-1
    do jat=iat+1,n
      step=step+1
      
      r  = dist
      c6 = c6a
      ! r2r4 stored in main as sqrt
      c8 = Three*c6*r2r4(iz)*r2r4(iz)

      ! energy for BJ damping
      tmp = sqrt(c8/c6)
      e6  = c6/(r**6+(a1*tmp+a2)**6)
      e8  = c8/(r**8+(a1*tmp+a2)**8)
      
      ! calculate gradients
      if (igrad .eq. 1) then
        ! grad for BJ damping
        e6dr = c6*(-6*r**5)/(r**6+(a1*tmp+a2)**6)**2
        e8dr = c8*(-8*r**7)/(r**8+(a1*tmp+a2)**8)**2
      endif

    enddo
  enddo

End Subroutine edisp
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine edisp_Grad(max_elem, maxc, n, dist, iz, rs6, rs8, e6, e8, e6dr, e8dr, c6a, igrad)
!**********************************************************************
! compute energy
!**********************************************************************

  integer                                             ,intent(in)  :: max_elem
  integer                                             ,intent(in)  :: maxc
  integer                                             ,intent(in)  :: n
  real(rkp)                                           ,intent(in)  :: dist
  integer                                             ,intent(in)  :: iz
  real(rkp)                                           ,intent(in)  :: rs6
  real(rkp)                                           ,intent(in)  :: rs8
  real(rkp)                                           ,intent(out) :: e6
  real(rkp)                                           ,intent(out) :: e8
  real(rkp)                                           ,intent(out) :: e6dr
  real(rkp)                                           ,intent(out) :: e8dr
  real(rkp)                                           ,intent(in)  :: c6a
  integer                                             ,intent(in)  :: igrad

  integer                                                          :: iat,jat
  real(rkp)                                                        :: r,tmp,c6,c8,a1,a2
  real(rkp)                                                        :: damp6,damp8
  integer                                                          :: step

  !  r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
  !  the large number of digits is just to keep the results consistent
  !  with older versions. They should not imply any higher accuracy than
  !  the old values
  real(rkp) ,dimension(94)                              ,parameter :: r2r4 = [2.00734898d0,  1.56637132d0,  5.01986934d0,  3.85379032d0, &
                                                                              3.64446594d0,  3.10492822d0,  2.71175247d0,  2.59361680d0, &
                                                                              2.38825250d0,  2.21522516d0,  6.58585536d0,  5.46295967d0, &
                                                                              5.65216669d0,  4.88284902d0,  4.29727576d0,  4.04108902d0, &
                                                                              3.72932356d0,  3.44677275d0,  7.97762753d0,  7.07623947d0, &
                                                                              6.60844053d0,  6.28791364d0,  6.07728703d0,  5.54643096d0, &
                                                                              5.80491167d0,  5.58415602d0,  5.41374528d0,  5.28497229d0, &
                                                                              5.22592821d0,  5.09817141d0,  6.12149689d0,  5.54083734d0, &
                                                                              5.06696878d0,  4.87005108d0,  4.59089647d0,  4.31176304d0, &
                                                                              9.55461698d0,  8.67396077d0,  7.97210197d0,  7.43439917d0, &
                                                                              6.58711862d0,  6.19536215d0,  6.01517290d0,  5.81623410d0, &
                                                                              5.65710424d0,  5.52640661d0,  5.44263305d0,  5.58285373d0, &
                                                                              7.02081898d0,  6.46815523d0,  5.98089120d0,  5.81686657d0, &
                                                                              5.53321815d0,  5.25477007d0, 11.02204549d0,  0.15679528d0, &
                                                                              9.35167836d0,  9.06926079d0,  8.97241155d0,  8.90092807d0, &
                                                                              8.85984840d0,  8.81736827d0,  8.79317710d0,  7.89969626d0, &
                                                                              8.80588454d0,  8.42439218d0,  8.54289262d0,  8.47583370d0, &
                                                                              8.45090888d0,  8.47339339d0,  7.83525634d0,  8.20702843d0, &
                                                                              7.70559063d0,  7.32755997d0,  7.03887381d0,  6.68978720d0, &
                                                                              6.05450052d0,  5.88752022d0,  5.70661499d0,  5.78450695d0, &
                                                                              7.79780729d0,  7.26443867d0,  6.78151984d0,  6.67883169d0, &
                                                                              6.39024318d0,  6.09527958d0, 11.79156076d0, 11.10997644d0, &
                                                                              9.51377795d0,  8.67197068d0,  8.77140725d0,  8.65402716d0, &
                                                                              8.53923501d0,  8.85024712d0 ]

  ! these new data are scaled with k2=4./3.  and converted a_0 via
  ! B_To_Ang=0.52917726d0
  real(rkp) ,dimension(94)                              ,parameter :: rcov = [0.80628308d0, 1.15903197d0, 3.02356173d0, 2.36845659d0, &
                                                                              1.94011865d0, 1.88972601d0, 1.78894056d0, 1.58736983d0, &
                                                                              1.61256616d0, 1.68815527d0, 3.52748848d0, 3.14954334d0, &
                                                                              2.84718717d0, 2.62041997d0, 2.77159820d0, 2.57002732d0, &
                                                                              2.49443835d0, 2.41884923d0, 4.43455700d0, 3.88023730d0, &
                                                                              3.35111422d0, 3.07395437d0, 3.04875805d0, 2.77159820d0, &
                                                                              2.69600923d0, 2.62041997d0, 2.51963467d0, 2.49443835d0, &
                                                                              2.54483100d0, 2.74640188d0, 2.82199085d0, 2.74640188d0, &
                                                                              2.89757982d0, 2.77159820d0, 2.87238349d0, 2.94797246d0, &
                                                                              4.76210950d0, 4.20778980d0, 3.70386304d0, 3.50229216d0, &
                                                                              3.32591790d0, 3.12434702d0, 2.89757982d0, 2.84718717d0, &
                                                                              2.84718717d0, 2.72120556d0, 2.89757982d0, 3.09915070d0, &
                                                                              3.22513231d0, 3.17473967d0, 3.17473967d0, 3.09915070d0, &
                                                                              3.32591790d0, 3.30072128d0, 5.26603625d0, 4.43455700d0, &
                                                                              4.08180818d0, 3.70386304d0, 3.98102289d0, 3.95582657d0, &
                                                                              3.93062995d0, 3.90543362d0, 3.80464833d0, 3.82984466d0, &
                                                                              3.80464833d0, 3.77945201d0, 3.75425569d0, 3.75425569d0, &
                                                                              3.72905937d0, 3.85504098d0, 3.67866672d0, 3.45189952d0, &
                                                                              3.30072128d0, 3.09915070d0, 2.97316878d0, 2.92277614d0, &
                                                                              2.79679452d0, 2.82199085d0, 2.84718717d0, 3.32591790d0, &
                                                                              3.27552496d0, 3.27552496d0, 3.42670319d0, 3.30072128d0, &
                                                                              3.47709584d0, 3.57788113d0, 5.06446567d0, 4.56053862d0, &
                                                                              4.20778980d0, 3.98102289d0, 3.82984466d0, 3.85504098d0, &
                                                                              3.88023730d0, 3.90543362d0 ]

  e6   = Zero
  e8   = Zero

  e6dr = Zero
  e8dr = Zero

  a1   = rs6
  a2   = rs8

  ! DFT-D3
  step=0
  do iat=1,n-1
    do jat=iat+1,n
      step=step+1
      
      r  = dist
      c6 = c6a
      ! r2r4 stored in main as sqrt
      c8 = Three*c6*r2r4(iz)*r2r4(iz)

      ! energy for BJ damping
      tmp = sqrt(c8/c6)
      e6  = c6/(r**6+(a1*tmp+a2)**6)
      e8  = c8/(r**8+(a1*tmp+a2)**8)
      
      ! calculate gradients
      ! grad for BJ damping
      e6dr = c6*(-6*r**5)/(r**6+(a1*tmp+a2)**6)**2
      e8dr = c8*(-8*r**7)/(r**8+(a1*tmp+a2)**8)**2

    enddo
  enddo

End Subroutine edisp_Grad
!--------------------------------------------------------------------------------------------------------------------------------!


End Module