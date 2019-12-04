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

Module O2_UMN_DiatomicPotential_Class

! The LeRoy model is used.

  use Parameters_Module         ,only:  rkp, Zero, Half, One, Two, B_To_Ang, Kcm_To_Hartree, KcmAng_To_HartB
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  private
  public  ::    O2_UMN_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    O2_UMN_DiatomicPotential_Type
  contains
    procedure         ::    Initialize        =>    Initialize_O2_UMN_DiatomicPotential
    procedure         ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_O2_UMN
    procedure         ::    DiatomicPotential =>    DiatomicPotential_O2_UMN
  End Type

  character(*)  ,parameter  ::    Name_DiaPot    = 'UMN'
  logical       ,parameter  ::    i_Debug_Global = .False.
 
  contains

!________________________________________________________________________________________________________________________________!
Subroutine Initialize_O2_UMN_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(O2_UMN_DiatomicPotential_Type)  ,intent(out)  ::    This
  type(Input_Type)                      ,intent(in)   ::    Input
  character(:) ,allocatable             ,intent(in)   ::    SpeciesName
  integer                               ,intent(in)   ::    iMol
  real(rkp)                             ,intent(in)   ::    Mass1
  real(rkp)                             ,intent(in)   ::    Mass2
  logical ,optional                     ,intent(in)   ::    i_Debug
  
  logical                                             ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_O2_UMN_DiatomicPotential" )
  !i_Debug_Loc   =     Logger%On()
  
  allocate( This%Name        ,source = trim(Name_DiaPot) )
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
Elemental Function DiatomicPotential_O2_UMN( This, R ) result( V )

  class(O2_UMN_DiatomicPotential_Type)  ,intent(in)  :: This
  real(rkp)                             ,intent(in)  :: R
  real(rkp)                                          :: V
    
  real(rkp)                                    :: RAng
  real(rkp)                                    :: VDiat

  real(rkp)                      ,parameter    :: VRef     = 0.1915103559 ! Original was 240.486_rkp Kcal/mol - 0.19172848_rkp Hartree
  
  !Reference energy of infinitely separated O2 + O in hartree (taken from DSEC corrected calculations) = 0.19172848 Hartree = 120.311 Kcal/mol
  !The enthalpy of formation of ozone is 34.09895 kcal/mol = 0.05434012 Hartree
  ! TOT = 0.2460686 Hartree = 154.4104107 kcal/mol
  
  RAng = R * B_To_Ang                                                                                                           
        
  call Ev2gm2(RAng, VDiat)
  
  V = VDiat                                                   ! Simone Added 3.0                                             
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!  


!________________________________________________________________________________________________________________________________!
Elemental Subroutine Compute_Vd_dVd_O2_UMN( This, R, V, dV )

  class(O2_UMN_DiatomicPotential_Type)  ,intent(in)  :: This
  real(rkp)                             ,intent(in)  :: R
  real(rkp)                             ,intent(out) :: V
  real(rkp)                             ,intent(out) :: dV
  
  real(rkp)                                    :: RAng
  real(rkp)                                    :: VDiat
  real(rkp)                                    :: dVDiat

  real(rkp)                      ,parameter    :: VRef     = 0.1915103559 ! Original was 240.486_rkp Kcal/mol - 0.19172848_rkp Eh Hartree                                                       
  
  RAng = R * B_To_Ang                                                                                                             !this is a tempprary way of computing the correct values

  call Ev2gm2_Grad(RAng, VDiat, dVDiat)
  
  V  =  VDiat                                                  ! Simone Added 3.0                                                                       
  dV = dVDiat * B_To_Ang
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine Ev2gm2( RAng, V) 
! Compute the diatomic potential of groupd-state triplet O2
!
! References: J. Chem. Phys. 132, 074307 (2010)
!
! Input:  R      interatomic distance in Angstrom
! Output: V      potential in kcal/mol
!         grad   gradient (kcal/mol)/Angstrom

! Parameters of analytical even-tempered Gaussian expansions for the ground state potential energy curve of O2 CBS+SR+SO+CV. Units: alpha in
! Anstromgs**-2, beta=dimensionless, a_k in milihartree.

! Original parameters
!      alpha = 0.785_rkp
!      beta = 1.307_rkp
!      a(0) = -2388.5641690_rkp
!      a(1) = 18086.977116_rkp
!      a(2) = -71760.197585_rkp
!      a(3) = 154738.09175_rkp
!      a(4) = -215074.85646_rkp
!      a(5) = 214799.54567_rkp
!      a(6) = -148395.42850_rkp
!      a(7) = 73310.781453_rkp

  use Parameters_Module     ,only:  Zero, Two

  real(rkp)                                 ,intent(in)     ::    RAng                                                            ! Distances of atom-atom pairs [Angstrom]
  real(rkp)                                 ,intent(out)    ::    V   
  
  real(rkp)                                                 ::    VDisp
  integer                                                   ::    k
  
  real(rkp) ,dimension(0:7)                 ,parameter      ::    a     = [-1.488979427684798e3_rkp, 1.881435846488955e4_rkp,    &
                                                                           -1.053475425838226e5_rkp, 2.755135591229064e5_rkp,    &
                                                                           -4.277588997761775e5_rkp, 4.404104009614092e5_rkp,    & 
                                                                           -2.946204062950765e5_rkp, 1.176861219078620e5_rkp ]        
  real(rkp)                                 ,parameter      ::    alpha = 9.439784362354936e-1_rkp   
  real(rkp)                                 ,parameter      ::    beta  = 1.262242998506810_rkp     

  V = Zero
  do k = 0,7
    V = V + a(k)*exp(-alpha * beta**k * RAng**2)
  end do
  ! From milihartree to hartree
  V = V*1.e-3_rkp
  
  ! D3 dispersion correction  
  call d3disp(RAng, VDisp)
  V = V + VDisp                          
 
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------! 


!________________________________________________________________________________________________________________________________!
Pure Subroutine Ev2gm2_Grad( RAng, V, dV) 
! Compute the diatomic potential of groupd-state triplet O2
!
! References: J. Chem. Phys. 132, 074307 (2010)
!
! Input:  R      interatomic distance in Angstrom
! Output: V      potential in kcal/mol
!         grad   gradient (kcal/mol)/Angstrom

! Parameters of analytical even-tempered Gaussian expansions for the ground state potential energy curve of O2 CBS+SR+SO+CV. Units: alpha in
! Anstromgs**-2, beta=dimensionless, a_k in milihartree.

! Original parameters
!      alpha = 0.785_rkp
!      beta = 1.307_rkp
!      a(0) = -2388.5641690_rkp
!      a(1) = 18086.977116_rkp
!      a(2) = -71760.197585_rkp
!      a(3) = 154738.09175_rkp
!      a(4) = -215074.85646_rkp
!      a(5) = 214799.54567_rkp
!      a(6) = -148395.42850_rkp
!      a(7) = 73310.781453_rkp

  use Parameters_Module     ,only:  Zero, Two

  real(rkp)                                 ,intent(in)     ::    RAng                                                            ! Distances of atom-atom pairs [Angstrom]
  real(rkp)                                 ,intent(out)    ::    V   
  real(rkp)                                 ,intent(out)    ::    dV 
  
  real(rkp)                                                 ::    VDisp
  real(rkp)                                                 ::    dVDisp
  integer                                                   ::    k
  
  real(rkp) ,dimension(0:7)                 ,parameter      ::    a     = [-1.488979427684798e3_rkp, 1.881435846488955e4_rkp,    &
                                                                           -1.053475425838226e5_rkp, 2.755135591229064e5_rkp,    &
                                                                           -4.277588997761775e5_rkp, 4.404104009614092e5_rkp,    & 
                                                                           -2.946204062950765e5_rkp, 1.176861219078620e5_rkp ]        
  real(rkp)                                 ,parameter      ::    alpha = 9.439784362354936e-1_rkp   
  real(rkp)                                 ,parameter      ::    beta  = 1.262242998506810_rkp     

  V = Zero
  do k = 0,7
    V = V + a(k)*exp(-alpha * beta**k * RAng**2)
  end do
  ! From milihartree to hartree
  V = V*1.e-3_rkp
    
  dV = Zero
  do k = 0,7
    dV = dV - Two * a(k) * alpha * beta**k * RAng * exp(-alpha * beta**k * RAng**2)
  end do
  ! Convert from milihartree/A to hartree/A
  dV = dV*1.e-3_rkp

  call d3disp_Grad( RAng, VDisp, dVDisp)
  V  = V  + VDisp    
  dV = dV + dVDisp                              
      
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!     


!________________________________________________________________________________________________________________________________!
Pure Subroutine d3disp(RAng, VDisp)
! Dispersion correction based on Grimme's D3(BJ) calculation for diatomic pairs
!
! Several Subroutines of DFTD3 V3.1 Rev 1 by Grimme were merged into subroutine edisp and they have been heavily modified to calculate only
! that dispersion energy correction that is needed.
!
! S. Grimme, J. Antony, S. Ehrlich and H. Krieg, J. Chem. Phys, 132 (2010), 154104
! and 
! S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011), 1456-1465
!
! The C6 values are fixed.
  use Parameters_Module     ,only:  Zero, Two, One, B_To_Ang, Kcm_To_Hartree, KcmAng_To_HartB 
  
  real(rkp)                              ,intent(in)  ::    RAng
  real(rkp)                              ,intent(out) ::    VDisp    

  real(rkp)                                           ::    e6, e8
  real(rkp)                                           ::    RBohr

  ! Generalized parameters for BJ damping from P. Verma, B. Wang, L. E. Fernandez, and D. G. Truhlar, J. Phys. Chem. A 121, 2855 (2017).

  RBohr = RAng / B_To_Ang
  
  ! Calculate dispersion correction
  call edisp( RBohr, e6, e8 )
  
  VDisp = (-e6 -Two*e8) 

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine d3disp_Grad(RAng, VDisp, dVDisp)
! Dispersion correction based on Grimme's D3(BJ) calculation for diatomic pairs
!
! Several Subroutines of DFTD3 V3.1 Rev 1 by Grimme were merged into subroutine edisp and they have been heavily modified to calculate only
! that dispersion energy correction that is needed.
!
! S. Grimme, J. Antony, S. Ehrlich and H. Krieg, J. Chem. Phys, 132 (2010), 154104
! and 
! S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011), 1456-1465
!
! The C6 values are fixed.
  use Parameters_Module     ,only:  Zero, Two, One, B_To_Ang, Kcm_To_Hartree, KcmAng_To_HartB 
  
  real(rkp)                              ,intent(in)        ::    RAng
  real(rkp)                              ,intent(out)       ::    VDisp    
  real(rkp)                              ,intent(out)       ::    dVDisp  

  real(rkp)                                                 ::    e6, e8, e6dr, e8dr
  real(rkp)                                                 ::    RBohr

  ! Generalized parameters for BJ damping from P. Verma, B. Wang, L. E. Fernandez, and D. G. Truhlar, J. Phys. Chem. A 121, 2855 (2017).

  RBohr = RAng / B_To_Ang
  
  ! Calculate dispersion correction
  call edisp_Grad( RBohr, e6, e8, e6dr, e8dr )

  VDisp  = (-e6   -Two*e8)   
  dVDisp = (-e6dr -Two*e8dr) / B_To_Ang

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!



!________________________________________________________________________________________________________________________________!
Pure Subroutine edisp( RBohr, e6, e8 )
! compute energy
      
  use Parameters_Module     ,only:  Zero, One, Two, Three

  real(rkp)                               ,intent(in)       ::    RBohr  
  real(rkp)                               ,intent(out)      ::    e6
  real(rkp)                               ,intent(out)      ::    e8
  
  real(rkp)                                                 ::    tmp, c8Step
  
  real(rkp)                                   ,parameter    ::    r2r4Scalar = 2.59361680_rkp                      
  real(rkp)                                   ,parameter    ::    rs6        = 0.5299_rkp
  real(rkp)                                   ,parameter    ::    rs8        = 2.20_rkp
  real(rkp)                                   ,parameter    ::    c6         = 12.8_rkp

  ! DFT-D3
  c8Step = Three * c6 * r2r4Scalar**2

  ! energy for BJ damping
  tmp = sqrt(c8Step / c6)              
  e6  = c6     / (RBohr**6 + (rs6*tmp + rs8)**6)
  e8  = c8Step / (RBohr**8 + (rs6*tmp + rs8)**8)

End Subroutine 
!--------------------------------------------------------------------------------------------------------------------------------!    


!________________________________________________________________________________________________________________________________!
Pure Subroutine edisp_Grad( RBohr, e6, e8, e6dr, e8dr )
! compute energy
      
  use Parameters_Module     ,only:  Zero, One, Two, Three

  real(rkp)                               ,intent(in)       ::    RBohr    
  real(rkp)                               ,intent(out)      ::    e6
  real(rkp)                               ,intent(out)      ::    e8
  real(rkp)                               ,intent(out)      ::    e6dr
  real(rkp)                               ,intent(out)      ::    e8dr   
  
  real(rkp)                                                 ::    tmp, c8Step
  
  real(rkp)                                   ,parameter    ::    r2r4Scalar = 2.59361680_rkp                                                                                           
  real(rkp)                                   ,parameter    ::    rs6        = 0.5299_rkp
  real(rkp)                                   ,parameter    ::    rs8        = 2.20_rkp
  real(rkp)                                   ,parameter    ::    c6         = 12.8_rkp

  ! DFT-D3
  c8Step = Three * c6 * r2r4Scalar**2

  ! energy for BJ damping
  tmp = sqrt(c8Step / c6)              
  e6  = c6     / (RBohr**6 + (rs6*tmp + rs8)**6)
  e8  = c8Step / (RBohr**8 + (rs6*tmp + rs8)**8)

  ! calculate gradients
  e6dr =     c6 * (-6.0_rkp * RBohr**5) / (RBohr**6 + (rs6*tmp + rs8)**6)**2
  e8dr = c8Step * (-8.0_rkp * RBohr**7) / (RBohr**8 + (rs6*tmp + rs8)**8)**2

End Subroutine 
!--------------------------------------------------------------------------------------------------------------------------------!       

      
End Module
