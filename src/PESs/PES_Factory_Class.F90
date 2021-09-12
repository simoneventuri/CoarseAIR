! -*-F90-*-
!===============================================================================================================
! 
! Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
! 
! Copyright (C) 2018 Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
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

Module PES_Factory_Class

  use Parameters_Module     ,only:  rkp
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public  ::    PES_Factory_Type

  Type      ::    PES_Factory_Type
  contains
    private
    procedure ,nopass ,public ::  Construct_PES
  End Type

  logical   ,parameter    ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Construct_PES( Input, Atoms, iPES, PES, i_Debug )

  use Input_Class               ,only:    Input_Type
  use Atom_Class                ,only:    Atom_Type
  use PES_Class                 ,only:    PES_Type
  use O3_NN_PES_Class           ,only:    O3_NN_PES_Type
  use N3_NASA_PES_Class         ,only:    N3_NASA_PES_Type                                                                      
  use N3_UMN_PES_Class          ,only:    N3_UMN_PES_Type  
  use N3Half_NASA_PES_Class     ,only:    N3Half_NASA_PES_Type                                                                      
  use CO2_NASA_PES_Class        ,only:    CO2_NASA_PES_Type
  use COAr_NASA_PES_Class       ,only:    COAr_NASA_PES_Type
  use O3_UMN_PES_Class          ,only:    O3_UMN_PES_Type
  use O4_UMN_PES_Class          ,only:    O4_UMN_PES_Type
  use N2O_UMN_PES_Class         ,only:    N2O_UMN_PES_Type
  use N2O_Basel_PES_Class       ,only:    N2O_Basel_PES_Type
  use NO2_Basel_PES_Class       ,only:    NO2_Basel_PES_Type
  use N4_NASA_PES_Class         ,only:    N4_NASA_PES_Type
  use N4_UMN_PES_Class          ,only:    N4_UMN_PES_Type
  use N4_UMN_PIPNN_PES_Class    ,only:    N4_UMN_PIPNN_PES_Type
  use CHN_UIUC_PES_Class        ,only:    CHN_UIUC_PES_Type
  use CNO_Basel_PES_Class       ,only:    CNO_Basel_PES_Type
  use LEPS_PES_Class            ,only:    LEPS_PES_Type
  use NN_PES_Class              ,only:    NN_PES_Type
  use GP_PES_Class              ,only:    GP_PES_Type
  use BNN_PES_Class             ,only:    BNN_PES_Type
  use HO3_NNM_PES_Class         ,only:    HO3_NNM_PES_Type

  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms
  integer                                   ,intent(in)     ::    iPES
  class(PES_Type) ,allocatable              ,intent(out)    ::    PES
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Construct_PES")  !, Active = i_Debug_Loc )
  ! i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Write( "System for the PES to be constructed: Input%System = ", Input%System )
  if (adjustl(trim(Input%PES_Model(iPES))) == 'NN') then
    if (i_Debug_Loc) call Logger%Write( "Constructing a NN_PES_Type object" )
    allocate( NN_PES_Type :: PES )
  elseif (adjustl(trim(Input%PES_Model(iPES))) == 'BNN') then
    if (i_Debug_Loc) call Logger%Write( "Constructing a BNN_PES_Type object" )
    allocate( BNN_PES_Type :: PES )
  elseif (adjustl(trim(Input%PES_Model(iPES))) == 'GP') then
    if (i_Debug_Loc) call Logger%Write( "Constructing a NN_PES_Type object" )
    allocate( GP_PES_Type :: PES )
  elseif (adjustl(trim(Input%PES_Model(iPES))) == 'LEPS') then
    if (i_Debug_Loc) call Logger%Write( "Constructing a LEPS_PES_Type object" )
    allocate( LEPS_PES_Type :: PES )
  else
    select case (Input%System)    

      case('N3', 'NNN', 'NaNbNc')
        !Input%System = 'N3'
        select case (adjustl(trim(Input%PES_Model(iPES))))                                                                     
          case ('NASA')                                                                                                              
            if (i_Debug_Loc) call Logger%Write( "Constructing a N3_NASA_PES_Type object" )
            allocate( N3_NASA_PES_Type :: PES )
          case ('UMN')                                                                                                              
            if (i_Debug_Loc) call Logger%Write( "Constructing a N3_UMN_PES_Type object" )
            allocate( N3_UMN_PES_Type :: PES )
          case ('NASA_FromN4')                                                                                                              
            if (i_Debug_Loc) call Logger%Write( "Constructing a N3Half_NASA_PES_Type object" )
            allocate( N3Half_NASA_PES_Type :: PES )
          case default
            call Error( "PES Model not supported: Input%PES_Model(1) = " // Input%PES_Model(1) )
        end select

      case('CO2', 'O2C', 'COO', 'OCO', 'OOC')
        !Input%System = 'CO2'
        select case (adjustl(trim(Input%PES_Model(iPES))))                                                                         
          case ('NASA_13A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a CO2_NASA_PES_Type object; Surface Selected 1 3A'" )
            allocate( CO2_NASA_PES_Type :: PES )
          case ('NASA_13A2')                                                                                                              
            if (i_Debug_Loc) call Logger%Write( 'Constructing a CO2_NASA_PES_Type object; Surface Selected 1 3A"' )
            allocate( CO2_NASA_PES_Type :: PES )
          case ('NASA_23A2')                                                                                                              
            if (i_Debug_Loc) call Logger%Write( 'Constructing a CO2_NASA_PES_Type object; Surface Selected 2 3A"' )
            allocate( CO2_NASA_PES_Type :: PES )
          case default
            call Error( "PES Model not supported: Input%PES_Model(iPES) = " // Input%PES_Model(iPES) )
        end select

      case('O3', 'OOO')
        !Input%System = 'O3'
        select case (adjustl(trim(Input%PES_Model(iPES))))                                                                     
          case ('UMN_11A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O3_UMN_PES_Type object; Surface Selected 1 1A'" )
            allocate( O3_UMN_PES_Type :: PES )
          case ('UMN_11A2')        
            if (i_Debug_Loc) call Logger%Write( 'Constructing a O3_UMN_PES_Type object; Surface Selected 1 1A"' )
            allocate( O3_UMN_PES_Type :: PES )
          case ('UMN_13A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O3_UMN_PES_Type object; Surface Selected 1 3A'" )
            allocate( O3_UMN_PES_Type :: PES )         
          case ('UMN_13A2')        
            if (i_Debug_Loc) call Logger%Write( 'Constructing a O3_UMN_PES_Type object; Surface Selected 1 3A"' )
            allocate( O3_UMN_PES_Type :: PES )    
          case ('UMN_15A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O3_UMN_PES_Type object; Surface Selected 1 5A'" )
            allocate( O3_UMN_PES_Type :: PES )         
          case ('UMN_15A2')        
            if (i_Debug_Loc) call Logger%Write( 'Constructing a O3_UMN_PES_Type object; Surface Selected 1 5A"' )
            allocate( O3_UMN_PES_Type :: PES )   
          case ('UMN_21A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O3_UMN_PES_Type object; Surface Selected 2 1A'" )
            allocate( O3_UMN_PES_Type :: PES )         
          case ('UMN_23A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O3_UMN_PES_Type object; Surface Selected 2 3A'" )
            allocate( O3_UMN_PES_Type :: PES )
          case ('UMN_25A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O3_UMN_PES_Type object; Surface Selected 2 5A'" )
            allocate( O3_UMN_PES_Type :: PES )      

          case ('NN_11A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O3_NN_PES_Type object; Surface Selected 1 1A'" )
            allocate( O3_NN_PES_Type :: PES )
          case ('NN_11A2')        
            if (i_Debug_Loc) call Logger%Write( 'Constructing a O3_NN_PES_Type object; Surface Selected 1 1A"' )
            allocate( O3_NN_PES_Type :: PES )
          case ('NN_13A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O3_NN_PES_Type object; Surface Selected 1 3A'" )
            allocate( O3_NN_PES_Type :: PES )         
          case ('NN_13A2')        
            if (i_Debug_Loc) call Logger%Write( 'Constructing a O3_NN_PES_Type object; Surface Selected 1 3A"' )
            allocate( O3_NN_PES_Type :: PES )    
          case ('NN_15A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O3_NN_PES_Type object; Surface Selected 1 5A'" )
            allocate( O3_NN_PES_Type :: PES )         
          case ('NN_15A2')        
            if (i_Debug_Loc) call Logger%Write( 'Constructing a O3_NN_PES_Type object; Surface Selected 1 5A"' )
            allocate( O3_NN_PES_Type :: PES )   
          case ('NN_21A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O3_NN_PES_Type object; Surface Selected 2 1A'" )
            allocate( O3_NN_PES_Type :: PES )         
          case ('NN_23A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O3_NN_PES_Type object; Surface Selected 2 3A'" )
            allocate( O3_NN_PES_Type :: PES )
          case ('NN_25A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O3_NN_PES_Type object; Surface Selected 2 5A'" )
            allocate( O3_NN_PES_Type :: PES )      

          case default
            call Error( "PES Model not supported: Input%PES_Model(iPES) = " // Input%PES_Model(iPES) )
        end select

      case('N2O', 'ON2', 'NNO', 'NON', 'ONN')
        !Input%System = 'N2O'
        select case (adjustl(trim(Input%PES_Model(iPES))))                                                                     
          case ('UMN_13A1')        
            if (i_Debug_Loc) call Logger%Write( 'Constructing a N2O_UMN_PES_Type object; Surface Selected 1 3A"' )   
            allocate(N2O_UMN_PES_Type :: PES ) 
          case ('UMN_13A2')        
            if (i_Debug_Loc) call Logger%Write( 'Constructing a N2O_UMN_PES_Type object; Surface Selected 1 3A"' )   
            allocate(N2O_UMN_PES_Type :: PES ) 
          case ('Basel_13A1')        
            if (i_Debug_Loc) call Logger%Write( 'Constructing a N2O_Basel_PES_Type object; Surface Selected 1 3A"' )   
            allocate(N2O_Basel_PES_Type :: PES ) 
          case ('Basel_13A2')        
            if (i_Debug_Loc) call Logger%Write( 'Constructing a N2O_Basel_PES_Type object; Surface Selected 1 3A"' )   
            allocate(N2O_Basel_PES_Type :: PES ) 
          case default
            call Error( "PES Model not supported: Input%PES_Model(iPES) = " // Input%PES_Model(iPES) )
        end select    

      case('NO2', 'O2N', 'NOO', 'OON', 'ONO' )
        !Input%System = 'N2O'
        select case (adjustl(trim(Input%PES_Model(iPES))))                                                                     
          case ('Basel_2A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a N2O_Basel_PES_Type object; Surface Selected 2A'" )   
            allocate(NO2_Basel_PES_Type :: PES )
          case ('Basel_2A2')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a N2O_Basel_PES_Type object; Surface Selected 2A''" )   
            allocate(NO2_Basel_PES_Type :: PES ) 
          case ('Basel_4A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a N2O_Basel_PES_Type object; Surface Selected 4A'" )   
            allocate(NO2_Basel_PES_Type :: PES ) 
          case default
            call Error( "PES Model not supported: Input%PES_Model(iPES) = " // Input%PES_Model(iPES) )
        end select    

      case('COAr', 'CArO', 'ArCO', 'ArOC', 'OCAr', 'OArC')
        !Input%System = 'COAr'
        select case (Input%PES_Model(iPES))                                                                                  
          case ('NASA')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a COAr_NASA_PES_Type object" )
            allocate( COAr_NASA_PES_Type :: PES )
          case default
            call Error( "PES Model not supported: Input%PES_Model(iPES) = " // Input%PES_Model(iPES) )
        end select

      case('CHN', 'CNH', 'NCH', 'NHC', 'HCN', 'HNC')
        !Input%System = 'CHN'
        select case (Input%PES_Model(iPES))                                                                                  
          case ('UIUC')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a CHN_UIUC_PES_Type object" )
            allocate( CHN_UIUC_PES_Type :: PES )
          case default
            call Error( "PES Model not supported: Input%PES_Model(iPES) = " // Input%PES_Model(iPES) )
        end select

       case('CON', 'CNO', 'NCO', 'NOC', 'OCN', 'ONC')
        !Input%System = 'CON'
        select case (Input%PES_Model(iPES))                                                                                  
          case ('Basel_2A1')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a CNO_Basel_PES_Type object; Surface Selected 2A'" )
            allocate( CNO_Basel_PES_Type :: PES )
          case ('Basel_2A2')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a CNO_Basel_PES_Type object; Surface Selected 2A''" )
            allocate( CNO_Basel_PES_Type :: PES )
          case ('Basel_4A2')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a CNO_Basel_PES_Type object; Surface Selected 4A''" )
            allocate( CNO_Basel_PES_Type :: PES )
          case default
            call Error( "PES Model not supported: Input%PES_Model(iPES) = " // Input%PES_Model(iPES) )
        end select

      case('N4', 'N2N2', 'NNNN', 'NaNbNcNd')
        !Input%System = 'N4'
        select case (adjustl(trim(Input%PES_Model(iPES))))                                                                     
          case ('NASA')                                                                                                              
            if (i_Debug_Loc) call Logger%Write( "Constructing a N4_NASA_PES_Type object" )
            allocate( N4_NASA_PES_Type :: PES )
          case ('UMN')                                                                                                              
            if (i_Debug_Loc) call Logger%Write( "Constructing a N4_UMN_PES_Type object" )
            allocate( N4_UMN_PES_Type :: PES )

          case ('UMN_PIPNN')                                                                                                              
            if (i_Debug_Loc) call Logger%Write( "Constructing a N4_UMN_PIPNN_PES_Type object" )
            allocate( N4_UMN_PIPNN_PES_Type :: PES )
            PES%UseSurface = [1,1,1]
            if (i_Debug_Loc) call Logger%Write( "PES%UseSurface = ", PES%UseSurface )
          case ('UMN_PIPNN_A')                                                                                                              
            if (i_Debug_Loc) call Logger%Write( "Constructing a N4_UMN_PIPNN_PES_Type object" )
            allocate( N4_UMN_PIPNN_PES_Type :: PES )
            PES%UseSurface = [1,0,0]
            if (i_Debug_Loc) call Logger%Write( "PES%UseSurface = ", PES%UseSurface )
          case ('UMN_PIPNN_B')                                                                                                              
            if (i_Debug_Loc) call Logger%Write( "Constructing a N4_UMN_PIPNN_PES_Type object" )
            allocate( N4_UMN_PIPNN_PES_Type :: PES )
            PES%UseSurface = [0,1,0]
            if (i_Debug_Loc) call Logger%Write( "PES%UseSurface = ", PES%UseSurface )
          case ('UMN_PIPNN_C')                                                                                                              
            if (i_Debug_Loc) call Logger%Write( "Constructing a N4_UMN_PIPNN_PES_Type object" )
            allocate( N4_UMN_PIPNN_PES_Type :: PES )
            PES%UseSurface = [0,0,1]
            if (i_Debug_Loc) call Logger%Write( "PES%UseSurface = ", PES%UseSurface )

          case default
            call Error( "PES Model not supported: Input%PES_Model(1) = " // Input%PES_Model(1) )
        end select

      case('O4', 'O2O2', 'OOOO', 'OaObOcOd')
        !Input%System = 'N4'
        select case (adjustl(trim(Input%PES_Model(iPES))))                                                                     
          case ('UMN_Singlet')        
            if (i_Debug_Loc) call Logger%Write( "Constructing a O4_UMN_PES_Type object" )
            allocate( O4_UMN_PES_Type :: PES )
          case ('UMN_Triplet')          
            if (i_Debug_Loc) call Logger%Write( "Constructing a O4_UMN_PES_Type object" )
            allocate( O4_UMN_PES_Type :: PES )
          case ('UMN_Quintet')      
            if (i_Debug_Loc) call Logger%Write( "Constructing a O4_UMN_PES_Type object" )
            allocate( O4_UMN_PES_Type :: PES )
          case default
            call Error( "PES Model not supported: Input%PES_Model(1) = " // Input%PES_Model(1) )
        end select

      case('HO3')
        select case (adjustl(trim(Input%PES_Model(iPES))))
          case ('NNM')
            if (i_Debug_Loc) call Logger%Write( "Constructing a HO3_NNM_Pes_Type object")
            allocate( HO3_NNM_PES_Type :: PES )
          case default
            call Error( "PES Model not supported: Input%PES_Model(1) = " // Input%PES_Model(1) )
      end select

      case default
        call Error( "PES not supported: Input%System = " // Input%System )
    end select
  end if
  
  
  if (i_Debug_Loc) call Logger%Write( "Calling PES%Initialize" )
  call PES%Initialize( Input, Atoms, iPES, i_Debug )

  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
