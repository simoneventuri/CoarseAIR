Module Logger_Variable_Class

  implicit none

  private
  public        ::    Logger_Variable_Type

  Type                                          ::    Logger_Variable_Type
    integer                                     ::    VarDim                                          !> Dimension of the variable: 0, 1, 2
    class(*)      ,allocatable                  ::    Var0d, Var1d(:), Var2d(:,:)                     !> Value of the variable (The exact value ised is given by VarDim)
    character(:)  ,allocatable                  ::    VarFmt                                          !> The Format of the variable
    character(:)  ,allocatable                  ::    VarStr                                          !> The string of the variable
  contains                                                                                                      ! Starting derived-type procedure declarations
    generic     ,public                         ::    ConvertToString => ConvertToString_0D, ConvertToString_1D
    procedure   ,private                        ::    ConvertToString_0D
    procedure   ,private                        ::    ConvertToString_1D
    procedure   ,private                        ::    ConvertToStrings_1D
  End Type

  integer       ,parameter                      ::    Indentation_Step          =     2
  character(*)  ,parameter                      ::    Default_Name_LogFile      =     "logfile.log"
  character(*)  ,parameter                      ::    Default_Logical_Format    =     "l1"
  character(*)  ,parameter                      ::    Default_Integer_Format    =     "i0"
  character(*)  ,parameter                      ::    Default_Real_Format       =     "g0"
  character(*)  ,parameter                      ::    Default_Character_Format  =     "a"
  character(*)  ,parameter                      ::    Default_Spacing2_Format   =     "3x"

  logical ,parameter :: Debug=.False.

  Interface             ConvertToString
    Module Procedure    ConvertToString_0D
    Module Procedure    ConvertToString_1D
    Module Procedure    ConvertToStrings_1D
  End Interface

  Interface             Convert___To___String
    Module Procedure    Convert_Logical_To_String_0D , Convert_Logical_To_String_1D , Convert_Logical_To_Strings_1D
    Module Procedure    Convert_Integer_To_String_0D , Convert_Integer_To_String_1D , Convert_Integer_To_Strings_1D
    Module Procedure    Convert_Real4_To_String_0D   , Convert_Real4_To_String_1D   , Convert_Real4_To_Strings_1D
    Module Procedure    Convert_Real8_To_String_0D   , Convert_Real8_To_String_1D   , Convert_Real8_To_Strings_1D
    Module Procedure    Convert_String_To_String_0D  , Convert_String_To_String_1D  , Convert_String_To_Strings_1D
  End Interface

  Interface             Logger_Variable_Type
    Module Procedure    Set_LogVar_0d
  End Interface

!   Interface             Set_Vector_Format
!     Module Procedure    Set_Vector_Format_0D, Set_Vector_Format_1D
!   End Interface

  contains


!   Usage:
!   Var1  = Logger_Variable_Type( V1, F1, Fd )
! !                               |   |   |
! !                               |   |   v
! !                               |   |   Default format for all variables
! !                               |   v
! !                               |   Format for the current variable
! !                               v
! !                               Variable to be converted
! !
! the call to:  Var1    =   Logger_Variable_Type( V1, F1, Fv )
! set the components:
!  - Var1%Str         String holding the value of the variable
Subroutine Write_4xV0(  This,                                               &
                        V1,  V2,  V3,  V4,                                  &
                        Unused, i_Prefix, NewLine, Advance, Backspace, Status,         &
                        F1,  F2,  F3,  F4                                   )
  class(Logger_Type)                                    ,intent(inout)  ::    This                            !< Passed-object dummy argument corresponding to the Logger object
  class(*)                                              ,intent(in)     ::    V1,  V2,  V3,  V4
  logical       ,optional   ,dimension(:,:,:,:)         ,intent(in)     ::    Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
  logical       ,optional                               ,intent(in)     ::    i_Prefix                        !< Indicator of the prefix presence
  logical       ,optional                               ,intent(in)     ::    NewLine                         !< Indicator whether or not a black line must be written before writing the variables
  logical       ,optional                               ,intent(in)     ::    Advance                         !< Indicator whether or not the line should be advanced
  logical       ,optional                               ,intent(in)     ::    Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
  integer       ,optional                               ,intent(out)    ::    Status
  character(*)  ,optional                               ,intent(in)     ::    Fv, F1,  F2,  F3,  F4
  character(:)  ,allocatable                                            ::    S1,  S2,  S3,  S4
  character(:)  ,allocatable                                            ::    Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::    Fmt
  character(:)  ,allocatable                                            ::    Local_Prefix                    ! Local prefix variable
  integer                                                               ::    Local_Status, ios


  type(Logger_Variable_Type)  ::  Var1, Var2, Var3, Var4

  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0

  Var1    =   Logger_Variable_Type( V1, F1, Fv ) ! ; call ConvertToString( V1 , S1 , F1 , Status=ios ); if ( ios /= 0 ) Local_Status = ios
  Var2    =   Logger_Variable_Type( V2, F2, Fv ) ! ; call ConvertToString( V2 , S2 , F2 , Status=ios ); if ( ios /= 0 ) Local_Status = ios
  Var3    =   Logger_Variable_Type( V3, F3, Fv ) ! ; call ConvertToString( V3 , S3 , F3 , Status=ios ); if ( ios /= 0 ) Local_Status = ios
  Var4    =   Logger_Variable_Type( V4, F4, Fv ) ! ; call ConvertToString( V4 , S4 , F4 , Status=ios ); if ( ios /= 0 ) Local_Status = ios

  VarAll  =   Logger_Variable_Type( [Var1, Var2,Var3,Var4] )    ! Fmt           =       This%Get_Prefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
!   VarAll%Str
!   VarAll%Fmt


  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status

  write( This%Unit, VarAll%Fmt, Advance=Adv, IOStat=ios ) VarAll%Str
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Subroutine


! This procedure opens a log file.
! If the unit argument is provided, then the Log file unit is set to it... it's still under development...
! One should check that the file is actually opened
Subroutine Set_LogVar_0d( This, Var, Fmt, FmfDef )

  class(Logger_Variable_Type)                           ,intent(out)    ::    This                            !< Passed-object dummy argument corresponding
  class(*)                                              ,intent(in)     ::    Var
  character(*)  ,optional ,target                       ,intent(in)     ::    Fmt     ! Format specific for current variable (use this one if present)
  character(*)  ,optional ,target                       ,intent(in)     ::    FmfDef  ! Format common to a set of variables  (use this one if the specific format is not present )

  character(*)  ,pointer                                                ::    FmfLoc

  FmfLoc    =>    null()
  if ( present(FmfDef) ) FmfLoc => FmfDef
  if ( present(Fmt   ) ) FmfLoc => Fmt

    call ConvertToStrings_1D( Var, This%VarStr, FmfLoc, Status )

  Type                                          ::    Logger_Variable_Type
    integer                                     ::    VarDim                                          !> Dimension of the variable: 0, 1, 2
    class(*)      ,allocatable                  ::    Var0d, Var1d(:), Var2d(:,:)                     !> Value of the variable (The exact value ised is given by VarDim)
    character(:)  ,allocatable                  ::    VarFmt                                          !> The Format of the variable
    character(:)  ,allocatable                  ::    VarStr                                          !> The string of the variable
  contains
Subroutine ConvertToString_0D( Variable, String, Format, Status )
  class(*)                                              ,intent(in)     ::    Variable
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  integer                                                               ::    Local_Status
  select type (Variable)

End Subroutine

Pure Function Set_Generic_Scalar_Format( Default_Format, Optional_Format ) result(Format)
  character(*)                                          ,intent(in)     ::    Default_Format
  character(*)  ,optional                               ,intent(in)     ::    Optional_Format
  character(:)  ,allocatable                                            ::    Format
  Format        =       Default_Format
  if ( present(Optional_Format) ) then
    if ( Is_Valid_Format(Optional_Format) ) then              ! PROBLEM: This is specific to Character variable: Pass a procedure
      Format = Optional_Format
    else        ! What should we do here ??? (1) an error or (2) something a bit less critical like taking the default format
    end if
  end if
End Function

Pure Function Set_Generic_Vector_Format( Default_Format, Optional_Format ) result(Format)
  character(*)                                          ,intent(in)     ::    Default_Format
  character(*)  ,optional                               ,intent(in)     ::    Optional_Format
  character(:)  ,allocatable                                            ::    Format
  Format        = "*(" // Default_Format // ",3x)"
  if ( present(Optional_Format) ) then
    if ( Is_Valid_Format(Optional_Format) ) then              ! PROBLEM: This is specific to Character variable: Pass a procedure
      Format = "*(" // Optional_Format // ",3x)"
    else        ! What should we do here ??? (1) an error or (2) something a bit less critical like taking the default format
    end if
  end if
End Function

! TODO: Implement the algo telling if a given string is a valid fortran format.
Pure Elemental Function Is_Valid_Format( String ) result(Valid_Format)
  character(*)                                          ,intent(in)     ::    String
  logical                                                               ::    Valid_Format
  character(:)  ,allocatable                                            ::    String_Loc
  String_Loc = String
  Valid_Format  =       .True.
End Function


! **************************************************************************************************************
! **************************************************************************************************************
!                                           TOOLS
! **************************************************************************************************************
! **************************************************************************************************************

! This function set and check a valid open status.
! If a valid optional open status is passed, then it is set other wise the default open status is taken
Function Get_OptOrDef_Value( Default_Value, Valid_Values, Optional_Value ) result( Output_Value )
  character(*)                                          ,intent(in)     ::    Default_Value                   !< Default value used if no optional value
  character(*)  ,dimension(:)                           ,intent(in)     ::    Valid_Values                    !< Valid values used to check validity of optional values if present
  character(*)  ,optional                               ,intent(in)     ::    Optional_Value                  !< Optional values used if present and valid
  character(:)  ,allocatable                                            ::    Output_Value                    !< Output values
  Output_Value  =       Default_Value
  if ( present(Optional_Value) ) then
    if ( Is_Valid(Optional_Value,Valid_Values) ) then
      Output_Value = Optional_Value
    else
      call Error_Unvalid_Value( "", Optional_Value, Valid_Values )
    end if
  end if
End Function

Pure Function Is_Valid( Value, Valid_Values ) result(Valid)
  implicit none
  character(*)                                  ,intent(in)     ::    Value                                   !< Value to be checked for validity
  character(*)          ,dimension( : )         ,intent(in)     ::    Valid_Values                            !< Valid values used for validity check
  logical                                                       ::    Valid                                   !< Indicator of input object validity
  integer                                                       ::    i                                       ! Index of valid strings
  Valid         =       .False.                                                                                 ! Initialization of the object validity indicator to false
  do i = 1,size(Valid_Values)                                                                                   ! Loop on all valid strings
    if ( trim(Value) == trim(Valid_Values(i)) ) Valid = .True.                                                  ! If the object if found in the list of valid strings, then setting validity indicator to true
  end do                                                                                                        ! End do loop on valid strings
End Function

Subroutine Add_Element_To_Array( Element, Array )
  implicit none
  character(*)                                          ,intent(in)     ::    Element
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::    Array
  integer                                                               ::    Length
  character(:)  ,dimension(:)   ,allocatable                            ::    Array_tmp
#ifdef WORKAROUND_GFORTRAN_SOURCE_ALLOCATION
  integer                                                               ::    i
#endif
  if ( .not. allocated(Array) ) allocate( character(0) :: Array(0) )
  Length        =       max( len(Array), len(Element) )
  allocate( character(Length) :: Array_tmp(size(Array)+1) )
#ifdef WORKAROUND_GFORTRAN_SOURCE_ALLOCATION
!   ------------------------------
  do i = 1,size(Array)
    Array_tmp(i)    =       Array(i)
  end do
!   ------------------------------
#else
  Array_tmp(1:size(Array))    =       Array     ! COMPILER_BUG:GFORTRAN
#endif
  Array_tmp(size(Array)+1)    =       Element
  call move_alloc( Array_tmp, Array )
End Subroutine

Subroutine Remove_Element_From_Array( Array )
  implicit none
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::    Array
  integer                                                               ::    NElements
  integer                                                               ::    Length
  character(:)  ,dimension(:)   ,allocatable                            ::    Array_tmp
#ifdef WORKAROUND_GFORTRAN_SOURCE_ALLOCATION
  integer                                                       ::    n, i
#endif
  if ( .not. allocated(Array) ) allocate( character(0) :: Array(0) )
  NElements     =       size(Array)
  Length        =       len(Array(1:NElements-1))
#ifdef WORKAROUND_GFORTRAN_SOURCE_ALLOCATION
!   ------------------------------
  n = NElements-1
  allocate( Character(Length) :: Array_tmp(n) )
!   Array_tmp = Array(1:NElements-1)    ! COMPILER_BUG:GFORTRAN
  do i = 1,n
    Array_tmp(i) = Array(i)
  end do
!   ------------------------------
#else
  allocate( Array_tmp, source=Array(1:NElements-1) )
#endif
  call move_alloc( Array_tmp, Array )
End Subroutine

Function VecTrim( Input_String ) result(Output_String)
  character(*)  ,dimension(:)                           ,intent(in)     ::    Input_String
  character(:)  ,dimension(:)           ,allocatable                    ::    Output_String
  integer                                                               ::    i
  integer                                                               ::    Length
  integer                                                               ::    Size1
  Length        =       Max_Len_Trim(Input_String)
  Size1         =       size(Input_String,1)
  allocate( character(Length) :: Output_String(Size1) )
  do i = 1,Size1
    Output_String(i)    =       trim( Input_String(i) )
  end do
End Function

Pure Function Max_Len_Trim( Strings ) result( Length )
  character(*)  ,dimension(:)                   ,intent(in)             ::    Strings                         !< Array of character string
  integer                                                               ::    Length                          !< Maximum length without trailling blanks along all elements of the input string array
  integer                                                               ::    i                               ! Index of string' elements
  Length        =       0                                                                                       ! Initialization of maximum length of string
  do i = 1,size(Strings,1)                                                                                      ! Loop on all elements
    Length      =       max( Length, len_trim(Strings(i)) )                                                     ! Setting the maximum length
  end do                                                                                                        ! End loop on all elements
End Function



! **************************************************************************************************************
! **************************************************************************************************************
!                                       CONVERTION PROCEDURES
! **************************************************************************************************************
! **************************************************************************************************************

Subroutine ConvertToString_0D( Variable, String, Format, Status )
  class(*)                                              ,intent(in)     ::    Variable
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  integer                                                               ::    Local_Status
  select type (Variable)
    type is (logical);      call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    type is (integer);      call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    type is (real(4));      call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    type is (real(8));      call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    type is (character(*)); call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    class default;       ! Error
      String        =       "?"
      Local_Status  =       -1
  end select
  if ( present(Status) ) Status = Local_Status
End Subroutine

Subroutine ConvertToString_1D( Variable, String, Format, Status )
  class(*)      ,dimension(:)                           ,intent(in)     ::    Variable
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  integer                                                               ::    Local_Status
  select type (Variable)
    type is (logical);      call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    type is (integer);      call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    type is (real(4));      call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    type is (real(8));      call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    type is (character(*)); call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    class default;       ! Error
      String        =       "?"
      Local_Status  =       -1
  end select
  if ( present(Status) ) Status = Local_Status
End Subroutine

Subroutine ConvertToStrings_1D( Variable, String, Format, Status )
  class(*)      ,dimension(:)                           ,intent(in)     ::    Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  integer                                                               ::    Local_Status
  select type (Variable)
    type is (logical);      call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    type is (integer);      call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    type is (real(4));      call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    type is (real(8));      call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    type is (character(*)); call Convert___To___String( Variable, String, Format=Format, Status=Local_Status )
    class default;       ! Error
      String        =       "?"
      Local_Status  =       -1
  end select
  if ( present(Status) ) Status = Local_Status
End Subroutine

Subroutine Convert_Logical_To_String_0D( Variable, String, Format, Status )
  logical                                               ,intent(in)     ::    Variable
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)                                                      ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios
  Local_Format    =       Set_Generic_Scalar_Format( Default_Logical_Format, Format )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Integer_To_String_0D( Variable, String, Format, Status )
  integer                                               ,intent(in)     ::    Variable
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)                                                      ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios
  Local_Format    =       Set_Generic_Scalar_Format( Default_Integer_Format, Format )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Real4_To_String_0D( Variable, String, Format, Status )
  real(4)                                               ,intent(in)     ::    Variable                             !< Real number to be converted into a string
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)                                                      ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios
  Local_Format    =       Set_Generic_Scalar_Format( Default_Real_Format, Format )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Real8_To_String_0D( Variable, String, Format, Status )
  real(8)                                               ,intent(in)     ::    Variable                             !< Real number to be converted into a string
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)                                                      ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios
  Local_Format    =       Set_Generic_Scalar_Format( Default_Real_Format, Format )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_String_To_String_0D( Variable, String, Format, Status )
  character(*)                                          ,intent(in)     ::    Variable                             !< Real number to be converted into a string
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)                                                      ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios
  if ( present(Format) ) then
    Local_Format    =       Set_Generic_Scalar_Format( Default_Character_Format, Format )
    Local_Format    =       "(" // Local_Format // ")"
    write( Long_String , Local_Format , iostat=ios ) Variable
    if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
    String          =       trim(Long_String)
  else
    String          =       Variable
    ios             =       0
  end if
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Logical_To_String_1D( Variable, String, Format, Status )
  logical       ,dimension(:)                           ,intent(in)     ::    Variable
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)                                                      ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios
  Local_Format    =       Set_Generic_Vector_Format( Default_Logical_Format, Format )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Integer_To_String_1D( Variable, String, Format, Status )
  integer       ,dimension(:)                           ,intent(in)     ::    Variable
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)                                                      ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios
  Local_Format    =       Set_Generic_Vector_Format( Default_Integer_Format, Format )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Real4_To_String_1D( Variable, String, Format, Status )
  real(4)       ,dimension(:)                           ,intent(in)     ::    Variable
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)                                                      ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios
  Local_Format    =       Set_Generic_Vector_Format( Default_Real_Format, Format )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Real8_To_String_1D( Variable, String, Format, Status )
  real(8)       ,dimension(:)                           ,intent(in)     ::    Variable
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)                                                      ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios
  Local_Format    =       Set_Generic_Vector_Format( Default_Real_Format, Format )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_String_To_String_1D( Variable, String, Format, Status )
  character(*)  ,dimension(:)                           ,intent(in)     ::    Variable
  character(:)  ,allocatable                            ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)                                                      ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios
  integer                                                               ::    i
  if ( present(Format) ) then
    Local_Format    =       Set_Generic_Vector_Format( Default_Character_Format, Format )
    Local_Format    =       "(" // Local_Format // ")"
    write( Long_String , Local_Format , iostat=ios ) Variable
    if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
    String          =       trim(Long_String)
  else
    String          =       ""
    do i = 1,size(Variable)
      String        =       String // Variable(i) // "   "
    end do
    String          =       trim(String)
    ios             =       0
  end if
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Logical_To_Strings_1D( Variable, String, Format, Status )
  logical       ,dimension(:)                           ,intent(in)     ::    Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)  ,dimension( size(Variable) )                        ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios, Stat, i
  Local_Format    =       Set_Generic_Scalar_Format( Default_Logical_Format, Format )
  Local_Format    =       "(" // Local_Format // ")"
  Stat            =       0
  do i = 1,size(Variable)
    write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
    if ( ios /= 0 ) then
      write(Long_String(i) ,"(g0)") Variable(i)
      Stat        =       ios
    end if
  end do
  allocate( String , source = VecTrim(Long_String), Stat=ios )
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Subroutine Convert_Integer_To_Strings_1D( Variable, String, Format, Status )
  integer       ,dimension(:)                           ,intent(in)     ::    Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)  ,dimension( size(Variable) )                        ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios, Stat, i
  Local_Format    =       Set_Generic_Scalar_Format( Default_Integer_Format, Format )
  Local_Format    =       "(" // Local_Format // ")"
  Stat            =       0
  do i = 1,size(Variable)
    write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
    if ( ios /= 0 ) then
      write(Long_String(i) ,"(g0)") Variable(i)
      Stat        =       ios
    end if
  end do
  allocate( String , source = VecTrim(Long_String), Stat=ios )
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Subroutine Convert_Real4_To_Strings_1D( Variable, String, Format, Status )
  real(4)       ,dimension(:)                           ,intent(in)     ::    Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)  ,dimension( size(Variable) )                        ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios, Stat, i
  Local_Format    =       Set_Generic_Scalar_Format( Default_Real_Format, Format )
  Local_Format    =       "(" // Local_Format // ")"
  Stat            =       0
  do i = 1,size(Variable)
    write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
    if ( ios /= 0 ) then
      write(Long_String(i) ,"(g0)") Variable(i)
      Stat        =       ios
    end if
  end do
  allocate( String , source = VecTrim(Long_String), Stat=ios )
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Subroutine Convert_Real8_To_Strings_1D( Variable, String, Format, Status )
  real(8)       ,dimension(:)                           ,intent(in)     ::    Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)  ,dimension( size(Variable) )                        ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios, Stat, i
  Local_Format    =       Set_Generic_Scalar_Format( Default_Real_Format, Format )
  Local_Format    =       "(" // Local_Format // ")"
  Stat            =       0
  do i = 1,size(Variable)
    write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
    if ( ios /= 0 ) then
      write(Long_String(i) ,"(g0)") Variable(i)
      Stat        =       ios
    end if
  end do
  allocate( String , source = VecTrim(Long_String), Stat=ios )
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Subroutine Convert_String_To_Strings_1D( Variable, String, Format, Status )
  character(*)  ,dimension(:)                           ,intent(in)     ::    Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::    String
  character(*)  ,optional                               ,intent(in)     ::    Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(10000)  ,dimension( size(Variable) )                        ::    Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::    Local_Format
  integer                                                               ::    ios, Stat, i
  if ( present(Format) ) then
    Local_Format    =       Set_Generic_Scalar_Format( Default_Real_Format, Format )
    Local_Format    =       "(" // Local_Format // ")"
    Stat            =       0
    do i = 1,size(Variable)
      write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
      if ( ios /= 0 ) then
        write(Long_String(i) ,"(g0)") Variable(i)
        Stat        =       ios
      end if
    end do
    allocate( String , source = VecTrim(Long_String), Stat=ios )
    if ( ios /= 0 ) Stat = ios
  else
    allocate( String , source = VecTrim(Variable), Stat=Stat )
  end if
  if ( present(Status) ) Status = Stat
End Subroutine


! **************************************************************************************************************
! **************************************************************************************************************
!                                   PROCEDURES RELATED TO FORMATS
! **************************************************************************************************************
! **************************************************************************************************************

Subroutine Get_Integer_Format( Variable, Input_Format, Output_Format, Status )
  integer                                               ,intent(in)     ::    Variable
  character(*)  ,optional                               ,intent(in)     ::    Input_Format
  character(:)  ,allocatable                            ,intent(out)    ::    Output_Format
  integer       ,optional                               ,intent(out)    ::    Status
  character(:)  ,allocatable                                            ::    String
  integer                                                               ::    NDigits
  integer                                                               ::    ios
  if ( present(Input_Format) ) then
    Output_Format   =     Input_Format
    ios             =     0
  else
    NDigits         =       floor( log10( real( abs(Variable) ) ) ) + 1                                                     ! Getting the number of digits of the input integer number
    call Convert___To___String( NDigits, String, Status=ios )
    Output_Format   =       "(i" // trim(String) // ")"
  end if
  if ( present(Status) ) Status = ios
End Subroutine

! Function Set_Scalar_Format_0D( Variable, Optional_Format ) result(Format)
!   class(*)                                              ,intent(in)     ::    Variable
!   character(*)  ,optional                               ,intent(in)     ::    Optional_Format
!   character(:)  ,allocatable                                            ::    Format
!   select type (Variable)
!     type is (logical);          Format = Set_Generic_Scalar_Format( Default_Logical_Format,    Optional_Format )
!     type is (integer);          Format = Set_Generic_Scalar_Format( Default_Integer_Format,    Optional_Format )
!     type is (real(4));          Format = Set_Generic_Scalar_Format( Default_Real_Format,       Optional_Format )
!     type is (real(8));          Format = Set_Generic_Scalar_Format( Default_Real_Format,       Optional_Format )
!     type is (character(*));     Format = Set_Generic_Scalar_Format( Default_Character_Format,  Optional_Format )
!     class default       ! Error
!   end select
! End Function
!
! Function Set_Scalar_Format_1D( Variable, Optional_Format ) result(Format)
!   class(*)      ,dimension(:)                           ,intent(in)     ::    Variable
!   character(*)  ,optional                               ,intent(in)     ::    Optional_Format
!   character(:)  ,allocatable                                            ::    Format
!   select type (Variable)
!     type is (logical);          Format = Set_Generic_Scalar_Format( Default_Logical_Format,    Optional_Format )
!     type is (integer);          Format = Set_Generic_Scalar_Format( Default_Integer_Format,    Optional_Format )
!     type is (real(4));          Format = Set_Generic_Scalar_Format( Default_Real_Format,       Optional_Format )
!     type is (real(8));          Format = Set_Generic_Scalar_Format( Default_Real_Format,       Optional_Format )
!     type is (character(*));     Format = Set_Generic_Scalar_Format( Default_Character_Format,  Optional_Format )
!     class default       ! Error
!   end select
! End Function

! Function Set_Vector_Format_0D( Variable, Optional_Format ) result(Format)
!   class(*)                                              ,intent(in)     ::    Variable
!   character(*)  ,optional                               ,intent(in)     ::    Optional_Format
!   character(:)  ,allocatable                                            ::    Format
!   select type (Variable)
!     type is (logical);          Format = Set_Generic_Vector_Format( Default_Logical_Format,    Optional_Format )
!     type is (integer);          Format = Set_Generic_Vector_Format( Default_Integer_Format,    Optional_Format )
!     type is (real(4));          Format = Set_Generic_Vector_Format( Default_Real_Format,       Optional_Format )
!     type is (real(8));          Format = Set_Generic_Vector_Format( Default_Real_Format,       Optional_Format )
!     type is (character(*));     Format = Set_Generic_Vector_Format( Default_Character_Format,  Optional_Format )
!     class default       ! Error
!   end select
! End Function
!
! Function Set_Vector_Format_1D( Variable, Optional_Format ) result(Format)
!   class(*)      ,dimension(:)                           ,intent(in)     ::    Variable
!   character(*)  ,optional                               ,intent(in)     ::    Optional_Format
!   character(:)  ,allocatable                                            ::    Format
!   select type (Variable)
!     type is (logical);          Format = Set_Generic_Vector_Format( Default_Logical_Format,    Optional_Format )
!     type is (integer);          Format = Set_Generic_Vector_Format( Default_Integer_Format,    Optional_Format )
!     type is (real(4));          Format = Set_Generic_Vector_Format( Default_Real_Format,       Optional_Format )
!     type is (real(8));          Format = Set_Generic_Vector_Format( Default_Real_Format,       Optional_Format )
!     type is (character(*));     Format = Set_Generic_Vector_Format( Default_Character_Format,  Optional_Format )
!     class default       ! Error
!   end select
! End Function


End Module
