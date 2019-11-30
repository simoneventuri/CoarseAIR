Module String_Module

  use, intrinsic :: ISO_Fortran_Env ,only: Output_Unit

  implicit none

  private
  public        ::      Parse, Split
  public        ::      RemoveLeftChar, RemoveSpace
  public        ::      LowerCase, UpperCase
  public        ::      Tab2Space
  public        ::      Reorder
  public        ::      Compact
  public        ::      Presence
  public        ::      FindPos
  public        ::      Count_Presence
  public        ::      Convert_To_String
  public        ::      Convert_To_Real
  public        ::      Convert_To_Integer
  public        ::      Convert_To_Logical
  public        ::      Add_Line_To_String
  public        ::      Get_Number_Of_Digits
  public        ::      Remove_Duplicate
  public        ::      Remove_Trailing_Zeros
  public        ::      Max_Len_Trim, VecTrim
  public        ::      Replace_Character
  public        ::      Add_Element_If_Absent
  public        ::      Get_Characters_AtRightOf_Numbers
  public        ::      Get_Numbers_AtLeftOf_Characters
  public        ::      Set_Enhanced_Name
  public        ::      Remove_Last_Directory_From_Path
  public        ::      Set_Length
  public        ::      Inline
  public        ::      Get_SubString, Get_SubString_Indexes
  public        ::      Exchange_Element_In_Array
  public        ::      Add_Prefix
  public        ::      Add_Suffix
  public        ::      Set_LogUnit
  public        ::      Convert_Ratio
  public        ::      Real_Format_From_Length
  public        ::      Is_A_Number
  public        ::      Remove_Character
  public        ::      Equal
  public        ::      Is_Numeric
  public        ::      Is_Letter

  Interface             UpperCase
    Module Procedure    UpperCase_0D
    Module Procedure    UpperCase_1D
  End Interface

  Interface             Get_Numbers_AtLeftOf_Characters
    Module Procedure    Get_Numbers_AtLeftOf_Characters_0D
    Module Procedure    Get_Numbers_AtLeftOf_Characters_1D
  End Interface

  Interface             Reorder
    Module Procedure    Reorder_Character, Reorder_Real_8, Reorder_Real_4
  End Interface

  Interface             Presence
    Module Procedure    Presence_Character_0D, Presence_Integer_0D, Presence_Real4_0D, Presence_Real8_0D!, Presence_Character_0D_Deferred_Length
    Module Procedure    Presence_Character_1D, Presence_Integer_1D, Presence_Real4_1D, Presence_Real8_1D
  End Interface

  Interface             Count_Presence
    Module Procedure    Count_Presence_Character_0D, Count_Presence_Integer_0D, Count_Presence_Real4_0D, Count_Presence_Real8_0D
    Module Procedure    Count_Presence_Character_1D, Count_Presence_Integer_1D, Count_Presence_Real4_1D, Count_Presence_Real8_1D
  End Interface

  Interface             FindPos
    Module Procedure    FindPos_0D, FindPos_1D
  End Interface

  Interface             Convert_To_String
    Module Procedure    Logical_To_String_0D
    Module Procedure    Integer_To_String_0D
    Module Procedure    Real4_To_String_0D
    Module Procedure    Real8_To_String_0D
    Module Procedure    Real16_To_String_0D
    Module Procedure    String_To_String_0D
    Module Procedure    Integer1D_To_String0D
  End Interface

  Interface             Convert_To_Real
    Module Procedure    String_To_Real8_0D
    Module Procedure    String_To_Real8_1D
  End Interface

  Interface             Convert_To_Integer
    Module Procedure    String_To_Integer_0D, String_To_Integer_1D
  End Interface

  Interface             Convert_To_Logical
    Module Procedure    String_To_Logical_0D
  End Interface

  Interface             Add_Line_To_String
    Module Procedure    Add_Line_To_String, Add_Lines_To_String
  End Interface

  Interface             Max_Len_Trim
    Module Procedure    Max_Len_Trim_0D, Max_Len_Trim_1D, Max_Len_Trim_2D
  End Interface

  Interface             VecTrim
    Module Procedure    VecTrim_0D, VecTrim_1D, VecTrim_2D
  End Interface

  Interface             Get_Characters_AtRightOf_Numbers
    Module Procedure    Get_Characters_AtRightOf_Numbers_0d
    Module Procedure    Get_Characters_AtRightOf_Numbers_1d
  End Interface

  Interface             Add_Element_If_Absent
    Module Procedure    Add_Element_If_Absent_Integer
    Module Procedure    Add_Element_If_Absent_C0,       Add_Element_If_Absent_C1
  End Interface

  Interface             Replace_Character
    Module Procedure    Replace_Character
    Module Procedure    Replace_Characters
  End Interface

  Interface             Inline
    Module Procedure    Inline_Strings_1D
    Module Procedure    Inline_Reals_1D
    Module Procedure    Inline_Integers_1D
  End Interface

  Interface             Add_Prefix
    Module Procedure    Add_Prefix_0d
    Module Procedure    Add_Prefix_1d
  End Interface

  Interface             Add_Suffix
    Module Procedure    Add_Suffix_0d
    Module Procedure    Add_Suffix_1d
  End Interface

  Interface             Remove_Character
    Module Procedure    Remove_Character
    Module Procedure    Remove_Characters
  End Interface

  integer                       ::      LogUnit         =       Output_Unit                                     !< Unit used for logs
  integer       ,parameter      ::      i_ASCII_VertTab =       9
  integer       ,parameter      ::      i_ASCII_Space   =       32
  integer       ,parameter      ::      ierr_NO_ERROR   =       0

  contains

! This procedure sets the log unit used by all procedures in this module.
! The default log unit is the intrinsic Output_Unit.
! Note that no checking is done to ensure that the input Unit is actually connected to a opened file.
Subroutine Set_LogUnit( Unit )
  integer                                               ,intent(in)     ::      Unit
  LogUnit       =       Unit
End Subroutine

! **************************************************************************************************************
! **************************************************************************************************************
!       PROCEDURES FOR CHECKING WHETHER A VARIABLE IS PRESENT IN A LIST OF VARIABLES
! **************************************************************************************************************
! **************************************************************************************************************

Pure Function Presence_Character_0D( Var, List_Var, CaseSensitive ) result(Is_Present)
  character(*)                                          ,intent(in)     ::      Var                             !< Variable to be checked for presence in the list of variable
  character(*)          ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variable used for presence checking
  logical               ,optional                       ,intent(in)     ::      CaseSensitive                  !< Indicator whether the search should be case sensitive
  logical                                                               ::      Is_Present                      !< Presence indicator of the input variable in the list of variables
  logical                                                               ::      CaseSensitive_
  integer                                                               ::      i                            ! Index of elements in the list of variable
  character(:)  ,allocatable                                            ::      Str1, Str2

  CaseSensitive_ =       .true.
  if ( present(CaseSensitive) ) CaseSensitive_ = CaseSensitive

  Str1       =       trim(Var)
  if ( .Not.CaseSensitive_ ) Str1 = UpperCase( Str1 )

  Is_Present    =       .false.                                                                                 ! Initialization of the variable presence indicator to false
  do i = 1,size(List_Var)                                                                                       ! Loop on all variables elements
    Str2        =       trim(List_Var(i))
    if ( .Not.CaseSensitive_ ) Str2 = UpperCase( Str2 )
    if ( Str1 /= Str2 ) cycle                                                                                   ! If current variable from the list is different from the searched variable, then going to the next element
    Is_Present  =       .true.                                                                                  ! Setting presence indicator to true
    return                                                                                                      ! Exiting the procedure if no counting of occurence is required
  end do                                                                                                        ! End do loop on Vars

End Function

Pure Subroutine Exchange_Element_In_Array( Array, Old_Element, New_Element )
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      Array                           !<
  character(*)                                          ,intent(in)     ::      Old_Element                     !<
  character(*)                                          ,intent(in)     ::      New_Element                     !<
  integer                                                               ::      Element_Index, i
  integer                                                               ::      Length
  character(:)  ,dimension(:)   ,allocatable                            ::      Array_tmp
  Element_Index         =       FindPos_0D( Old_Element, Array )                                                   ! Finding the position of the element to be replaced in the array
  if ( Element_Index == 0 ) return                                                                              ! Exiting the procedure if the element has not beend found
  Length                =       max( len(Array), len(New_Element) )                                                     ! Getting the new length of the array
  allocate( character(Length) :: Array_tmp(size(Array)) )
  do i = 1,size(Array)
    Array_tmp(i)        =       Array(i)
    if ( i == Element_Index ) Array_tmp(i) = New_Element
  end do
  call move_alloc( Array_tmp, Array )
End Subroutine

!  Pure Function Presence_Character_0D_Deferred_Length( Var, List_Var ) result(Is_Present)
!   character(*)                                          ,intent(in)     ::      Var                             !< Variable to be checked for presence in the list of variable
!   character(:)  ,dimension(:)   ,allocatable            ,intent(in)     ::      List_Var                        !< List of variable used for presence checking
!   logical                                                               ::      Is_Present                      !< Presence indicator of the input variable in the list of variables
!   integer                                                               ::      iStr                            ! Index of elements in the list of variable
!   Is_Present    =       .false.                                                                                 ! Initialization of the variable presence indicator to false
!   do iStr = 1,size(List_Var)                                                                                    ! Loop on all variables elements
!     if ( trim(Var) /= trim(List_Var(iStr)) ) cycle                                                              ! If current variable from the list is different from the searched variable, then going to the next element
!     Is_Present  =       .true.                                                                                  ! Setting presence indicator to true
!     return                                                                                                      ! Exiting the procedure if no counting of occurence is required
!   end do                                                                                                        ! End do loop on Vars
! End Function

Pure Function Presence_Integer_0D( Var, List_Var ) result(Is_Present)
  integer                                               ,intent(in)     ::      Var                             !< Variable to be checked for presence in the list of variable
  integer               ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variable used for presence checking
  logical                                                               ::      Is_Present                      !< Presence indicator of the variable in the list of variables
  integer                                                               ::      iVar                            ! Index of elements in the list of variables
  Is_Present    =       .false.                                                                                 ! Initialization of the variable presence indicator to false
  do iVar = 1,size(List_Var)                                                                                    ! Loop on all elements
    if ( Var /= List_Var(iVar) ) cycle                                                                          ! If current variable from the list is different from the searched variable, then going to the next element
    Is_Present  =       .true.                                                                                  ! Setting presence indicator to true
    return                                                                                                      ! Exiting the procedure if no counting of occurence is required
  end do                                                                                                        ! End do loop on strings
End Function

Pure Function Presence_Real4_0D( Var, List_Var ) result(Is_Present)
  real(4)                                               ,intent(in)     ::      Var                             !< Variable to be checked for presence in the list of variable
  real(4)               ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variable used for presence checking
  logical                                                               ::      Is_Present                      !< Presence indicator of the variable in the list of variables
  integer                                                               ::      iVar                            ! Index of elements in the list of variables
  Is_Present    =       .false.                                                                                 ! Initialization of the variable presence indicator to false
  do iVar = 1,size(List_Var)                                                                                    ! Loop on all elements
    if ( Var /= List_Var(iVar) ) cycle                                                                          ! If current variable from the list is different from the searched variable, then going to the next element
    Is_Present  =       .true.                                                                                  ! Setting presence indicator to true
    return                                                                                                      ! Exiting the procedure if no counting of occurence is required
  end do                                                                                                        ! End do loop on strings
End Function

Pure Function Presence_Real8_0D( Var, List_Var ) result(Is_Present)
  real(8)                                             ,intent(in)     ::      Var                             !< Variable to be checked for presence in the list of variable
  real(8)               ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variable used for presence checking
  logical                                                               ::      Is_Present                      !< Presence indicator of the variable in the list of variables
  integer                                                               ::      iVar                            ! Index of elements in the list of variables
  Is_Present    =       .false.                                                                                 ! Initialization of the variable presence indicator to false
  do iVar = 1,size(List_Var)                                                                                    ! Loop on all elements
    if ( Var /= List_Var(iVar) ) cycle                                                                          ! If current variable from the list is different from the searched variable, then going to the next element
    Is_Present  =       .true.                                                                                  ! Setting presence indicator to true
    return                                                                                                      ! Exiting the procedure if no counting of occurence is required
  end do                                                                                                        ! End do loop on strings
End Function



Pure Function Presence_Character_1D( Var, List_Var, CaseSensitive ) result(Is_Present)
  character(*)          ,dimension(:)                   ,intent(in)     ::      Var                             !< Array of variables to be checked for presence in the list of variable
  character(*)          ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variable used for presence checking
  logical               ,optional                       ,intent(in)     ::      CaseSensitive                  !< Indicator whether the search should be case sensitive
  logical               ,dimension( size(Var) )                         ::      Is_Present                      !< Array of presence indicator of the variables in the list of variables
  integer                                                               ::      i                               ! Index of element of the array of variables
  forall(i=1:size(Var)) Is_Present(i) = Presence( Var(i), List_Var, CaseSensitive )                            ! Setting the presence indicators for all variables to be checked
End Function

Pure Function Presence_Integer_1D( Var, List_Var ) result(Is_Present)
  integer               ,dimension(:)                   ,intent(in)     ::      Var                             !< Array of variables to be checked for presence in the list of variable
  integer               ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variable used for presence checking
  logical               ,dimension( size(Var) )                         ::      Is_Present                      !< Array of presence indicator of the variables in the list of variables
  integer                                                               ::      i                               ! Index of element of the array of variables
  forall(i=1:size(Var)) Is_Present(i) = Presence( Var(i), List_Var )                                            ! Setting the presence indicators for all variables to be checked
End Function

Pure Function Presence_Real4_1D( Var, List_Var ) result(Is_Present)
  real(4)               ,dimension(:)                   ,intent(in)     ::      Var                             !< Array of variables to be checked for presence in the list of variable
  real(4)               ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variable used for presence checking
  logical               ,dimension( size(Var) )                         ::      Is_Present                      !< Array of presence indicator of the variables in the list of variables
  integer                                                               ::      i                               ! Index of element of the array of variables
  forall(i=1:size(Var)) Is_Present(i) = Presence( Var(i), List_Var )                                            ! Setting the presence indicators for all variables to be checked
End Function

Pure Function Presence_Real8_1D( Var, List_Var ) result(Is_Present)
  real(8)               ,dimension(:)                   ,intent(in)     ::      Var                             !< Array of variables to be checked for presence in the list of variable
  real(8)               ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variable used for presence checking
  logical               ,dimension( size(Var) )                         ::      Is_Present                      !< Array of presence indicator of the variables in the list of variables
  integer                                                               ::      i                               ! Index of element of the array of variables
  forall(i=1:size(Var)) Is_Present(i) = Presence( Var(i), List_Var )                                            ! Setting the presence indicators for all variables to be checked
End Function


! **************************************************************************************************************
! **************************************************************************************************************
!       PROCEDURES FOR COUNTING THE NUMBER OF TIME A VARIABLE IS PRESENT IN A LIST OF VARIABLES
! **************************************************************************************************************
! **************************************************************************************************************

Pure Function Count_Presence_Character_0D( Var, List_Var, Trimed, CaseSensitive ) result(NCounts)
  character(*)                                          ,intent(in)     ::      Var                             !< Variable whose presence in the list of variables is to be counted
  character(*)          ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variable used for counting
  logical               ,optional                       ,intent(in)     ::      Trimed                          !< Indicator whether
  logical               ,optional                       ,intent(in)     ::      CaseSensitive                  !< Indicator whether the search should be case sensitive
  integer                                                               ::      NCounts                         !< Number of occurence of the input variable in the list of variables
  integer                                                               ::      iVar                            ! Index of elements in the list of variables
  NCounts       =       0                                                                                       ! Initialization of the count
  do iVar = 1,size(List_Var)                                                                                    ! Loop on all elements in the list of variables
!     if ( trim(Var) == trim(List_Var(iVar)) ) NCounts = NCounts + 1                                              ! If current variable from the list is identical to the searched variable, then incrementing of the counts
    if ( Equal( Var, List_Var(iVar), Trimed=Trimed, CaseSensitive=CaseSensitive ) ) NCounts = NCounts + 1     ! If current variable from the list is identical to the searched variable, then incrementing of the counts
  end do                                                                                                        ! End do loop on elements
End Function

Pure Function Count_Presence_Integer_0D( Var, List_Var ) result(NCounts)
  integer                                               ,intent(in)     ::      Var                             !< Variable whose presence in the list of variables is to be counted
  integer               ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variable used for counting
  integer                                                               ::      NCounts                         !< Number of occurence of the input variable in the list of variables
  integer                                                               ::      iVar                            ! Index of elements in the list of variables
  NCounts       =       0                                                                                       ! Initialization of the count
  do iVar = 1,size(List_Var)                                                                                    ! Loop on all elements in the list of variables
    if ( Var == List_Var(iVar) ) NCounts = NCounts + 1                                                          ! If current variable from the list is identical to the searched variable, then incrementing of the counts
  end do                                                                                                        ! End do loop on elements
End Function

Pure Function Count_Presence_Real4_0D( Var, List_Var ) result(NCounts)
  real(4)                                               ,intent(in)     ::      Var                             !< Variable whose presence in the list of variables is to be counted
  real(4)               ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variable used for counting
  integer                                                               ::      NCounts                         !< Number of occurence of the input variable in the list of variables
  integer                                                               ::      iVar                            ! Index of elements in the list of variables
  NCounts       =       0                                                                                       ! Initialization of the count
  do iVar = 1,size(List_Var)                                                                                    ! Loop on all elements in the list of variables
    if ( Var == List_Var(iVar) ) NCounts = NCounts + 1                                                          ! If current variable from the list is identical to the searched variable, then incrementing of the counts
  end do                                                                                                        ! End do loop on elements
End Function

Pure Function Count_Presence_Real8_0D( Var, List_Var ) result(NCounts)
  real(8)                                               ,intent(in)     ::      Var                             !< Variable whose presence in the list of variables is to be counted
  real(8)               ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variable used for counting
  integer                                                               ::      NCounts                         !< Number of occurence of the input variable in the list of variables
  integer                                                               ::      iVar                            ! Index of elements in the list of variables
  NCounts       =       0                                                                                       ! Initialization of the count
  do iVar = 1,size(List_Var)                                                                                    ! Loop on all elements in the list of variables
    if ( Var == List_Var(iVar) ) NCounts = NCounts + 1                                                          ! If current variable from the list is identical to the searched variable, then incrementing of the counts
  end do                                                                                                        ! End do loop on elements
End Function

Pure Function Count_Presence_Character_1D( Var, List_Var, Trimed, CaseSensitive ) result(NCounts)
  character(*)          ,dimension(:)                   ,intent(in)     ::      Var                             !< Array of variables whose presence in the list of variable is to be counted
  character(*)          ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variables used for counting
  logical               ,optional                       ,intent(in)     ::      Trimed                          !< Indicator whether
  logical               ,optional                       ,intent(in)     ::      CaseSensitive                  !< Indicator whether the search should be case sensitive
  integer               ,dimension( size(Var) )                         ::      NCounts                         !< Number of occurence of each element of the array of variable in the list of variables
  integer                                                               ::      i                               ! Index of element in the array of variable
  forall(i=1:size(Var)) NCounts(i) = Count_Presence( Var(i), List_Var, Trimed, CaseSensitive )                 ! Counting the occurence of current element in the list of variable
End Function

Pure Function Count_Presence_Integer_1D( Var, List_Var ) result(NCounts)
  integer               ,dimension(:)                   ,intent(in)     ::      Var                             !< Array of variables whose presence in the list of variable is to be counted
  integer               ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variables used for counting
  integer               ,dimension( size(Var) )                         ::      NCounts                         !< Number of occurence of each element of the array of variable in the list of variables
  integer                                                               ::      i                               ! Index of element in the array of variable
  forall(i=1:size(Var)) NCounts(i) = Count_Presence( Var(i), List_Var )                                         ! Counting the occurence of current element in the list of variable
End Function

Pure Function Count_Presence_Real4_1D( Var, List_Var ) result(NCounts)
  real(4)               ,dimension(:)                   ,intent(in)     ::      Var                             !< Array of variables whose presence in the list of variable is to be counted
  real(4)               ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variables used for counting
  integer               ,dimension( size(Var) )                         ::      NCounts                         !< Number of occurence of each element of the array of variable in the list of variables
  integer                                                               ::      i                               ! Index of element in the array of variable
  forall(i=1:size(Var)) NCounts(i) = Count_Presence( Var(i), List_Var )                                         ! Counting the occurence of current element in the list of variable
End Function

Pure Function Count_Presence_Real8_1D( Var, List_Var ) result(NCounts)
  real(8)               ,dimension(:)                   ,intent(in)     ::      Var                             !< Array of variables whose presence in the list of variable is to be counted
  real(8)               ,dimension(:)                   ,intent(in)     ::      List_Var                        !< List of variables used for counting
  integer               ,dimension( size(Var) )                         ::      NCounts                         !< Number of occurence of each element of the array of variable in the list of variables
  integer                                                               ::      i                               ! Index of element in the array of variable
  forall(i=1:size(Var)) NCounts(i) = Count_Presence( Var(i), List_Var )                                         ! Counting the occurence of current element in the list of variable
End Function









Pure Function FindPos_0D( String, List_String, Back, CaseSensitive ) result(iLoc)

  character(*)                                          ,intent(in)     ::      String                          !< String to be checked for presence
  character(*)          ,dimension(:)                   ,intent(in)     ::      List_String                     !< List of string used for presence checking
  logical               ,optional                       ,intent(in)     ::      Back                            !< Indicator whether the position index is counted from the end of the list of strings
  logical               ,optional                       ,intent(in)     ::      CaseSensitive                  !< Indicator whether the search should be case sensitive
  integer                                                               ::      iLoc                            !< Index of the position of the input string in the list of strings

  integer                                                               ::      iStr                            ! Index of string of the input string in the list of strings
  integer                                                               ::      iStart                          ! Starting index used to found the string position
  integer                                                               ::      iStop                           ! Stopping index used to found the string position
  integer                                                               ::      iStep                           ! Step in index incrementation used to found the string position
  logical                                                               ::      CaseSensitive_
  character(:)  ,allocatable                                            ::      Target_String
  character(:)  ,allocatable                                            ::      Current_String

  CaseSensitive_  =       .true.
  if ( present(CaseSensitive) ) CaseSensitive_ = CaseSensitive
  iLoc          =       0                                                                                       ! Initialization of the position index
  iStart        =       1                                                                                       ! Setting the starting index to unity
  iStop         =       size(List_String)                                                                       ! Setting the stopping index to the dimension of the list of strings
  iStep         =       1                                                                                       ! Setting index step to one

  if ( present(Back) ) then                                                                                     ! If the "Back" optional argument is present
    if (Back) then                                                                                              ! If the "Back" argument is true
      iStart    =        size(List_String)                                                                      ! Setting the starting index to the dimension of the list of strings
      iStop     =        1                                                                                      ! Setting the stopping index to unity
      iStep     =       -1                                                                                      ! Setting index step to -1 to go backward
    end if                                                                                                      ! End if case on the argument value
  end if                                                                                                        ! End if case on optional argument presence

  Target_String         =       trim(String)
  if ( .Not. CaseSensitive_ ) Target_String = UpperCase(Target_String)
  do iStr = iStart,iStop,iStep                                                                                  ! Loop on all strings
    Current_String      =       trim(List_String(iStr))
    if ( .Not. CaseSensitive_ ) Current_String = UpperCase(Current_String)
    if ( Target_String /= Current_String ) cycle                                                                ! If current string from the list is different from the searched string, then going to the next element
    iLoc                =       iStr                                                                            ! Setting the index if the string
    return                                                                                                      ! Exiting the procedure if no counting of occurence is required
  end do                                                                                                        ! End do loop on strings

End Function

Pure Function FindPos_1D( String, List_String, Back, CaseSensitive ) result(iLoc)

  character(*)          ,dimension(:)                   ,intent(in)     ::      String                          !< String to be checked for presence
  character(*)          ,dimension(:)                   ,intent(in)     ::      List_String                     !< List of string used for presence checking
  logical               ,optional                       ,intent(in)     ::      Back                            !< Indicator than the position index is counted from the end of the list of strings
  logical               ,optional                       ,intent(in)     ::      CaseSensitive                  !< Indicator whether the search should be case sensitive
  integer               ,dimension( size(String) )                      ::      iLoc                            !< Index of positions of input strings in the list of strings

  integer                                                               ::      iStr                            ! Index of string of the input string in the list of strings

  do iStr = 1,size(String)                                                                                      ! Loop on all strings to be cheched for presence
     iLoc(iStr) =       FindPos( String(iStr), List_String, Back, CaseSensitive )                              ! Finding the position of the current string in the list of strings
  end do                                                                                                        ! End loop on strings to be cheched for presence

End Function


!<==============================================================================================================
!> @brief       Orders string array elements according to a reference stroig array
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!> @todo        String with different should be treated (see KinSpe_Class::Check_SpeCat for possible implementation)
!<==============================================================================================================
!> @details
!! This procedure orders string array elements according to a reference strong array. \n
!! In the following, the string to be ordered is denoted "object string" whereas the string from which the
!! ordering is performed is called the "reference string". \n
!! It is assumed that lower dimension of both the object and reference strings start at 1. \n
!! If the two strings have different dimension, then the output optional error indicator is set to one and
!! the procedure is exited whitout any reordering. \n
!! An insertion sorting method is used to sort the array elements.
!<==============================================================================================================
Subroutine Reorder_Character( Array, ValRef, IdxRef, Mapping, i_Debug )

  character(*)                  ,dimension(:)           ,intent(inout)  ::      Array                           !< Input array on which the ordering is performed
  character(*)  ,optional       ,dimension(:)           ,intent(in)     ::      ValRef                          !< Reference value array from which the ordering is performed
  integer       ,optional       ,dimension(:)           ,intent(in)     ::      IdxRef                          !< Reference index array from which the ordering is performed
  integer       ,optional       ,dimension(:)           ,intent(out)    ::      Mapping                         !< Index mapping from old to new order
  logical       ,optional                               ,intent(in)     ::      i_Debug                         !< Debugging indicator

  logical                                                               ::      i_Debug_Loc                     ! Local debugging indicator
  logical                                                               ::      i_Value_Reordering              ! Indicator of value-based reordering
  logical                                                               ::      i_Index_Reordering              ! Indicator of index-based reordering
  integer                                                               ::      iNew, iOld, jOld                ! New/Old element index
  integer                       ,dimension( size(Array) )               ::      iOld_to_iNew                    ! Index correspondance from old to new index
  integer                       ,dimension( size(Array) )               ::      Mapping_Loc                     ! Local index mapping from old to new order
  character(:)  ,dimension(:)   ,allocatable                            ::      Array_Loc                       ! Local copy of the Input-array
  character(:)  ,dimension(:)   ,allocatable                            ::      Array_Ref                       ! Local reference value array from which the ordering is performed (called Reference-array)
  character(*)  ,parameter                                              ::      ProcName='Reorder_Character'    ! Procedure name

  i_Debug_Loc = .false.; if ( present(i_Debug) ) i_Debug_Loc = i_Debug                                                                 ! Setting local debugging indicator
  if (i_Debug_Loc) write(LogUnit,"(2x,'[Reorder_Character]: Entering')")                                        ! Debugging

! ==============================================================================================================
!    TREATING OPTIONAL INPUT ARGUMENTS
! ==============================================================================================================
  i_Value_Reordering      =       .false.                                                                       ! Initialization of the value-based reordering indicator to false
  i_Index_Reordering      =       .false.                                                                       ! Initialization of the index-based reordering indicator to false
  if ( present(ValRef) ) i_Value_Reordering = .true.                                                            ! If the reference-value array is present in the list of arguments, then setting the value-based reordering indicator to true
  if ( present(IdxRef) ) i_Index_Reordering = .true.                                                            ! If the reference-index array is present in the list of arguments, then setting the index-based reordering indicator to true
  if (i_Debug_Loc) write(LogUnit,"(2x,'[Reorder_Character]: i_Value_Reordering = ',l3)") i_Value_Reordering
  if (i_Debug_Loc) write(LogUnit,"(2x,'[Reorder_Character]: i_Index_Reordering = ',l3)") i_Index_Reordering

! ==============================================================================================================
!    CHECKING DATA CONSISTENCY
! ==============================================================================================================
  if ( size(Array) == 0 ) then
    if (i_Debug_Loc) write(LogUnit,"(2x,'[Reorder_Character]: Input argument is a zero-dimension array => Exiting')")
    return
  end if

  if ( i_Value_Reordering .and. i_Index_Reordering ) &                                                          ! If both index-based and value-based reordering are considered, then ...
  call Error_Message( ErrMsg="Both the index-based and value-based reordering options cannot be specified together", ProcName=ProcName ) ! ... printing an error message and stopping the code

  if ( i_Index_Reordering ) then                                                                                ! If index-based reordering is considered, then checking that the array and index variable has the same dimension
    if ( size(Array) /= size(IdxRef) ) &                                                                        ! If size of the reference index variable and the array variables which need to be reordered does not match, then ...
    call Error_Message( ErrMsg="The size of the reference index variable and the array variables which need to be reordered does not match", ProcName=ProcName ) ! ... printing an error message and stopping the code
  end if                                                                                                        ! End if case on index-based reordering


! ==============================================================================================================
!    IF VALUE-BASED ORDERING, SETTING THE ARRAY OF ORDERED Input-array INDEX
! ==============================================================================================================
  if ( i_Value_Reordering ) then                                                                                ! If value-based reordering
    if (i_Debug_Loc) write(LogUnit,"(2x,'[Reorder_Character]: Value-based ordering')")                          ! Debugging
#ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
    allocate( character(len(ValRef)) :: Array_Ref(size(ValRef)) )
    Array_Ref = ValRef
#else
    allocate( Array_Ref, source=ValRef )                                                                        ! Allocating the local Reference-array
#endif
    do iOld = 1,size(Array)                                                                                     ! Loop on Input-array's elements
      do iNew = 1,size(Array_Ref)                                                                               ! Loop on Reference-array's elements
        if ( Array(iOld) == Array_Ref(iNew) ) iOld_to_iNew(iOld) = iNew                                         ! Setting index correspondance from Input-array to Reference-array
      end do                                                                                                    ! End loop on Reference-array's elements
    end do                                                                                                      ! End loop on Input-array's elements
    Mapping_Loc =       1                                                                                       ! Initialization of ordered Input-array index
    do iOld = 1,size(Array)                                                                                     ! Loop on Input-array's elements
      do jOld = 1,size(Array)                                                                                   ! Loop on Input-array's elements
        if ( (iOld/=jOld) .and. (iOld_to_iNew(iOld)>iOld_to_iNew(jOld)) ) Mapping_Loc(iOld) = Mapping_Loc(iOld) + 1   ! If the the two Input-array elements are different and if the first index if greater than the second, then Setting the index of the ordered Input-array
      end do                                                                                                      ! End loop on Input-array's elements
    end do                                                                                                      ! End loop on Input-array's elements
  end if                                                                                                        ! End of value-based reordering

! ==============================================================================================================
!    IF INDEX-BASED ORDERING, COPYING THE ARRAY OF ORDERED Input-array INDEX FROM Input
! ==============================================================================================================
  if ( i_Index_Reordering ) then                                                                                ! If index-based reordering
    if (i_Debug_Loc) write(LogUnit,"(2x,'[Reorder_Character]: Index-based ordering')")                          ! Debugging
    Mapping_Loc =       IdxRef                                                                                  ! Copying the array of ordered Input-array index from Input
  end if                                                                                                        ! End of index-based reordering

! ==============================================================================================================
!    SETTING THE ORDERED VALUES INTO THE OUTPUT ARRAY
! ==============================================================================================================
#ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
  allocate( character(len(Array)) :: Array_Loc(size(Array)) )
  Array_Loc = Array
#else
  allocate( Array_Loc, source=Array )                                                                           ! Creating a temporary copy of the Input-array
#endif
  do iOld = 1,size(Array)                                                                                       ! Loop on Input-array's elements
    Array( Mapping_Loc(iOld) )  =       Array_Loc(iOld)                                                         ! Setting the ordered Input-array value
  end do                                                                                                        ! End loop on Input-array's elements

! ==============================================================================================================
!    SETTING OPTIONAL OUTPUT ARGUMENTS
! ==============================================================================================================
  if ( present(Mapping) ) Mapping = Mapping_Loc                                                                 ! Storing the array of ordered Input-array index in the output variable

  if (i_Debug_Loc) then
    do iOld = 1,size(Array)
      write(LogUnit,"(2x,'[Reorder_Character]: iOld=',i3,3x,'New = Mapping_Loc(iOld)=',i3,3x,'Array(iOld)=',a,3x,'Array(iNew)=',a)")    &
      iOld, Mapping_Loc(iOld), Array_Loc(iOld), Array(iOld)
    end do
    write(LogUnit,"(2x,'[Reorder_Character]: Exiting')")
  end if

End Subroutine

Subroutine Reorder_Real_8( Inp, ValRef, IdxRef, Mapping )
  real(8)                       ,dimension(:)           ,intent(inout)  ::      Inp                             !< Input array on which the ordering is performed (called Inp-array)
  real(8)       ,optional       ,dimension(:)           ,intent(in)     ::      ValRef                          !< Reference value array from which the ordering is performed
  integer       ,optional       ,dimension(:)           ,intent(in)     ::      IdxRef                          !< Reference index array from which the ordering is performed
  integer       ,optional       ,dimension(:)           ,intent(out)    ::      Mapping                             !< Output array giving the reordering index

  logical                                                               ::      i_Value_Reordering                        ! Indicator of value-based reordering
  logical                                                               ::      i_Index_Reordering                        ! Indicator of index-based reordering
  integer                                                               ::      iInp, jInp                      ! Index element in the Inp-array
  real(8)                       ,dimension(size(Inp))                   ::      Array_Loc                          ! Local copy of the Inp-array
  integer                       ,dimension(size(Inp))                   ::      Mapping_Loc                     ! Array of ordered Inp-array index


  if ( size(Inp) == 0 ) then
    write(LogUnit,"(2x,'[Reorder_Real_8]: Input argument is a zero-dimension array')")
    write(LogUnit,"(2x,'[Reorder_Real_8]: Exiting')")
    return
  end if


! ==============================================================================================================
!    TREATING OPTIONAL INPUT ARGUMENTS
! ==============================================================================================================
  i_Value_Reordering      =       .false.
  i_Index_Reordering      =       .false.
  if    ( present(ValRef) )     i_Value_Reordering        =       .true.                                                  ! If the reference value array is provided in input, then setting the corresponding indicator
  if    ( present(IdxRef) )     i_Index_Reordering        =       .true.                                                  ! If the reference index array is provided in input, then setting the corresponding indicator
  if    (i_Value_Reordering.and.i_Index_Reordering) then
    write(*,"(2x,'[Reorder_Real_8]: ERROR: Both value-based and index-based reordering is forbidden')")
    stop
  end if

! ==============================================================================================================
!    IF VALUE-BASED ORDERING, SETTING THE ARRAY OF ORDERED INP-ARRAY INDEX
! ==============================================================================================================
  if    ( i_Value_Reordering )    then                                                                                    ! If value-based reordering
    write(*,"(2x,'[Reorder_Real_8]: ERROR: Value-based ordering not implemented for real numbers')")
    stop
  end if                                                                                                        ! End of value-based reordering

! ==============================================================================================================
!    IF INDEX-BASED ORDERING, COPYING THE ARRAY OF ORDERED INP-ARRAY INDEX FROM INPUT
! ==============================================================================================================
  if    ( i_Index_Reordering )    then                                                                                    ! If index-based reordering
    if  (size(Inp) /= size(IdxRef))     then
      write(*,"(2x,'[Reorder_Real_8]: ERROR: size of Inp and IdxRef are different')")
      stop
    end if
    Mapping_Loc =       IdxRef                                                                                  ! Copying the array of ordered inp-array index from input
  end if                                                                                                        ! End of index-based reordering

! ==============================================================================================================
!    SETTING THE ORDERED VALUES INTO THE OUTPUT ARRAY
! ==============================================================================================================
  Array_Loc        =       Inp                                                                                     ! Creating a temporary copy of the Inp-array
  do iInp = 1,size(Inp)                                                                                         ! Loop on Inp-array's elements
    jInp        =       Mapping_Loc(iInp)                                                                       ! Getting the index of the ordered Inp-array element
    Inp(jInp)   =       Array_Loc(iInp)                                                                            ! Setting the ordered Inp-array value
  end do                                                                                                        ! End loop on Inp-array's elements

! ==============================================================================================================
!    SETTING OPTIONAL OUTPUT ARGUMENTS
! ==============================================================================================================
  if    (present(Mapping))          Mapping     =       Mapping_Loc                                                     ! Storing the array of ordered inp-array index in the output variable

End Subroutine

Subroutine Reorder_Real_4( Inp, ValRef, IdxRef, Mapping )
  real(4)                       ,dimension(:)           ,intent(inout)  ::      Inp                             !< Input array on which the ordering is performed (called Inp-array)
  real(4)       ,optional       ,dimension(:)           ,intent(in)     ::      ValRef                          !< Reference value array from which the ordering is performed
  integer       ,optional       ,dimension(:)           ,intent(in)     ::      IdxRef                          !< Reference index array from which the ordering is performed
  integer       ,optional       ,dimension(:)           ,intent(out)    ::      Mapping                             !< Output array giving the reordering index

  logical                                                               ::      i_Value_Reordering                        ! Indicator of value-based reordering
  logical                                                               ::      i_Index_Reordering                        ! Indicator of index-based reordering
  integer                                                               ::      iInp, jInp                      ! Index element in the Inp-array
  real(4)                       ,dimension(size(Inp))                   ::      Array_Loc                          ! Local copy of the Inp-array
  integer                       ,dimension(size(Inp))                   ::      Mapping_Loc                     ! Array of ordered Inp-array index


  if ( size(Inp) == 0 ) then
    write(LogUnit,"(2x,'[Reorder_Real_4]: Input argument is a zero-dimension array')")
    write(LogUnit,"(2x,'[Reorder_Real_4]: Exiting')")
    return
  end if


! ==============================================================================================================
!    TREATING OPTIONAL INPUT ARGUMENTS
! ==============================================================================================================
  i_Value_Reordering      =       .false.
  i_Index_Reordering      =       .false.
  if    ( present(ValRef) )     i_Value_Reordering        =       .true.                                                  ! If the reference value array is provided in input, then setting the corresponding indicator
  if    ( present(IdxRef) )     i_Index_Reordering        =       .true.                                                  ! If the reference index array is provided in input, then setting the corresponding indicator
  if    (i_Value_Reordering.and.i_Index_Reordering) then
    write(*,"(2x,'[Reorder_Real_4]: ERROR: Both value-based and index-based reordering is forbidden')")
    stop
  end if

! ==============================================================================================================
!    IF VALUE-BASED ORDERING, SETTING THE ARRAY OF ORDERED INP-ARRAY INDEX
! ==============================================================================================================
  if    ( i_Value_Reordering )    then                                                                                    ! If value-based reordering
    write(*,"(2x,'[Reorder_Real_4]: ERROR: Value-based ordering not implemented for real numbers')")
    stop
  end if                                                                                                        ! End of value-based reordering

! ==============================================================================================================
!    IF INDEX-BASED ORDERING, COPYING THE ARRAY OF ORDERED INP-ARRAY INDEX FROM INPUT
! ==============================================================================================================
  if    ( i_Index_Reordering )    then                                                                                    ! If index-based reordering
    if  (size(Inp) /= size(IdxRef))     then
      write(*,"(2x,'[Reorder_Real_4]: ERROR: size of Inp and IdxRef are different')")
      stop
    end if
    Mapping_Loc =       IdxRef                                                                                  ! Copying the array of ordered inp-array index from input
  end if                                                                                                        ! End of index-based reordering

! ==============================================================================================================
!    SETTING THE ORDERED VALUES INTO THE OUTPUT ARRAY
! ==============================================================================================================
  Array_Loc        =       Inp                                                                                     ! Creating a temporary copy of the Inp-array
  do iInp = 1,size(Inp)                                                                                         ! Loop on Inp-array's elements
    jInp        =       Mapping_Loc(iInp)                                                                       ! Getting the index of the ordered Inp-array element
    Inp(jInp)   =       Array_Loc(iInp)                                                                            ! Setting the ordered Inp-array value
  end do                                                                                                        ! End loop on Inp-array's elements

! ==============================================================================================================
!    SETTING OPTIONAL OUTPUT ARGUMENTS
! ==============================================================================================================
  if    (present(Mapping))          Mapping     =       Mapping_Loc                                                     ! Storing the array of ordered inp-array index in the output variable

End Subroutine

! Subroutine Reorder_Character( Str_Obj, Str_Ref, ierr )
! !   character(*)                        ,dimension(:)           ,intent(inout)  ::      Str_Obj                         !< Object string on which the ordering is performed
!   character(*)                        ,dimension(:)           ,intent(in)     ::      Str_Ref                         !< Reference string from which the ordering is performed
!   integer             ,optional                       ,intent(out)    ::      ierr                            !< Error indicator (on output ierr=0 if no error)
!   character(len(Str_Obj))     ,dimension(size(Str_Obj))               ::      Str_Tem                         ! Temorary string
!   integer                     ,dimension(size(Str_Obj))               ::      Ind_Obj                         ! Index order of components of the object string in the frame of the reference string
!   integer                                                             ::      iobj                            ! Index of components of the object string
!   integer                                                             ::      iNew                            ! Index of components of the reference string
!   integer                                                             ::      I_wrk                           ! Local index
!   character(len(Str_Obj))                                             ::      S_wrk                           ! Local string component
!   if  (size(Str_Obj) /= size(Str_Ref))        then
!     if        (present(ierr)) ierr    =       1
! !     return
!   end if
!   Str_Tem     =       Str_Obj                                                                                 ! Setting the temporary string to the object string values
!   Ind_Obj     =       0                                                                                       ! Initialization of the index order
!   do iobj = 1,size(Str_Obj)                                                                                   ! Loop on all object components
!   do iNew = 1,size(Str_Ref)                                                                                   ! Loop on all reference components
!     if        (Str_Obj(iobj) == Str_Ref(iNew))        Ind_Obj(iobj)   =       iNew                                    ! If the two components are identical then setting the index order
!   end do                                                                                                      ! End loop on reference components
!   end do                                                                                                      ! End loop on object components
!   do iobj = 2,size(Str_Obj)
!     I_wrk     =       Ind_Obj(iobj)
!     S_wrk     =       Str_Tem(iobj)
!     if        (I_wrk >= Ind_Obj(iobj-1))      cycle
!     Ind_Obj(iobj)     =           Ind_Obj(iobj-1)
!     Str_Tem(iobj)     =           Str_Tem(iobj-1)
!     do iNew = iobj-2,1,-1
!       if      (I_wrk >= Ind_Obj(iNew))        exit
!       Ind_Obj(iNew+1) =       Ind_Obj(iNew)
!       Str_Tem(iNew+1) =       Str_Tem(iNew)
!     end do
!     Ind_Obj(iNew+1)   =       I_wrk
!     Str_Tem(iNew+1)   =       S_wrk
!   end do
!   Str_Obj     =       Str_Tem                                                                                 ! Saving the temporary object into the final object
!   if  (present(ierr))         ierr    =       ierr_NO_ERROR
! End Subroutine

!<==============================================================================================================
!> @brief       Indicates whether of not a string corresponds to a digit
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure indicates whether of not a string corresponds to a digit. \n
!! The function's result a logical variable which has the following values:
!!  - true if the input string is a digit (0,1,...,9)
!!  - false otherwise
!<==============================================================================================================
Function is_digit(str)  result(res)
  character     ,intent(in)     ::      str                                                                     !< Input character string to be checked for didgits
  logical                       ::      res                                                                     !< Output indicator whether of not the input character is a number
  select case   (str)
  case  ('0':'9')
    res =       .true.
  case default
    res =       .false.
  end select
End Function

!<==============================================================================================================
!> @brief       Indicates whether of not a string corresponds to a letter
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure indicates whether of not a string corresponds to a letter, either lower or upper case. \n
!! The function's result a logical variable which has the following values:
!!  - true if the input string is a letter (a,...,z,A,...,Z)
!!  - false otherwise
!<==============================================================================================================
Pure Elemental Function Is_Letter(str) result(res)
  character     ,intent(in)     ::      str                                                                     !< Input character string to be checked for letters
  logical                       ::      res                                                                     !< Output indicator whether of not the input character is a letter
  select case   (str)
  case  ('A':'Z','a':'z')
    res =       .true.
  case default
    res =       .false.
  end select
End Function

!
!
! Subroutine Convert( From, To )
!   class(*)                                              ,intent(in)     ::      From
!   class(*)  ,allocatable                                ,intent(out)    ::      To
!   real(8)                                                               ::      RealNumber                      ! Real number corresponding to the input string
!   select type (From)
!   type is ()
!   type is ()
!   type is ()
!   end select
! End Function


!<==============================================================================================================
!> @brief       Converts a number string to an integer value
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure converts a number string to an integer value. \n
!! It is called using the generic interface \c Convert_To_Integer. \n
!<==============================================================================================================
Function String_To_Integer_0D( String ) result(Number)
  character(*)                                          ,intent(in)     ::      String                          !< Input character string to be converted
  integer                                                               ::      Number                          !< Integer number corresponding to the input string
  real(8)                                                               ::      RealNumber                      ! Real number corresponding to the input string
  RealNumber    =       Convert_To_Real( String )                                                               ! Converting the input string into a real number
  Number        =       nint( RealNumber )                                                                      ! Converting the real number into a integer number
End Function

Function String_To_Integer_1D( String ) result(Number)
  character(*)  ,dimension(:)                           ,intent(in)     ::      String                          !< Input character string to be converted
  integer       ,dimension( size(String) )                              ::      Number                          !< Integer number corresponding to the input string
  integer                                                               ::      i                               ! Index of array elements
  do i = 1,size(String)                                                                                         ! Loop on all array's elements
    Number(i)   =       Convert_To_Integer( String(i) )                                                         ! Converting the input string into an integer number
  end do                                                                                                        ! End loop on array's elements
End Function

!<==============================================================================================================
!> @brief       Converts a number string to an real value
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure converts a number string to an real value. \n
!! It is called using the generic interface \c Convert_To_Real. \n
!<==============================================================================================================
Pure Function String_To_Real8_0D( String ) result(Number)
  character(*)                                          ,intent(in)     ::      String                          !< Input character string to be converted
  real(8)                                                               ::      Number                          !< Real number corresponding to the input string
  read(String,*) Number                                                                                         ! Converting the input string into a real number
End Function

!  Pure Function String_To_Real4_0D( String ) result(Number)
!   character(*)                                          ,intent(in)     ::      String                          !< Input character string to be converted
!   real(4)                                                               ::      Number                          !< Real number corresponding to the input string
!   read(String,*) Number                                                                                         ! Converting the input string into a real number
! End Function

Function String_To_Real8_1D( String ) result(Number)
  character(*)  ,dimension(:)                           ,intent(in)     ::      String                          !< Input character string to be converted
  real(8)       ,dimension( size(String) )                              ::      Number                          !< Real number corresponding to the input string
  integer                                                               ::      i                               ! Index of array elements
  do i = 1,size(String)                                                                                         ! Loop on all array's elements
    read(String(i),*) Number(i)                                                                                 ! Converting the input string into a real number
  end do                                                                                                        ! End loop on array's elements
End Function

! Function String_To_Real4_1D( String ) result(Number)
!   character(*)  ,dimension(:)                           ,intent(in)     ::      String                          !< Input character string to be converted
!   real(4)       ,dimension( size(String) )                              ::      Number                          !< Real number corresponding to the input string
!   integer                                                               ::      i                               ! Index of array elements
!   do i = 1,size(String)                                                                                         ! Loop on all array's elements
!     read(String(i),*) Number(i)                                                                                 ! Converting the input string into a real number
!   end do                                                                                                        ! End loop on array's elements
! End Function

!<==============================================================================================================
!> @brief       Converts a logical string to an logical value
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure converts a logical string to an logical value. \n
!! It is called using the generic interface \c Convert_To_Logical. \n
!! The output result is:
!!  - true  if the logical string  is either "TRUE", "ON" or 1
!!  - false if the logical string  is either "FALSE", "OFF" or 0
!<==============================================================================================================
Function String_To_Logical_0D( String ) result(Var)
  character(*)                                          ,intent(in)     ::      String                          !< Input character string to be converted
  logical                                                               ::      Var                             !< Logical variable corresponding to the input string
  character(*)  ,dimension(6)   ,parameter              ::      Valid_TRUE =[ 'TRUE  ', 'ON    ', 'YES   ', '1     ', 'T     ', '.TRUE.' ]        ! Keyword indicating a true value
  Var   =       Presence( UpperCase(adjustl(trim(String))) ,Valid_TRUE)
End Function

!<==============================================================================================================
!> @brief       Extracts the number located at the LHS of a scalar string
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure extracts the number located at the LHS of a scalar string. \n
!! It is called using the generic interface \c Get_Numbers_AtLeftOf_Characters. \n
!<==============================================================================================================
Subroutine Get_Numbers_AtLeftOf_Characters_0D ( Input_String, Output_Number, default )
  character(*)                                          ,intent(in)     ::      Input_String                    !< String whose all numbers at the LHS of the first set of letters have to be extracted
  real(8)                                               ,intent(out)    ::      Output_Number                   !<
  character(*)  ,optional                               ,intent(in)     ::      default                         !<
  integer                                                               ::      i
  character(len(Input_String))                                          ::      StrNum
  i             =       0
  do
    i           =       i + 1
    if  (i>len_trim(Input_String))       exit
    if  (Is_Letter(Input_String(i:i)))   exit
  end do
  if ( i == 1 ) then
    if  (present(default))      then
      StrNum    =       default
    else
      StrNum    =       '1'
    end if
  else
    StrNum      =       Input_String(1:i-1)
  end if
  Output_Number       =       Convert_To_Real(StrNum)
End Subroutine

Subroutine Get_Numbers_AtLeftOf_Characters_1D ( Input_String, Output_Number, default )
  character(*)  ,dimension(:)                           ,intent(in)     ::      Input_String                    !< Array of string whose all numbers at the LHS of the first set of letters have to be extracted
  real(8)       ,dimension(:)   ,allocatable            ,intent(out)    ::      Output_Number                   !<
  character(*)  ,optional                               ,intent(in)     ::      default                         !<
  integer                                                               ::      i                               !
  allocate( Output_Number(size(Input_String)) )
  do i = 1,size(Input_String)
    call Get_Numbers_AtLeftOf_Characters( Input_String(i), Output_Number(i), default )
  end do
End Subroutine

!! This procedure extracts all characters from a string which are located at the right of the first set of letters.
!! If the input string starts with a letter, then the ouput and input string are identical
!! if the input string has only numbers, then the ouput is an empty string.
! String_NumbersLetters = "15N2"        =>      String_Letters = "N2"
! String_NumbersLetters = "N2"          =>      String_Letters = "N2"
! String_NumbersLetters = "154"         =>      String_Letters = ""
Pure Subroutine Get_Characters_AtRightOf_Numbers_1d( Input_String, Output_String )

  character(*)  ,dimension(:)                           ,intent(in)     ::      Input_String                    !< Array of string whose all letters of the RHS of the first set of numbers have to be extracted
  character(:)  ,dimension(:)   ,allocatable            ,intent(out)    ::      Output_String                   !<

  integer                                                               ::      i, k
  integer                                                               ::      Length                          ! Length of the output string
  integer                                                               ::      LenLoc                          ! Local length a element
  integer                                                               ::      NElements                       ! Number of elements in the input string
  character(:)  ,allocatable                                            ::      String
  character(:)  ,allocatable                                            ::      String_Letters
  character(:)  ,allocatable                                            ::      String_Numbers

  NElements             =       size(Input_String)                                                              ! Getting the number of elements of the input string
  Length                =       0                                                                               ! Initializing the length of the output string (Required because of iterative computation)
  do i = 1,NElements                                                                                            ! Loop on all the elements to be processed
    String              =       Input_String(i)                                                                 ! Copŷing current string in a new variable (Required in order to process each character one-by-one)
    LenLoc              =       len_trim(String)                                                                ! Getting the length of current element (without trailing blancs)
    do k = 1,LenLoc                                                                                             ! Loop on all character of current element in order to find the index of te first character corresponding to a letter
      if ( Is_Letter(String(k:k)) ) exit                                                                        ! If current character corresponds to a letter, then exiting
    end do                                                                                                      ! End loop on characters of current element
    if ( k == 1 ) then                                                                                          ! If the first character corresponds to a letter, then the full string corresponds to the string we are looking for
      String_Numbers    =       ""                                                                              ! Extracting the string of numbers
      String_Letters    =       trim(String)                                                                    ! Extracting the string of letters
    else                                                                                                        ! Splitting the string into the number and the letter part.
      String_Numbers    =       String(1:k-1)                                                                   ! Extracting the string of numbers
      String_Letters    =       String(k:LenLoc)                                                                ! Extracting the string of letters
    end if                                                                                                      ! End if case on the first lettr index
    Length              =       max( Length, len_trim(String_Letters) )                                         ! Getting the maximum length
  end do

  allocate( character(Length) :: Output_String(NElements) )                                                     ! Allocating the output string
  Output_String         =       ""                                                                              ! Initializing the output string values (Required if somme elements have only numbers, then the output string is an empty string)
  do i = 1,NElements                                                                                            ! Loop on all the elements to be processed
    String              =       Input_String(i)                                                                 ! Copŷing current string in a new variable (Required in order to process each character one-by-one)
    LenLoc              =       len_trim(String)                                                                ! Getting the length of current element (without trailing blancs)
    do k = 1,LenLoc                                                                                             ! Loop on all character of current element in order to find the index of te first character corresponding to a letter
      if ( Is_Letter(String(k:k)) ) exit                                                                        ! If current character corresponds to a letter, then exiting
    end do                                                                                                      ! End loop on characters of current element
    if ( k == 1 ) then                                                                                          ! If the first character corresponds to a letter, then the full string corresponds to the string we are looking for
      String_Numbers    =       ""                                                                              ! Extracting the string of numbers
      String_Letters    =       trim(String)                                                                    ! Extracting the string of letters
    else                                                                                                        ! Splitting the string into the number and the letter part.
      String_Numbers    =       String(1:k-1)                                                                   ! Extracting the string of numbers
      String_Letters    =       String(k:LenLoc)                                                                ! Extracting the string of letters
    end if                                                                                                      ! End if case on the first lettr index
    Output_String(i)    =       trim(String_Letters)                                                            ! Setting the output string
  end do

End Subroutine

Pure Subroutine Get_Characters_AtRightOf_Numbers_0d( Input_String, Output_String )

  character(*)                                          ,intent(in)     ::      Input_String                    !< Array of string whose all letters of the RHS of the first set of numbers have to be extracted
  character(:)  ,allocatable                            ,intent(out)    ::      Output_String                   !<

  integer                                                               ::      k
  integer                                                               ::      Length                          ! Length of the output string
  integer                                                               ::      LenLoc                          ! Local length a element

  character(:)  ,allocatable                                            ::      String
  character(:)  ,allocatable                                            ::      String_Letters
  character(:)  ,allocatable                                            ::      String_Numbers

    String              =       Input_String                                                                    ! Copying current string in a new variable (Required in order to process each character one-by-one)
    LenLoc              =       len_trim(String)                                                                ! Getting the length of current element (without trailing blancs)
    do k = 1,LenLoc                                                                                             ! Loop on all character of current element in order to find the index of te first character corresponding to a letter
      if ( Is_Letter(String(k:k)) ) exit                                                                        ! If current character corresponds to a letter, then exiting
    end do                                                                                                      ! End loop on characters of current element
    if ( k == 1 ) then                                                                                          ! If the first character corresponds to a letter, then the full string corresponds to the string we are looking for
      String_Numbers    =       ""                                                                              ! Extracting the string of numbers
      String_Letters    =       trim(String)                                                                    ! Extracting the string of letters
    else                                                                                                        ! Splitting the string into the number and the letter part.
      String_Numbers    =       String(1:k-1)                                                                   ! Extracting the string of numbers
      String_Letters    =       String(k:LenLoc)                                                                ! Extracting the string of letters
    end if                                                                                                      ! End if case on the first lettr index
    Length              =       len_trim(String_Letters)                                                        ! Getting the maximum length


  allocate( character(Length) :: Output_String )                                                     ! Allocating the output string
  Output_String         =       ""                                                                              ! Initializing the output string values (Required if somme elements have only numbers, then the output string is an empty string)

    String              =       Input_String                                                                 ! Copŷing current string in a new variable (Required in order to process each character one-by-one)
    LenLoc              =       len_trim(String)                                                                ! Getting the length of current element (without trailing blancs)
    do k = 1,LenLoc                                                                                             ! Loop on all character of current element in order to find the index of te first character corresponding to a letter
      if ( Is_Letter(String(k:k)) ) exit                                                                        ! If current character corresponds to a letter, then exiting
    end do                                                                                                      ! End loop on characters of current element
    if ( k == 1 ) then                                                                                          ! If the first character corresponds to a letter, then the full string corresponds to the string we are looking for
      String_Numbers    =       ""                                                                              ! Extracting the string of numbers
      String_Letters    =       trim(String)                                                                    ! Extracting the string of letters
    else                                                                                                        ! Splitting the string into the number and the letter part.
      String_Numbers    =       String(1:k-1)                                                                   ! Extracting the string of numbers
      String_Letters    =       String(k:LenLoc)                                                                ! Extracting the string of letters
    end if                                                                                                      ! End if case on the first lettr index
    Output_String       =       trim(String_Letters)                                                            ! Setting the output string


End Subroutine

!
! Subroutine GetStrRHS_1D ( Str, Str_RHS )
!   character(*)  ,dimension(:)                   ,intent(in)     ::      Str                                                     !<
! !   character(*)  ,dimension(:)                   ,intent(out)    ::      Str_RHS                                                 !<
!   character(:)  ,dimension(:)   ,allocatable    ,intent(out)    ::      Str_RHS                                                 !<
!   integer                                                       ::      i
!   do i = 1,size(Str)
!     call GetStrRHS( Str(i), Str_RHS(i) )
!   end do
! End Subroutine

!<==============================================================================================================
!> @brief       Removes characters at the LHS of a given character
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure removes all characters at the LHS of a given character. \n
!! The result is stored in a new output string called \c StrOut. \n
!! An intermediary local string \c StrLoc is used because the \c Split procedure has side effects on the its
!! first argument. \n
!! After the call to the \c Split procedure, the \c StrLoc and \c StrOut variables have all characters located
!! to the right and left side of the \c Separator character respectively.
!<==============================================================================================================
Function RemoveLeftChar( String, Separator, i_Debug ) result(StrOut)
  character(*)                                          ,intent(in)     ::      String                          !< Input character string
  character(*)                                          ,intent(in)     ::      Separator                          !< Separation character string
  logical       ,optional                               ,intent(in)     ::      i_Debug
  character(:)  ,allocatable                                            ::      StrOut                          ! Output character string
  character(:)  ,allocatable                                            ::      StrLoc                          ! Local copy of the input character string
  logical                                                               ::      i_Debug_Loc
  i_Debug_Loc = .false.; if ( present(i_Debug) ) i_Debug_Loc = i_Debug
  StrLoc        =       String
  StrOut        =       String
  if ( len_trim(String) == 0 ) return                                                                           ! If empty input string, then exiting the procedure
  call Split( String, StrOut, Separator, StrLoc, i_Debug=i_Debug )
End Function

!<==============================================================================================================
!> @brief       Parses a string into multiple string separated by a given character
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure parses the input string into sub-string based on the delimiters contained in a separation string. \n
!! Those sub-string are stored in the output vector string which is allocated inside the procedure. \n
!! The optional input variable \c EscRHS corresponds to a RHS escape character. \n
!<==============================================================================================================
Subroutine Parse( String_Input, Separator, String_Output, EscRHS, Ignore_Between )
  character(*)                                          ,intent(in)     ::      String_Input                    !< Input scalar character string to be proceeded
  character(*)                                          ,intent(in)     ::      Separator                       !< Separation character string
  character(:)  ,dimension(:)   ,allocatable            ,intent(out)    ::      String_Output                   !< Output vector character string containing sub-strings
  character(1)  ,optional       ,dimension(:)           ,intent(in)     ::      EscRHS                          !< RHS escape character
  character(*)  ,optional       ,dimension(:)           ,intent(in)     ::      Ignore_Between
  character(:)  ,allocatable                                            ::      String                    ! Local copy of the input string
  character(:)  ,allocatable                                            ::      String_LHS, String_RHS
  character(:)  ,dimension(:)   ,allocatable                            ::      String_Tmp                      ! Temporary string array
  integer ::  Length
  String        =       Compact(String_Input)                                                                   ! Storing local value of the input character string and compacting it
  if ( len_trim(String) == 0 ) return                                                                           ! If empty input string, then exiting the procedure
  do                                                                                                            ! Infinite loop on all characters of the input string
    if ( len_trim(String)==0 ) exit                                                                             ! If the current local character string is empty, then exiting the procedure
    call Split( String, String_LHS, Separator, String_RHS, EscRHS=EscRHS, Ignore_Between=Ignore_Between )       ! Spliting the string into a left and right string separated by the Separator character
    String      =       String_RHS                                                                              ! Setting the remaining string to be proceed to the RHS string
    if ( .not.allocated(String_Output) ) then                                                                   ! If un-allocated output String_Output variable, then the first sub-string to being proceeded
#ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
      allocate( character(len(String_LHS)) :: String_Output(1) )
      String_Output = [ String_LHS ]
#else
      allocate( String_Output, source = [ String_LHS ] )                                                         ! And so allocating the local list of string to unity and storing only the current value (since no previous values exists)
#endif
    else                                                                                                        ! If allocated output String_Output variable, then sub-strings have already been proceeded
#ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
      Length    =   max(len(String_Output),len(String_LHS))
      allocate( character(Length) :: String_Tmp(size(String_Output)+1) )
      String_Tmp = [ String_Output, String_LHS ]
      call move_alloc( String_Tmp, String_Output )                                                               ! Transfering allocation from temporary to final variable
#else
      allocate( String_Tmp, source = [ String_Output, String_LHS ] )                                            ! And so allocating the local list of string and storing both the previous and new values
      call move_alloc( String_Tmp, String_Output )                                                               ! Transfering allocation from temporary to final variable
#endif
    end if                                                                                                      ! End if cas on allocation status
  end do                                                                                                        ! End of loop
  String_Output(:)      =       VecTrim_1D(String_Output)
End Subroutine


!<==============================================================================================================
!> @brief       Converts multiple spaces and tabs to single spaces and removes left-hand spaces
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure converts multiple spaces and tabs to single spaces and removes left-hand spaces. \n
!<==============================================================================================================
Pure Function Compact( StrInp ) result(StrOut)
  character(*)                                          ,intent(in)     ::      StrInp                          !< Input character string
  character(len(StrInp))                                                ::      StrOut                          ! Output character string
  character(1)                                                          ::      char1                           ! String single character
  logical                                                               ::      i_space                         ! Spacing indicator
  integer                                                               ::      i, k                            ! Character index
  integer                                                               ::      i_ASCII                         ! Interger ASCII code for a given characters
  integer                                                               ::      LenStr                          ! Length of input string (without trailing blank characters)
  LenStr        =       len_trim(StrInp)                                                                        ! Getting the length of the input string (without trailing blank characters)
  StrOut        =       ''                                                                                      ! Initialization of the output string
  i_space       =       .false.                                                                                 ! Initialization of the space indicator
  k             =       0                                                                                       ! Initialization of the character index
  do i = 1,LenStr                                                                                               ! Loop on all string characters
    char1       =       StrInp(i:i)                                                                             ! Store current character
    i_ASCII     =       iachar(char1)                                                                           ! Getting the code for the ASCII character of the current character
    select case(i_ASCII)                                                                                        ! Select case according ASCII character code
    case (i_ASCII_VertTab,i_ASCII_Space)                                                                        ! If the character is a vertical tabulation or a space
      if (.not.i_space) then                                                                                    ! If no space have yet been treated
        k       =       k + 1                                                                                   ! Incrementing the output character string index
        StrOut(k:k)     =       ' '                                                                             ! Setting an space in the current character
      end if                                                                                                    ! End of if case on spacing indicator
      i_space   =       .true.                                                                                  ! Setting the spacing indicato
    case (33:)                                                                                                  ! If the character is not a control character
      k         =       k + 1                                                                                   ! Incrementing the output character string index
      StrOut(k:k)       =       char1                                                                           ! Setting the current charatcer to the input character
      i_space   =       .false.                                                                                 ! Setting the spacing indicato
    end select                                                                                                  ! End select on ASCII character code
  end do                                                                                                        ! End loop on string characters
  StrOut        =       adjustl(StrOut)                                                                         ! Removing left-hand spaces
End Function

!<==============================================================================================================
!> @brief       Converts tabs to spaces
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure converts tabs to spaces. \n
!<==============================================================================================================
Pure Function Tab2Space ( StrInp )      result(StrOut)
  character(*)                                          ,intent(in)     ::      StrInp                          !< Input character string
  character(len(StrInp))                                                ::      StrOut                          ! Output character string
  integer                                                               ::      i                               ! Character index
  StrOut        =       StrInp                                                                                  ! Initialization of the output string
  forall (i=1:len_trim(StrInp),iachar(StrOut(i:i))==i_ASCII_VertTab)    StrOut(i:i) = ' '                       ! If the current character is a vertical tabulation, then convert tabulation into space character
End Function

!<==============================================================================================================
!> @brief       Removes spaces and tabulation
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure removes spaces and tabulation. \n
!<==============================================================================================================
Pure Function RemoveSpace ( StrInp )    result(StrOut)
  character(*)                                          ,intent(in)     ::      StrInp                          !< Input character string
  character(len(StrInp))                                                ::      StrOut                          ! Output character string
  integer                                                               ::      i, k                            ! Character index
  integer                                                               ::      i_ASCII                         ! Interger ASCII code for a given characters
  StrOut        =       ''                                                                                      ! Initialization of the output string
  k             =       0                                                                                       ! Initialization of the character index
  do i = 1,len_trim(StrInp)                                                                                     ! Loop on all string characters
    i_ASCII     =       iachar(StrInp(i:i))                                                                     ! Getting the code for the ASCII character of the current character
    if  (i_ASCII == i_ASCII_Space)      cycle                                                                   ! If the character is a space, then going to the next character
    k           =       k + 1                                                                                   ! Incrementing the output character string index
    StrOut(k:k) =       StrInp(i:i)                                                                             ! Setting the current charatcer to the input character
  end do                                                                                                        ! End loop on string characters
End Function

!<==============================================================================================================
!> @brief       Splits a string into a RHS and A LHS part separated by a given string
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure splits a string into a LHS and A RHS part separated by a given string. \n
!! The first instance of a character from separation string \c Separator is located in the input string \c RHS. \n
!! The characters before the found delimiter are output in \c LHS, it is the Left-Hand Side string. \n
!! The characters after  the found delimiter are output in \c RHS, it is the Right-Hand-Side string. \n
!! The optional input strings contain the LHS and RHS escape characters which can be used for escaping
!! characters.
!! @TODO: It maigh be a good idea to invert the arguments !!!
!<==============================================================================================================
Subroutine Split( Input_String, LHS, Separator, RHS, EscLHS, EscRHS, i_Length_Eq_Inp, Ignore_Between, i_Debug )

  character(*)                                          ,intent(in)     ::      Input_String                          !< Input character string to be splitted into a LHS and a RHS parts
  character(:)  ,allocatable                            ,intent(out)    ::      LHS                             !< Left-hand-side part of the splitted string
  character(*)                                          ,intent(in)     ::      Separator                       !< Separation character string
  character(:)  ,allocatable                            ,intent(out)    ::      RHS                             !< Right-hand-side part of the splitted string
  character(1)  ,optional                               ,intent(in)     ::      EscLHS                          !< Left-hand-side escape character string
  character(1)  ,optional       ,dimension(:)           ,intent(in)     ::      EscRHS                          !< Right-hand-side escape character string
  logical       ,optional                               ,intent(in)     ::      i_Length_Eq_Inp                 ! Indicator to force the length of the output string to equal the one of the input string (required because otherwise the error:: fort: (4): Variable STRLHS has substring ending point 1 which is greater than the variable length of 0) Mayu be possible to
  character(*)  ,optional       ,dimension(:)           ,intent(in)     ::      Ignore_Between
  logical       ,optional                               ,intent(in)     ::      i_Debug

  logical                                                               ::      i_Debug_Loc
  logical                                                               ::      i_EscLHS                        ! Escape character indicator
  logical                                                               ::      i_EscRHS                        ! Escape character indicator
  logical                                                               ::      i_protect                       ! Character protection indicator
  logical                                                               ::      i_Length_Eq_Inp_Local
  logical                                                               ::      i_Ignore_Between
!   character(:)  ,allocatable    ,dimension(:)                           ::      Char_Ignore_Between
  integer                                                               ::      iStrInp                         ! Index of the input string
  integer                                                               ::      iRHS                            ! Index of the RHS string
  integer                                                               ::      iLHS                            ! Index of the LHS string
  integer                                                               ::      iSepIni                         ! Initial index of the separation character in the input string
  integer                                                               ::      iSepFin                         ! Final index of the separation character in the input string
  character(1)                                                          ::      Char1                           ! String single character
  character(1)                                                          ::      Char1_Next                      ! Next character
  character(len=len_trim(Input_String))                                 ::      StrInp                          ! Local copy of the input string
  integer                                                               ::      iEscRHS                         ! Index of element in the array of RHS escape character

  character(:)  ,allocatable                                            ::      String
  character(:)  ,allocatable                                            ::      Ignore_CharIni                  ! Initial character from which the separation character should be ignored
  character(:)  ,allocatable                                            ::      Ignore_CharFin                  ! final character   to   which the separation character should be ignored

  integer                                                               ::      Position_Ignore_CharIni
  integer                                                               ::      Position_Ignore_CharFin
  integer                                                               ::      i

  i_Debug_Loc = .false.; if ( present(i_Debug) ) i_Debug_Loc = i_Debug
!   i_Debug_Loc   =       .true.
!  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Entering')")
!  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: String to be splitted: ',a)") trim(Input_String)
!  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Separation character:  ',a)") trim(Separator)


  RHS   =       Input_String

! ==============================================================================================================
!       PROCESSING OPTIONAL INPUT ARGUMENTS
! ==============================================================================================================
  i_EscRHS              =       .false.                                                                         ! Initialization of the RHS escape character indicator
  i_EscLHS              =       .false.                                                                         ! Initialization of the LHS escape character indicator
  i_Length_Eq_Inp_Local =       .false.
  i_Ignore_Between      =       .false.
  if ( present(EscRHS) )                i_EscRHS = .true.                                                       ! If a RHS escape character is provided, the activating the RHS escape indicator
  if ( present(EscLHS) )                i_EscLHS = .true.                                                       ! If a LHS escape character is provided, the activating the LHS escape indicator
  if ( present(i_Length_Eq_Inp) )       i_Length_Eq_Inp_Local = i_Length_Eq_Inp
  if ( present(Ignore_Between) )        i_Ignore_Between = .true.
! !  if (i_Debug_Loc) then
!     write(LogUnit,"(12x,'[Split]: Dealing with optional input arguments')")
!     write(LogUnit,"(12x,'[Split]:  -> i_EscRHS              = ',g0)") i_EscRHS
!     write(LogUnit,"(12x,'[Split]:  -> i_EscLHS              = ',g0)") i_EscLHS
!     write(LogUnit,"(12x,'[Split]:  -> i_Length_Eq_Inp_Local = ',g0)") i_Length_Eq_Inp_Local
!     write(LogUnit,"(12x,'[Split]:  -> i_Ignore_Between      = ',g0)") i_Ignore_Between
!   end if
! ==============================================================================================================




  if ( i_EscRHS ) then
    RHS      =       RemoveSpace( adjustl(RHS) )                                                          ! Removing blanks and adjusting the to LHS
  else
!     RHS    =       Compact( adjustl(RHS) )                                                              ! Removing blanks and adjusting the to LHS
  end if



  allocate( character(len(RHS)) :: LHS )
  LHS(:)        =       ''                                                                                   ! Initialization of the LHS character string



  iSepIni   =       index( RHS, Separator )                                                                 ! Index of initial position of the separation character in the input character
  iSepFin   =       iSepIni + len_trim(Separator) - 1                                                      ! Index of final   position of the separation character in the input character

!  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: iSepIni = ',i0)") iSepIni
!  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: iSepFin = ',i0)") iSepFin

  if ( i_Ignore_Between ) then
  !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Ignoring the separation character in-between some given characters')")
  !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: size(Ignore_Between) = ',g0)") size(Ignore_Between)
  !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Ignore_Between) = ',*(a,3x))") Ignore_Between
    if ( size(Ignore_Between) == 1 ) then
      Ignore_CharIni    =       Ignore_Between(1)
      Ignore_CharFin    =       Ignore_Between(1)
    else if ( size(Ignore_Between) >= 2 ) then
      Ignore_CharIni    =       Ignore_Between(1)
      Ignore_CharFin    =       Ignore_Between(2)
    end if
  !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Ignore_CharIni = ',g0)") Ignore_CharIni
  !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Ignore_CharFin = ',g0)") Ignore_CharFin
!
!     Position_Ignore_CharIni     =       index( RHS, Ignore_CharIni )
!   !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Position of Ignore_CharIni in the string: Position_Ignore_CharIni = ',g0)") Position_Ignore_CharIni
!
!     Position_Ignore_CharFin     =       index( RHS, Ignore_CharFin )
!   !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Position of Ignore_CharFin in the string: Position_Ignore_CharFin = ',g0)") Position_Ignore_CharFin


!     This only works for ignore characters of 1 character
    Position_Ignore_CharIni     =       0
    Position_Ignore_CharFin     =       0
    do i = 1,len_trim(RHS)
      Char1     =       RHS(i:i)
      if ( (Position_Ignore_CharIni == 0) .and.  (Char1 == Ignore_CharIni) ) Position_Ignore_CharIni = i
      if ( (Position_Ignore_CharFin == 0) .and.  (Char1 == Ignore_CharFin) .and. (Position_Ignore_CharIni/=i) ) Position_Ignore_CharFin = i
    end do

  !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Position of Ignore_CharIni in the string: Position_Ignore_CharIni = ',g0)") Position_Ignore_CharIni
  !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Position of Ignore_CharFin in the string: Position_Ignore_CharFin = ',g0)") Position_Ignore_CharFin


    if ( ( Position_Ignore_CharIni < iSepIni ) .and. ( iSepIni < Position_Ignore_CharFin ) ) then
    !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: The separation character is inside the ignore set of characters => Finding the next splitting character')")

      String            =       RHS(Position_Ignore_CharFin+1:)
    !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: New string where to find the position of the separation character: String = ',g0)") String

      iSepIni   =       index( String, Separator )                                                                  ! Index of initial position of the separation character in the input character
      iSepFin   =       iSepIni + len_trim(Separator) - 1                                                      ! Index of final   position of the separation character in the input character

    !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Position in String: iSepIni = ',i0)") iSepIni
    !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Position in String: iSepFin = ',i0)") iSepFin

!       if iSepIni == 0 , then the separation character is not present in the string (outside the ignore character)
!       So the zero value is kept so that the procedure is exited below.
      if ( iSepIni /= 0 ) then
        iSepIni   =       iSepIni + (Position_Ignore_CharFin)
        iSepFin   =       iSepIni + len_trim(Separator) - 1                                                      ! Index of final   position of the separation character in the input character
      end if


    !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Original string where to find the position of the separation character: RHS = ',g0)") RHS
    !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Position in RHS: iSepIni = ',i0)") iSepIni
    !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Position in RHS: iSepFin = ',i0)") iSepFin

!     !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: String = ',g0)") RHS()
    end if



!     if ( Inf < Val )
!     if ( Val < Sup )
!     if ( Inf < Val ) .and. ( Val < Sup )        !     Inf < Val < Sup
!     Val = [Inf:Sup]
!     indicator   =       val .inside. [Inf,Sup]


!     if ( iSepIni )


  end if


! ==============================================================================================================
!       TREATING CASES FOR WHICH THE INPUT STRING DOES NOT NEED ANY SPLITTING
! ==============================================================================================================
! This section deals with the cases when no splitting is required. There are 3 situations for which the input
! string does not need any splitting. For all these cases the input string is stored in the output LHS string
! variable 'LHS' and the output RHS string variable 'RHS' is set to an empty string.
! This ensure the "Parse" calling procedure to stop its splitting iteration, the iteration being stopped when
! RHS string variable 'RHS' has a zero length. The 3 situations are:
!  - When the input character correspond to an empty string (obvious)
!  - When the separation character is absent from the input string
!  - When the separation character corresponds to the last character of the input string and that the
!    separation character is an RHS escape character (yes, it's a tricky one !)
! ==============================================================================================================
  if ( len_trim(RHS) == 0 ) then                                                                                ! If the input string is empty, then exiting the procedure without changing the input string
  !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Input string is empty => Exiting the procedure')")           ! Debugging
    return                                                                                                      ! Exiting the procedure
  end if                                                                                                        ! End if case on input string length
! --------------------------------------------------------------------------------------------------------------
  if ( iSepIni == 0 ) then                                                                                      ! If the separation character is absent from the input string, then
  !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: The separation character is absent from the input string => Exiting the procedure')") ! Debugging
    LHS      =       RHS                                                                                        ! Storing the input string in the LHS string
    RHS      =       ''                                                                                         ! Setting the RHS string to an empty string
    return                                                                                                      ! Exiting the procedure
  end if                                                                                                        ! End of if case on initial separation character index
! --------------------------------------------------------------------------------------------------------------
  if ( (i_EscRHS) .and. (iSepIni == len_trim(RHS)) )  then                                                      ! If a RHS escape character exists and if the last character of the input string corresponds to the separation character
  Char1 =       RHS(iSepIni:iSepIni)                                                                            ! Storing the last character
    do iEscRHS = 1,size(EscRHS)                                                                                 ! Loop on all RHS escape character
      if        (Char1 /= EscRHS(iEscRHS))      cycle                                                           ! Going to the next RHS escape character if the last character from the input string does not correspond to the current RHS escape character
        LHS  =       RHS                                                                                        ! Storing the input string in the LHS string
        RHS  =       ''                                                                                         ! Setting the RHS string to an empty string
        return                                                                                                  ! Exiting the procedure
    end do                                                                                                      ! End loop on RHS escape character
  end if                                                                                                        ! End of if case on initial separation character index
! ==============================================================================================================





  StrInp        =       RHS                                                                                     ! Storing the input string in a new variable
  RHS(:)        =       ''                                                                                      ! Initialization of RHS string
  iLHS          =       0                                                                                       ! Initialization of LHS string index
  iRHS          =       0                                                                                       ! Initialization of RHS string index
  iStrInp       =       0                                                                                       ! Initialization of input string index
  i_protect     =       .false.                                                                                 ! Initialization of the character protection indicator

!  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Loop on all characters inside the string: ',a)") StrInp
  do                                                                                                            ! Infinit loop on all characters of the input string

    iStrInp     =       iStrInp + 1                                                                             ! Incrementing the input string index
    if  (iStrInp > len_trim(StrInp))    exit                                                                    ! Exiting the loop if the the index is greater than the string length
    Char1       =       StrInp(iStrInp:iStrInp)                                                                 ! Getting the current character
  !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]:   iStrInp = ',i3,3x,'Char1 = ',a)") iStrInp, Char1           ! Debugging


!   Case of a LHS escape character
!   ------------------------------
    if (i_EscLHS) then                                                                                          ! If an escape character exists, then
    !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]:   Case of a LHS escape character')")                       ! Debugging
      if (i_protect) then                                                                                       ! If the current character is protected, then
        i_protect       =       .false.                                                                         ! Setting off the character protection mode
        iLHS            =       iLHS + 1                                                                        ! Index incrementation
        LHS(iLHS:iLHS)  =       Char1                                                                           ! Storing the current character in the LHS string
        cycle                                                                                                   ! Going to the next character
      end if                                                                                                    ! End of if case on character protection
      if ( Char1 == EscLHS ) then                                                                               ! If the current character corresponds to the escape character, then
        iLHS            =       iLHS + 1                                                                        ! Index incrementation
        LHS(iLHS:iLHS)  =       Char1                                                                           ! Storing the current character in the LHS string
        i_protect       =       .true.                                                                          ! Activating the character protection mode
        cycle                                                                                                   ! Going to the next character
      end if                                                                                                    ! End of if case on escape character
    end if                                                                                                      ! End of if case on existence of an escape character

!   Case of a RHS escape character
!   ------------------------------
    if (i_EscRHS) then                                                                                          ! If an RHS escape character exists, then
    !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]:   Case of a RHS escape character')")                       ! Debugging
      if ( iStrInp == iSepIni ) then                                                                            ! If the current input string index correspond to the initial separation index
        Char1_Next      =       StrInp(iStrInp+1:iStrInp+1)                                                     ! Storing the value of the next character
        do iEscRHS = 1,size(EscRHS)                                                                             ! Loop on all RHS escape characters
          if    (Char1_Next == EscRHS(iEscRHS)) then                                                            ! If the next character corresponds to the RHS escape character, then
             iSepIni        =       iSepIni + index( StrInp(iStrInp+1:), Separator )                            ! Modyfiying the index of initial position of the separation character in the input character: Moving to the next separation character index
             iSepFin        =       iSepIni + len_trim(Separator) - 1                                           ! Modyfiying the index of final   position of the separation character in the input character
          end if                                                                                                ! End of if case on escape character
        end do                                                                                                  ! End loop on RHS escape characters
      end if                                                                                                    ! End if case on input string index
    end if                                                                                                      ! End of if case on RHS escape indicator

    if ( iStrInp < iSepIni ) then                                                                               ! If the current input string index is lower than the initial separation index
    !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: iStrInp < iSepIni')")                                      ! Debugging
      iLHS              =       iLHS + 1                                                                        ! Index incrementation
      LHS(iLHS:iLHS)    =       Char1                                                                           ! Storing the current character in the LHS string
      cycle                                                                                                     ! Going to the next character
    end if                                                                                                      ! End if case on input string index

    if (iStrInp > iSepFin) then                                                                                 ! If the current input string index is greater than the final separation index
    !  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: iStrInp > iSepFin')")                                      ! Debugging
      iRHS              =       iRHS + 1                                                                        ! Index incrementation
      RHS(iRHS:iRHS)    =       Char1                                                                           ! Storing the current character in the LHS string
      cycle                                                                                                     ! Going to the next character
    end if                                                                                                      ! End if case on input string index

  end do                                                                                                        ! End loop on character of the input string
!  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: End loop')")

!  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: LHS = ',a)") LHS
!  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: RHS = ',a)") RHS

  RHS           =       adjustl(RHS)                                                                         ! Removing initial spaces
  RHS           =       trim(RHS)
  LHS           =       trim(LHS)

!  if (i_Debug_Loc) write(LogUnit,"(12x,'[Split]: Exiting')")

End Subroutine

!<==============================================================================================================
!> @brief       Converts a string into upper case
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure converts a string into upper case. \n
!! @COMPILER_BUG: The "elemental" attribute is not used for the "UpperCase" procedure since it prevents to use this procedure
!! with allocatable defered-length vector character string variables.
!! At least, with ifort version 16.0.0 an ICE is obtained.
!! This bug should be reported to the intel forum.
!! A simple workaround is to defined 2 procedures, one for scalar input variable and one for vectors, both being accessible through a generic binding name
!! Note that the ICE is only for ifort v16, everything goes smoothly with v15.
!<==============================================================================================================
! Pure Elemental Function UpperCase(StrInp) result(StrOut)
Pure Function UpperCase_0D(StrInp) result(StrOut)
  character(*)                  ,intent(in)     ::      StrInp                                                  !<
  character(len(StrInp))                        ::      StrOut
  integer                                       ::      i, ilen, ioffset, iquote, iav, iqc
  ilen          =       len_trim(StrInp)
  ioffset       =       iachar('A') - iachar('a')
  iquote        =       0
  StrOut        =       StrInp
  do i = 1,ilen
    iav=iachar(StrInp(i:i))
    if(iquote==0 .and. (iav==34 .or.iav==39)) then
      iquote    =       1
      iqc       =       iav
      cycle
    end if
    if(iquote==1 .and. iav==iqc) then
      iquote    =       0
    cycle
    end if
    if (iquote==1) cycle
    if(iav >= iachar('a') .and. iav <= iachar('z')) then
      StrOut(i:i)       =       achar(iav+ioffset)
    else
      StrOut(i:i)       =       StrInp(i:i)
    end if
  end do
End Function

Pure Function UpperCase_1D(StrInp) result(StrOut)
  character(*)                  ,dimension(:)                   ,intent(in)     ::      StrInp                                                  !<
  character(len(StrInp))        ,dimension(size(StrInp))                        ::      StrOut
  integer                                       ::      i
  do i = 1,size(StrInp)
    StrOut(i)   =       UpperCase(StrInp(i))
  end do
End Function

!<==============================================================================================================
!> @brief       Converts a string into lower case
!> @author      Bruno LOPEZ, blopez@ipfn.ist.utl.pt
!> @date        08/02/12 - Bruno LOPEZ - Initial creation of procedure
!<==============================================================================================================
!> @details
!! This procedure converts a string into lower case. \n
!<==============================================================================================================
Pure Function LowerCase (StrInp)        result(StrOut)
  character(*)                  ,intent(in)     ::      StrInp                                                  !<
  character(len(StrInp))                        ::      StrOut
  integer                                       ::      i, ilen, ioffset, iquote, iav, iqc
  ilen          =       len_trim(StrInp)
  ioffset       =       iachar('A')-iachar('a')
  iquote        =       0
  StrOut        =       StrInp
  do i = 1,ilen
    iav =       iachar(StrInp(i:i))
    if(iquote==0 .and. (iav==34 .or.iav==39)) then
      iquote    =       1
      iqc       =       iav
      cycle
    end if
    if(iquote==1 .and. iav==iqc) then
      iquote    =       0
      cycle
    end if
    if (iquote==1) cycle
    if(iav >= iachar('A') .and. iav <= iachar('Z')) then
      StrOut(i:i)       =       achar(iav-ioffset)
    else
      StrOut(i:i)       =       StrInp(i:i)
    end if
  end do
End Function
!
! Subroutine Convert_From_Character( InpVar, OutVar )
!   character(*)                                          ,intent(in)     ::      InpVar                          !< Input variable to be converted
!   class(*)                                              ,intent(out)    ::      OutVar                          !< Output variable resulting from the conversion
!   select type (OutVar)
!     type is (logical);    Value       =       Convert_To_Logical(InpVar)
!     type is (integer);    Value       =       Convert_To_Integer(InpVar)
!     type is (real(8));    Value       =       Convert_To_Real(InpVar)
!   end select
! End Subroutine


Pure Function String_To_String_0D( Var, Len, Pos, Fmt ) result(Str)
  character(*)                                          ,intent(in)     ::      Var                             !< String to be converted into a string (with possibly a different length
  integer       ,optional                               ,intent(in)     ::      Len                             !< Length of the output string
  character(1)  ,optional                               ,intent(in)     ::      Pos                             !< Position of the string: L:Left, C:Center, R:Right
  character(*)  ,optional                               ,intent(in)     ::      Fmt                             !< Format used to describe the number
  character(:)  ,allocatable                                            ::      Str                             !< Character corresponding to the input number
  character(1000)                                                       ::      Str_Loc                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Fmt_Loc                         ! Local value of the number format
  Fmt_Loc       =       "g0"                                                                                    ! Setting the default format
  if ( present(Fmt) ) Fmt_Loc = Fmt                                                                             ! Setting the input value for the format if present
  Fmt_Loc       =       "(" // Fmt_Loc // ")"                                                                   ! Adding the parentesis to the format
  write(Str_Loc,Fmt_Loc) Var                                                                                    ! Interger-to-string convertion
  Str   =       trim(Str_Loc)                                                                          ! Storing the string in the output variable with the right number of characters
!   if ( .Not. present(Fmt) ) Str   =       trim(adjustl(Str_Loc))                                                                          ! Storing the string in the output variable with the right number of characters
  if ( present(Len) ) then
    Str = Set_Length( Str, Len )
    if ( present(Pos) ) then
      if ( Pos == 'L') Str(:) = adjustl(Str)
!       if ( Pos == 'C')
      if ( Pos == 'R') Str(:) = adjustr(Str)
    end if
  end if
End Function

Pure Function Logical_To_String_0D( Var ) result(Str)
  logical                                               ,intent(in)     ::      Var                             !< Logical to be converted into a string
  character(:)  ,allocatable                                            ::      Str                             !< Character corresponding to the input number
  if ( Var ) then
    Str =       "T"
  else
    Str =       "F"
  end if
End Function


! REMARK:
! This procedure converts a integer number into a character string.
! Using character length allocation, the output character has a length wich exactly fits the number of digits
! of the input integer.
! An optional input argument <LeadingZeros> can be supplied to the procedure in order to force the number
! of leading zeros in the output character string.
! For example, using the Convert_To_String interface, one obtains the following output strings for the
! above calling sequences:
!       call Convert_To_String( 50, LeadingZeros=4 )    =>      "0050"
!       call Convert_To_String( 50 )                    =>      "50"
! @TODO: Check validity of input format
Pure Recursive Function Integer_To_String_0D( Var, LeadingZeros, Len, Pos, Fmt ) result(Str)
  integer                                               ,intent(in)     ::      Var                             !< Number to be converted into a string
  integer       ,optional                               ,intent(in)     ::      LeadingZeros                    !< Number of leading zeros
  integer       ,optional                               ,intent(in)     ::      Len                             !< Length of the output string
  character(1)  ,optional                               ,intent(in)     ::      Pos                             !< Position of the string: L:Left, C:Center, R:Right
  character(*)  ,optional                               ,intent(in)     ::      Fmt                             !< Format used to describe the number
  character(:)  ,allocatable                                            ::      Str                             !< Character corresponding to the input number
  character(1000)                                                       ::      Str_Loc                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Fmt_Loc                         ! Local value of the number format
  character(:)  ,allocatable                                            ::      NZeros
  NZeros        =       ""                                                                                      ! Initializing the number of leading zeros to an empty string => No leading zeros by default
  if ( present(LeadingZeros) ) NZeros = "." // Convert_To_String(LeadingZeros)                                  ! Setting the input value of leading zeros if present
  Fmt_Loc       =       "i0" // NZeros                                                                          ! Setting the default format
  if ( present(Fmt) ) Fmt_Loc = Fmt                                                                             ! Setting the input value for the format if present
  Fmt_Loc       =       "(" // Fmt_Loc // ")"                                                                   ! Adding the parentesis to the format
  write(Str_Loc,Fmt_Loc) Var                                                                                    ! Interger-to-string convertion
  Str   =       trim(Str_Loc)                                                                          ! Storing the string in the output variable with the right number of characters
  if ( .Not. present(Fmt) ) Str = trim(adjustl(Str_Loc))                                                                          ! Storing the string in the output variable with the right number of characters
  if ( present(Len) ) Str = Set_Length( Str, Len )
  if ( present(Pos) ) then
    if ( Pos == 'L') Str(:) = adjustl(Str)
    if ( Pos == 'R') Str(:) = adjustr(Str)
  end if
End Function

    Pure Function Real4_To_String_0D( Var, Len, Pos, Fmt ) result(String)
      real(4)                                               ,intent(in)     ::      Var                             !< Real number to be converted into a string
      integer       ,optional                               ,intent(in)     ::      Len                             !< Length of the output string
      character(1)  ,optional                               ,intent(in)     ::      Pos                             !< Position of the string: L:Left, C:Center, R:Right
      character(*)  ,optional                               ,intent(in)     ::      Fmt
      character(:)  ,allocatable                                            ::      String                          !< Character corresponding to the input number
  character(10000)                                                      ::      LongString                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      LocalFormat
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefaultFormat = "(*(g0,1x))"
  LocalFormat   =   "g0"                                                                                    ! Setting the default format
  if ( present(Fmt) ) LocalFormat = Fmt                                                                             ! Setting the input value for the format if present
  LocalFormat   =   "(" // LocalFormat // ")"                                                                   ! Adding the parentesis to the format
  write( LongString, LocalFormat, iostat=ios ) Var                                                                                        ! Interger-to-string convertion
  if ( ios /= 0 ) then
    LocalFormat =   DefaultFormat
    write( LongString, LocalFormat, iostat=ios ) Var
  end if
  String    =       trim(LongString)
  if ( .not. present(Fmt) ) String   =       Remove_Trailing_Zeros(trim(adjustl(LongString)))                                                   ! Storing the string in the output variable with the right number of characters
  if ( present(Len) )       String   =       Set_Length( String, Len )
  if ( present(Pos) ) then
    if ( Pos == 'L') String(:) = adjustl(String)
    if ( Pos == 'R') String(:) = adjustr(String)
  end if
    End Function

    Pure Function Real8_To_String_0D( Var, Len, Pos, Fmt ) result(String)
      real(8)                                               ,intent(in)     ::      Var                             !< Real number to be converted into a string
      integer       ,optional                               ,intent(in)     ::      Len                             !< Length of the output string
      character(1)  ,optional                               ,intent(in)     ::      Pos                             !< Position of the string: L:Left, C:Center, R:Right
      character(*)  ,optional                               ,intent(in)     ::      Fmt
      character(:)  ,allocatable                                            ::      String                          !< Character corresponding to the input number
  character(10000)                                                      ::      LongString                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      LocalFormat
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefaultFormat = "(*(g0,1x))"
  LocalFormat   =   "g0"                                                                                    ! Setting the default format
  if ( present(Fmt) ) LocalFormat = Fmt                                                                             ! Setting the input value for the format if present
  LocalFormat   =   "(" // LocalFormat // ")"                                                                   ! Adding the parentesis to the format
  write( LongString, LocalFormat, iostat=ios ) Var                                                                                        ! Interger-to-string convertion
  if ( ios /= 0 ) then
    LocalFormat =   DefaultFormat
    write( LongString, LocalFormat, iostat=ios ) Var
  end if
  String    =       trim(LongString)
  if ( .not. present(Fmt) ) String   =       Remove_Trailing_Zeros(trim(adjustl(LongString)))                                                   ! Storing the string in the output variable with the right number of characters
  if ( present(Len) )       String   =       Set_Length( String, Len )
  if ( present(Pos) ) then
    if ( Pos == 'L') String(:) = adjustl(String)
    if ( Pos == 'R') String(:) = adjustr(String)
  end if
    End Function

    Pure Function Real16_To_String_0D( Var, Len, Pos, Fmt ) result(String)
      real(16)                                              ,intent(in)     ::      Var                             !< Real number to be converted into a string
      integer       ,optional                               ,intent(in)     ::      Len                             !< Length of the output string
      character(1)  ,optional                               ,intent(in)     ::      Pos                             !< Position of the string: L:Left, C:Center, R:Right
      character(*)  ,optional                               ,intent(in)     ::      Fmt
      character(:)  ,allocatable                                            ::      String                             !< Character corresponding to the input number
  character(10000)                                                      ::      LongString                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      LocalFormat
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefaultFormat = "(*(g0,1x))"
  LocalFormat   =   "g0"                                                                                    ! Setting the default format
  if ( present(Fmt) ) LocalFormat = Fmt                                                                             ! Setting the input value for the format if present
  LocalFormat   =   "(" // LocalFormat // ")"                                                                   ! Adding the parentesis to the format
  write( LongString, LocalFormat, iostat=ios ) Var                                                                                        ! Interger-to-string convertion
  if ( ios /= 0 ) then
    LocalFormat =   DefaultFormat
    write( LongString, LocalFormat, iostat=ios ) Var
  end if
  String    =       trim(LongString)
  if ( .not. present(Fmt) ) String   =       Remove_Trailing_Zeros(trim(adjustl(LongString)))                                                   ! Storing the string in the output variable with the right number of characters
  if ( present(Len) )       String   =       Set_Length( String, Len )
  if ( present(Pos) ) then
    if ( Pos == 'L') String(:) = adjustl(String)
    if ( Pos == 'R') String(:) = adjustr(String)
  end if
    End Function

! Pure Function Real8_To_String_0D( Var, Len, Pos, Fmt ) result(Str)
!   real(8)                                               ,intent(in)     ::      Var                             !< Real number to be converted into a string
!   integer       ,optional                               ,intent(in)     ::      Len                             !< Length of the output string
!   character(1)  ,optional                               ,intent(in)     ::      Pos                             !< Position of the string: L:Left, C:Center, R:Right
!   character(*)  ,optional                               ,intent(in)     ::      Fmt
!   character(:)  ,allocatable                                            ::      Str                             !< Character corresponding to the input number
!   character(100)                                                        ::      Str_Loc                         ! Local character string required to store the number of a very large string before allocation of the output variable
!   character(:)  ,allocatable                                            ::      Fmt_Loc
!   Fmt_Loc       =       "g0"                                                                                    ! Setting the default format
!   if ( present(Fmt) ) Fmt_Loc = Fmt                                                                             ! Setting the input value for the format if present
!   Fmt_Loc       =       "(" // Fmt_Loc // ")"                                                                   ! Adding the parentesis to the format
!   write(Str_Loc,Fmt_Loc) Var                                                                                        ! Interger-to-string convertion
!   Str           =       trim(Str_Loc)
!   if ( .not. present(Fmt) ) Str   =       Remove_Trailing_Zeros(trim(adjustl(Str_Loc)))                                                   ! Storing the string in the output variable with the right number of characters
!   if ( present(Len) )       Str   =       Set_Length( Str, Len )
!   if ( present(Pos) ) then
!     if ( Pos == 'L') Str(:) = adjustl(Str)
!     if ( Pos == 'R') Str(:) = adjustr(Str)
!   end if
! End Function
!
! Pure Function Real_To_String_1D( Var, Fmt ) result(Str)
!   real(8)       ,dimension(:)                           ,intent(in)     ::      Var                             !< Real number to be converted into a string
!   character(*)  ,optional                               ,intent(in)     ::      Fmt
!   character(:)  ,allocatable    ,dimension(size(Var))                   ::      Str                             !< Character corresponding to the input number
!
!   do i = 1,size(Var)
!
! End Function

Pure Recursive Function Integer1D_To_String0D( Var, LeadingZeros, Len, Pos, Fmt ) result(Str)
  integer       ,dimension(:)                           ,intent(in)     ::      Var                             !< Number to be converted into a string
  integer       ,optional                               ,intent(in)     ::      LeadingZeros                    !< Number of leading zeros
  integer       ,optional                               ,intent(in)     ::      Len                             !< Length of the output string
  character(1)  ,optional                               ,intent(in)     ::      Pos                             !< Position of the string: L:Left, C:Center, R:Right
  character(*)  ,optional                               ,intent(in)     ::      Fmt                             !< Format used to describe the number
  character(:)  ,allocatable                                            ::      Str                             !< Character corresponding to the input number
  character(1000)                                                       ::      Str_Loc                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Fmt_Loc                         ! Local value of the number format
  character(:)  ,allocatable                                            ::      NZeros
  NZeros        =       ""                                                                                      ! Initializing the number of leading zeros to an empty string => No leading zeros by default
  if ( present(LeadingZeros) ) NZeros = "." // Convert_To_String(LeadingZeros)                                  ! Setting the input value of leading zeros if present
  Fmt_Loc       =       "*(i0" // NZeros // ",1x)"                                                              ! Setting the default format
  if ( present(Fmt) ) Fmt_Loc = Fmt                                                                             ! Setting the input value for the format if present
  Fmt_Loc       =       "(" // Fmt_Loc // ")"                                                                   ! Adding the parentesis to the format
  write(Str_Loc,Fmt_Loc) Var                                                                                    ! Interger-to-string convertion
  Str   =       trim(Str_Loc)                                                                          ! Storing the string in the output variable with the right number of characters
  if ( .Not. present(Fmt) ) Str = trim(adjustl(Str_Loc))                                                                          ! Storing the string in the output variable with the right number of characters
  if ( present(Len) ) Str = Set_Length( Str, Len )
  if ( present(Pos) ) then
    if ( Pos == 'L') Str(:) = adjustl(Str)
    if ( Pos == 'R') Str(:) = adjustr(Str)
  end if
End Function

Pure Function Set_Length( StrInp, Length ) result(StrOut)
  character(*)                                          ,intent(in)     ::      StrInp
  integer                                               ,intent(in)     ::      Length
  character(Length)                                                     ::      StrTmp
  character(:)  ,allocatable                                            ::      StrOut
!   allocate( character(Length) :: StrOut )
!   StrOut(:)     =       StrInp
  StrTmp        =       StrInp
  allocate( StrOut , source = StrTmp)
End Function

! @TODO: Implement the "Length_From_Format" procedure.
! This procedure returns a character string length from a given fortran format.
! For example,
!       Fmt                     Length
!       "(es15.8)"              15
!       "es15.8"                15
!       "(3(es15.8))"           45 ( = 3 * 15 )
Pure Function Length_From_Format( Fmt ) result(Length)
  character(*)                                          ,intent(in)     ::      Fmt
  integer                                                               ::      Length
  character(1)  ::      cdum
  Length        =       0
  cdum  =       Fmt
End Function


! This procedure returns a real format from the total length.
! Format = es<X>.<X-8>E3
! @TODO: Check that L-8 does not become negative
! Examples:
!   Length              =>              Format
!   15                                  es15.8E3
Pure Function Real_Format_From_Length( Length ) result( Fmt )
  integer                                               ,intent(in)     ::      Length
  character(:)  ,allocatable                                            ::      Fmt
  character(:)  ,allocatable                                            ::      L, Lm8
  L     =       Convert_To_String(Length)
  Lm8   =       Convert_To_String(Length - 8)
  Fmt   =       "es" // L // "." // Lm8 // "E3"
End Function



Pure Subroutine Add_Line_To_String( String, Line, At_End, At_Start, At_Position )

  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      String
  character(*)                                          ,intent(in)     ::      Line
  logical       ,optional                               ,intent(in)     ::      At_End  ! Default behavior
  logical       ,optional                               ,intent(in)     ::      At_Start
  integer       ,optional                               ,intent(in)     ::      At_Position

  integer                                                               ::      Length
  character(:)  ,dimension(:)   ,allocatable                            ::      String_tmp

  if ( .not. allocated(String) ) allocate( character(0) :: String(0) )
  Length        =       max( len(String), len(Line) )
  allocate( character(Length) :: String_tmp(size(String)+1) )

  if ( present(At_Start) ) then
  if ( At_Start ) then
    String_tmp(1)               =       Line
    String_tmp(2:size(String)+1)=       String
    call move_alloc( String_tmp, String )
    return
  end if
  end if

  String_tmp(1:size(String))    =       String
  String_tmp(size(String)+1)    =       Line
  call move_alloc( String_tmp, String )
End Subroutine

!
! Pure Subroutine Add_Line_To_String( String, Line, At_End, At_Start, At_Position )
!   character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      String
!   character(*)                                          ,intent(in)     ::      Line
!   logical       ,optional                               ,intent(in)     ::      At_End  ! Default behavior
!   logical       ,optional                               ,intent(in)     ::      At_Start
!   integer       ,optional                               ,intent(in)     ::      At_Position
!   integer                                                               ::      Length
!   character(:)  ,dimension(:)   ,allocatable                            ::      String_tmp
! #ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
!   integer                                                               ::      i
! #endif
!   if ( .not. allocated(String) ) allocate( character(0) :: String(0) )
!   Length        =       max( len(String), len(Line) )
!   allocate( character(Length) :: String_tmp(size(String)+1) )
!   if ( present(At_Start) ) then
!   if ( At_Start ) then
!     String_tmp(1)               =       Line
! #ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
!     do i = 2,size(String)+1
!       String_tmp(i)               =       String(i-1)
!     end do
! #else
!     String_tmp(2:size(String)+1)=       String
! #endif
!     call move_alloc( String_tmp, String )
!     return
!   end if
!   end if
! #ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
!   do i = 1,size(String)
!     String_tmp(i)               =       String(i)
!   end do
! #else
!   String_tmp(1:size(String))    =       String
! #endif
!   String_tmp(size(String)+1)    =       Line
!   call move_alloc( String_tmp, String )
! End Subroutine



Pure Subroutine Add_Lines_To_String( String, Lines, At_End, At_Start, At_Position )
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      String
  character(*)  ,dimension(:)                           ,intent(in)     ::      Lines
  logical       ,optional                               ,intent(in)     ::      At_End  ! Default behavior
  logical       ,optional                               ,intent(in)     ::      At_Start
  integer       ,optional                               ,intent(in)     ::      At_Position
  integer                                                               ::      Length
  character(:)  ,dimension(:)   ,allocatable                            ::      String_tmp
  if ( .not. allocated(String) ) allocate( character(0) :: String(0) )
  Length      =       max( len(String), len(Lines) )
  allocate( character(Length) :: String_tmp(size(String)+size(Lines)) )


  if ( present(At_Start) ) then
  if ( At_Start ) then
    String_tmp(1:size(Lines))                   =       Lines
    String_tmp(size(Lines)+1:size(String_tmp))  =       String
    call move_alloc( String_tmp, String )
    return
  end if
  end if

  String_tmp(1:size(String))                  =       String
  String_tmp(size(String)+1:size(String_tmp)) =       Lines
  call move_alloc( String_tmp, String )
End Subroutine




Subroutine Add_Element_If_Absent_Integer( Element, List )
  integer                                               ,intent(in)     ::      Element
  integer       ,dimension(:)   ,allocatable            ,intent(inout)  ::      List
  integer       ,dimension(:)   ,allocatable                            ::      List_tmp
!   integer       ,allocatable                            ::      List_tmp(:)
  if ( .not.allocated(List) ) then
    allocate( List(1) )
    List(1)       =       Element
  else
    if ( .not.Presence(Element,List) ) then
#ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
      allocate( List_tmp(size(List)+1) )
      List_tmp = [List,Element]
      call move_alloc( List_tmp, List )                                                                     ! Transfering the values from the temporary variable to the output variable
#else
      allocate( List_tmp, source=[List,Element] )
      call move_alloc( List_tmp, List )                                                                     ! Transfering the values from the temporary variable to the output variable
#endif
    end if
  end if
End Subroutine

Subroutine Add_Element_If_Absent_C0( Array, Element )
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      Array
  character(*)                                          ,intent(in)     ::      Element
  character(:)  ,dimension(:)   ,allocatable                            ::      List_Elements
  if ( .Not. allocated(Array) ) allocate( character(0) :: Array(0) )                                            ! Allocating the array if not allocated
  if ( .Not. Presence( UpperCase(Element), UpperCase(Array) ) ) then                                            ! If the element is not present in the array, then adding it
#ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
    allocate( character(max(len(Array),len(Element))) :: List_Elements(size(Array)+1) )
    List_Elements = [Array,Element]
    call move_alloc( List_Elements, Array )                                                                     ! Transfering the values from the temporary variable to the output variable
#else
    allocate( List_Elements, source = [Array,Element] )                                                         ! Source allocation of a temporary array containing the both the previous and the new elements
    call move_alloc( List_Elements, Array )                                                                     ! Transfering the values from the temporary variable to the output variable
#endif
  end if                                                                                                        ! End if case on presence of the element in the array
End Subroutine

Subroutine Add_Element_If_Absent_C1( Array, Elements )
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      Array
  character(*)  ,dimension(:)                           ,intent(in)     ::      Elements
  integer                                                               ::      i
  if ( .Not. allocated(Array) ) allocate( character(0) :: Array(0) )                                                            ! Allocating the array if not allocated
  do i = 1,size(Elements)
    call Add_Element_If_Absent( Array, Elements(i))
  end do                                                                       ! End if case on presence of the element in the array
End Subroutine




Pure Function Get_Number_Of_Digits( iNum ) result(NDig)
  integer                                                               ::      NDig                            !< Number of digits of the input integer number
  integer       ,intent(in)                                             ::      iNum                            !< Interger number whose number of digits is to be determined
  NDig  =       floor( log10( real( abs(iNum) ) ) ) + 1
End Function




Subroutine Remove_Duplicate( Array )
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      Array
  integer                                                               ::      i, j
  character(:)  ,dimension(:)   ,allocatable                            ::      List_Elements
  integer                                                               ::      Length, k
  i                     =       0
  do
    i                   =       i + 1
    if ( i >= size(Array) ) exit
    j                   =       i
    do
      j                 =       j + 1
      if ( j > size(Array) ) exit
      if ( Array(i) /= Array(j) ) cycle
!       write(*,"('[Remove_Duplicate]:  Remove element Array(j) = ',g0)")  Array(j)
      if ( size(Array) == 2 ) then
#ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
        allocate( character(len(Array)) :: List_Elements(1) )
        List_Elements(1) = Array(1)
        call move_alloc( List_Elements, Array )
#else
        allocate( List_Elements, source = [Array(1)] )
        call move_alloc( List_Elements, Array )
#endif
      else
#ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
        Length = 0
        do k = 1,size(Array)
          if ( k == j ) cycle
          Length = max( Length , len(Array) )
        end do
        allocate( character(Length) :: List_Elements(size(Array)-1) )
        List_Elements = [Array(1:j-1),Array(j+1:)]
        call move_alloc( List_Elements, Array )
#else
        allocate( List_Elements, source = [Array(1:j-1),Array(j+1:)] )
        call move_alloc( List_Elements, Array )
#endif
      end if
      j = j - 1
    end do
  end do


!   i             =       1                                                                                       ! Initializing the index of the lower elements
!   do                                                                                                            ! Looping over lower elements
!     j           =       i                                                                                       ! Initializing the index of the upper elements to the lower one
!     do                                                                                                          ! Looping over upper elements
!       j         =       j + 1                                                                                   ! Incrementing the index of the upper elements
!       if ( j > size(Array) ) exit                                                                               ! Exiting the loop on upper elements if all elements have been processed
!       if ( Array(i) /= Array(j) ) cycle                                                                         ! If the 2 elements are different, then skipping to the next upper elements
!       call Remove_Element( Array, j )                                                                           ! Removing the upper elements
!       j    =       j - 1                                                                                        ! Decrementing the index of the upper elements
!     end do                                                                                                      ! End loop on upper elements
!     i          =       i + 1                                                                                    ! Incrementing the index of the lower elements
!     if ( i >= size(Array) ) exit                                                                                ! Exiting the loop on lower elements if all elements have been processed
!   end do                                                                                                        ! End loop on lower elements


End Subroutine


! REMARK:
! The input argument corresponds to a character string corresponding to a real number.
! This procedure will remove the trailling zeros from the fractional part of the real number.
! If no fractional part is found, then the input character string is unchanged
!
! 0. 0.13400666E-11

! 1340066600000000
! E-11

Pure Function Remove_Trailing_Zeros( StrInp ) result(StrOut)
  character(*)                                          ,intent(in)     ::      StrInp
  character(:)  ,allocatable                                            ::      StrOut
  character(*)  ,parameter                                              ::      Separator_Int_Fra="."
  character(*)  ,parameter                                              ::      Separator_Fra_Exp="E"
  integer                                                               ::      iSpe
  character(:)  ,allocatable                                            ::      Integer_Part
  character(:)  ,allocatable                                            ::      Fractional_Part
  character(:)  ,allocatable                                            ::      Exponent_Part
  integer                                                               ::      i
  StrOut        =       UpperCase( StrInp )                                                                     ! The conversion to upper case is required to be sure that the exponent character is "E" and not "e"
  iSpe          =       index( StrOut, Separator_Int_Fra )                                                      ! Getting the index of the separation character between the integer and fractional parts (zero if any)
  if ( iSpe /= 0 ) then                                                                                         ! If a separation character is found
    Integer_Part        =       trim( StrOut(1:iSpe-1) )                                                        ! Getting the integer part of the real number, ie. all characters before the '.' character
    Fractional_Part     =       trim( StrOut(iSpe+1:)  )                                                        ! Getting the fractional part of the real number, ie. all characters after the '.' character
    Exponent_Part       =       ""                                                                              ! Initializing the exponent part
    iSpe          =       index( Fractional_Part, Separator_Fra_Exp )                                           ! Getting the index of the separation character between the fractional and exponent part (zero if any)
    if ( iSpe /= 0 ) then                                                                                       ! If a exponent part existe, then further splitting
      Exponent_Part     =       trim( Fractional_Part(iSpe+1:)  )                                               ! Getting the exponent part of the real number, ie. all characters after the 'E' character
      Fractional_Part   =       trim( Fractional_Part(1:iSpe-1)  )                                              ! Getting the fractional part of the real number, ie. all characters before the 'E' character
    end if                                                                                                      ! End if case on separation character presence
    if ( len_trim(Integer_Part) == 0 ) Integer_Part = "0"                                                       ! If no integer part (the input string has the format ".XXXX"), then explicitly setting the integer part is zero
    i   =       len(Fractional_Part) + 1                                                                        ! Initializing the character index to the length+1 of the fractional part string (beacause of backward processing of the string)
    do                                                                                                          ! Loop on all characters of the fractional part string
      i = i - 1                                                                                                 ! Incrementing the character index
      if ( i == 0 )                     exit                                                                    ! Exiting the loop if the 1st character the entire string has been processed, (that is, if i=0)
      if ( Fractional_Part(i:i)/="0" )  exit                                                                    ! Exiting the loop if current character does not correspond to a zero ("0")
      Fractional_Part(i:i)      =       " "                                                                     ! Replacing the zero character by a blank character
    end do                                                                                                      ! End loop on characters of the fractional part string
    Fractional_Part     =       trim( Fractional_Part )                                                         ! Removing trailling blanks
    if ( len_trim(Fractional_Part) == 0 ) Fractional_Part = "0"                                                 ! If the no fractional part (the input string has the format "X.000...000"), then setting it to zero so that it become "X.0"
    StrOut      =       Integer_Part                                                                            ! Reconstructing the character string from the integer ...
    if ( len_trim(Fractional_Part) /= 0 ) StrOut = StrOut // Separator_Int_Fra // Fractional_Part               ! ... the fractional part if required
    if ( len_trim(Exponent_Part) /= 0 )   StrOut = StrOut // Separator_Fra_Exp // Exponent_Part                 ! ... and the exponent part if required
  end if                                                                                                        ! End if case on separation character presence
End Function


Pure Function Max_Len_Trim_0D( Strings ) result( Length )
  character(*)                                  ,intent(in)             ::      Strings                         !< Character string
  integer                                                               ::      Length                          !< Maximum length without trailling blanks
  Length        =       len_trim(Strings)                                                                       ! Setting the maximum length
End Function

Pure Function Max_Len_Trim_1D( Strings ) result( Length )
  character(*)  ,dimension(:)                   ,intent(in)             ::      Strings                         !< Array of character string
  integer                                                               ::      Length                          !< Maximum length without trailling blanks along all elements of the input string array
  integer                                                               ::      i                               ! Index of string' elements
  Length        =       0                                                                                       ! Initialization of maximum length of string
  do i = 1,size(Strings,1)                                                                                      ! Loop on all elements
    Length      =       max( Length, len_trim(Strings(i)) )                                                     ! Setting the maximum length
  end do                                                                                                        ! End loop on all elements
End Function

Pure Function Max_Len_Trim_2D( Strings ) result( Length )
  character(*)  ,dimension( :, : )              ,intent(in)             ::      Strings                         !< Array of character string
  integer                                                               ::      Length                          !< Maximum length without trailling blanks along all elements of the input string array
  integer                                                               ::      i, j                            ! Index of string' elements
  Length        =       0                                                                                       ! Initialization of maximum length of string
  do j = 1,size(Strings,2)                                                                                      ! Loop on all elements
  do i = 1,size(Strings,1)                                                                                      ! Loop on all elements
    Length      =       max( Length, len_trim(Strings(i,j)) )                                                   ! Setting the maximum length
  end do                                                                                                        ! End loop on all elements
  end do                                                                                                        ! End loop on all elements
End Function
!
!  Pure Function VecTrim_0D( Inp ) result(Out)
!   character(*)                                          ,intent(in)     ::      Inp                             !< Input character string to be trimed
!   character(:)  ,allocatable                                            ::      Out                             !< Output character string
!   Out =       trim(Inp)                                                                                         ! Setting the maximum length
! End Function
!
!  Pure Function VecTrim_1D( Inp ) result(Out)
!   character(*)  ,dimension(:)                           ,intent(in)     ::      Inp                             !< Input character string to be trimed
!   character(Max_Len_Trim(Inp)) ,dimension(size(Inp))                    ::      Out                             !< Output character string of length Max_Len_Trim(Inp) and size size(Inp)
!   integer                                                               ::      i
!   do i = 1,size(Inp)
!     Out(i)    =       trim(Inp(i))
!   end do
! End Function
!
!  Pure Function VecTrim_2D( Inp ) result(Out)
!   character(*)  ,dimension( :, : )                      ,intent(in)     ::      Inp                             !< Input character string to be trimed
!   character(Max_Len_Trim(Inp)) ,dimension( size(Inp,1), size(Inp,2) )   ::      Out                             !< Output character string of length Max_Len_Trim(Inp)
!   integer                                                               ::      i, j
!   do j = 1,size(Inp,2)
!   do i = 1,size(Inp,1)
!     Out(i,j)    =       trim(Inp(i,j))
!   end do
!   end do
! End Function

Function VecTrim_0D( Input_String ) result(Output_String)
  character(:)                          ,allocatable    ,intent(in)     ::      Input_String
  character(:)                          ,allocatable                    ::      Output_String
  Output_String =       trim( Input_String )
End Function

Function VecTrim_1D( Input_String ) result(Output_String)
  character(:)  ,dimension(:)           ,allocatable    ,intent(in)     ::      Input_String
  character(:)  ,dimension(:)           ,allocatable                    ::      Output_String
  integer                                                               ::      i
  integer                                                               ::      Length
  integer                                                               ::      Size1
  Length        =       Max_Len_Trim(Input_String)
  Size1         =       size(Input_String,1)
  allocate( character(Length) :: Output_String(Size1) )
  do i = 1,Size1
    Output_String(i)    =       trim( Input_String(i) )
  end do
End Function

Function VecTrim_2D( Input_String ) result(Output_String)
  character(:)  ,dimension( :, : )      ,allocatable    ,intent(in)     ::      Input_String
  character(:)  ,dimension( :, : )      ,allocatable                    ::      Output_String
  integer                                                               ::      i, j
  integer                                                               ::      Length
  integer                                                               ::      Size1
  integer                                                               ::      Size2
  Length        =       Max_Len_Trim(Input_String)
  Size1         =       size(Input_String,1)
  Size2         =       size(Input_String,2)
  allocate( character(Length) :: Output_String(Size1,Size2) )
  do j = 1,Size2
  do i = 1,Size1
    Output_String(i,j)  =       trim( Input_String(i,j) )
  end do
  end do
End Function

! Function VecTrim_3D( Input_String ) result(Output_String)
!   character(:)  ,dimension(:,:,:)       ,allocatable    ,intent(in)     ::      Input_String
!   character(:)  ,dimension(:,:,:)       ,allocatable                    ::      Output_String
!   integer                                                               ::      i, j, k
!   integer                                                               ::      Length
!   integer                                                               ::      Size1
!   integer                                                               ::      Size2
!   integer                                                               ::      Size3
!   Length        =       Max_Len_Trim(Input_String)
!   Size1         =       size(Input_String,1)
!   Size2         =       size(Input_String,2)
!   Size3         =       size(Input_String,3)
!   allocate( character(Length) :: Output_String(Size1,Size2) )
!   do k = 1,Size3
!   do j = 1,Size2
!   do i = 1,Size1
!     Output_String(i,j,k)        =       trim( Input_String(i,j,k) )
!   end do
!   end do
! End Function













Pure Function Replace_Character( String_Inp, CharOld, CharNew ) result(String_Out)
  character(*)                                          ,intent(in)     ::      String_Inp                         !< String variable to be modified
  character(*)                                          ,intent(in)     ::      CharOld                        !<
  character(*)                                          ,intent(in)     ::      CharNew                        !<
  character(:)  ,allocatable                                            ::      String_Out                         !<
  integer                                                               ::      iOldIni                               ! Character index
  integer                                                               ::      iOldFin
  integer                                                               ::      iNewIni
  integer                                                               ::      iNewFin
  character(:)  ,allocatable                                            ::      String_Tmp
  String_Out    =       ""
  String_Tmp    =       String_Inp                                                                              ! Loading the input string in the temporary variable
  do                                                                                                            ! Loop on all character for character replacement
    iOldIni     =       index(String_Tmp,CharOld)                                                               ! Setting index of first (initial) character to be replaced in the old string
    iOldFin     =       iOldIni - 1 + len(CharOld)                                                              ! Setting index of last (final) character to be replaced in the old string
    iNewIni     =       iOldIni                                                                                 ! Setting index of first (initial) character to be replaced in the new string
    iNewFin     =       iOldIni - 1 + len(CharNew)                                                              ! Setting index of last (final) character to be replaced in the new string
    if ( iOldIni /= 0 ) then
      String_Out        =       String_Out // String_Tmp(1:iOldIni-1) // CharNew
      String_Tmp        =       String_Tmp(iOldFin+1:)
    else
      String_Out        =       String_Out // String_Tmp
      exit                                                                                    ! If character index is zero, then no character to be replace and so exiting the loop
    end if
  end do                                                                                                        ! End of loop on string characters
End Function

Pure Function Replace_Characters( String_Inp, CharOld, CharNew ) result(String_Out)
  character(*)                                          ,intent(in)     ::      String_Inp                      !< String variable to be modified
  character(*)  ,dimension(:)                           ,intent(in)     ::      CharOld                         !<
  character(*)  ,dimension(:)                           ,intent(in)     ::      CharNew                         !<
  character(:)  ,allocatable                                            ::      String_Out                      !<
  integer                                                               ::      i
  String_Out    =       String_Inp
  do i = 1,size(CharOld)
    String_Out  =       Replace_Character( String_Out, CharOld(i), CharNew(i) )
  end do
  String_Out    =       trim(String_Out)
End Function

! This procedure is used to enhanced species names so they have a nice rendering in LaTeX.
! Only works in a limited number of case:
! 1) the number must be lower than 9 (only one digits) (always the case for the considered species)
! 2) the fist character should not be a number (never the case for species name)
! 3) only work for one substitution (because of the allocate statement))
! However, that's fine for my needs...
Function Set_Enhanced_Name(Name) result(Enhanced_Name)
  character(*)                                          ,intent(in)     ::      Name                            !< Input name
  character(:)  ,allocatable                                            ::      Enhanced_Name                   !< Output name
  integer                                                               ::      i                               ! Index of character position
  character(:)  ,allocatable                                            ::      String_Old
  character(:)  ,allocatable                                            ::      String_New
  character(:)  ,allocatable                                            ::      String_L, String_R, String_C
  Enhanced_Name                    =       Name
  do i = 1,len(Enhanced_Name)
    String_Old                  =       Enhanced_Name(i:i)
    if ( .not. Is_Numeric(String_Old) ) cycle
    String_New                =       "_{" // String_Old // "}"
    String_L          =       Enhanced_Name(1:i-1)     ! String on the LHS of the substritution
    String_C          =       String_New
    String_R          =       Enhanced_Name(i+1:)      ! String on the RHS of the substritution
    Enhanced_Name     =       String_L // String_C // String_R
  end do
  do i = 1,len(Enhanced_Name)
    String_Old                  =       Enhanced_Name(i:i)
    if ( (String_Old/='+') .and. (String_Old/='-') ) cycle
    String_New                =       "^{" // String_Old // "}"
    String_L          =       Enhanced_Name(1:i-1)
    String_C          =       String_New
    String_R          =       Enhanced_Name(i+1:)
    Enhanced_Name     =       String_L // String_C // String_R
  end do
  Enhanced_Name         =       trim(Enhanced_Name)
End Function

Pure Function Is_Numeric( String )
  logical                                                       ::      Is_Numeric                              !< Indicator of numeric value stored in the input string
  character(*)                                  ,intent(in)     ::      String                                  !< String to be checked for validity as a numeric
  character(*)          ,parameter                              ::      Valid_Charaters=' 1234567890.'          ! List of characters wich are valid numbers
  Is_Numeric    =       ( verify(String,Valid_Charaters) == 0 )                                                 ! Setting the numeric indicator ( verify returns 0 if each character in String appears in Valid_Charaters)
End Function

Pure Function Is_A_Number( String )
  character(*)                                  ,intent(in)     ::      String                                  !< String to be checked for validity as a numeric
  logical                                                       ::      Is_A_Number                              !< Indicator of numeric value stored in the input string
  real(8)                                                       ::      Number
  integer                                                       ::      ios
  read(String,*,iostat=ios) Number                                                                                         ! Converting the input string into a real number
  Is_A_Number    =  ( ios == 0 )
End Function

Pure Function Remove_Last_Directory_From_Path( DirInp ) result(DirOut)
  character(*)                                          ,intent(in)     ::      DirInp
  character(:)  ,allocatable                                            ::      DirOut
  integer                                                               ::      Length
  integer                                                               ::      iSlash
  character(*)  ,parameter                                              ::      Slash_Char = '/'
  character(1)                                                          ::      Last_Char
  logical                                                               ::      Trailing_Slash
  DirOut        =       trim(DirInp)
  Length        =       len_trim(DirOut)
  Last_Char     =       DirOut(Length:Length)
  Trailing_Slash        =       ( Last_Char == Slash_Char )
  if (Trailing_Slash) DirOut = DirOut(1:Length-1)       ! Remove the trailing slash
  iSlash        =       index(DirOut,Slash_Char,back=.true.)
  if (Trailing_Slash) then
    DirOut      =       DirOut(1:iSlash)
  else
    DirOut      =       DirOut(1:iSlash-1)
  end if
End Function


! **************************************************************************************************************
! **************************************************************************************************************
!       PROCEDURES FOR INLINING AN ARRAY OF CHARACTER STRINGS INTO A SCALAR CHARACTER STRING
! **************************************************************************************************************
! **************************************************************************************************************

Pure Function Inline_Strings_1D( InpStr, Separator, Trimed, Invert, Fmt ) result(OutStr)

  character(*)  ,dimension(:)                           ,intent(in)     ::      InpStr                          !< Array of strings to be inlined
  character(*)  ,optional                               ,intent(in)     ::      Separator                       !< Separator used between elements
  character(:)  ,allocatable                                            ::      OutStr                          !< Scalar character string corresponding to the inlined input string
  character(*)  ,optional                               ,intent(in)     ::      Fmt                             !< Format used to describe the number
  logical       ,optional                               ,intent(in)     ::      Trimed                          !< Indicator whether the strings should be trimed
  logical       ,optional                               ,intent(in)     ::      Invert                          !< Indicator whether the strings elements should be inverted (NOT IMPLEMENTED)\

  integer                                                               ::      i                               ! Indexes of the input character string elements
  character(:)  ,allocatable                                            ::      String
  character(:)  ,allocatable                                            ::      Separator_Local                 ! Local separator used between elements
  logical                                                               ::      Trimed_Local                    ! Local indicator whether the strings should be trimed

  Separator_Local       =       " "                                                                             ! Initializing the separator to a single space character
  if ( present(Separator) ) Separator_Local = Separator                                                         ! Setting the separator to the input value if present

  Trimed_Local          =       .true.                                                                          ! Initializing the trimed indicator to true
  if ( present(Trimed) ) Trimed_Local = Trimed                                                                  ! Setting the trimed indicator to the input value if present

  OutStr        =       ""                                                                                      ! Initializing the ouput character to an empty string (required because of iterative computation)
  do i = 1,size(InpStr)                                                                                         ! Loop on all elements of the input string
    if ( present(Fmt) ) then
      String    =       Convert_To_String( InpStr(i), Fmt=Fmt )
    else
      String    =       InpStr(i)                                                                               ! Getting the current element
    end if
    if (Trimed_Local) String = trim(String)                                                                     ! Removing trailling blanks if required from optional input arguments (Remove them by default)
    OutStr      =       OutStr // String                                                                  ! Adding the current element to the outpout character
    if ( i /= size(InpStr) ) OutStr      =       trim(OutStr) // Separator_Local                                ! Adding the separation character if except for the last element
  end do                                                                                                        ! End loop on elements of the input string

End Function

! This procedure returns a scalar character string corresponding to the inlined input vector of reals
! Pure
Function Inline_Reals_1D( InpVar, Separator, Trimed, Invert, Fmt ) result(OutStr)
  real(8)       ,dimension(:)                           ,intent(in)     ::      InpVar                          !< Array of real variable to be converted to character and inlined
  character(*)  ,optional                               ,intent(in)     ::      Separator                       !< Separator used between elements
  logical       ,optional                               ,intent(in)     ::      Trimed                          !< Indicator whether the strings should be trimed
  logical       ,optional                               ,intent(in)     ::      Invert                          !< Indicator whether the strings elements should be inverted (NOT IMPLEMENTED)
  character(*)  ,optional                               ,intent(in)     ::      Fmt                             !< Format used to describe the number
  character(:)  ,allocatable                                            ::      OutStr                          !< Scalar character string corresponding to the inlined input string
  character(:)  ,dimension(:)   ,allocatable                            ::      InpStr                          !< Array of strings to be inlined
  integer                                                               ::      i                               ! Indexes of the input character string elements
  character(:)  ,allocatable                                            ::      String
  do i = 1,size(InpVar)
    String      =       Convert_To_String( InpVar(i), Fmt=Fmt )
    call Add_Line_To_String( InpStr, String )
  end do
  OutStr        =       Inline( InpStr, Separator, Trimed, Invert )
End Function

! This procedure returns a scalar character string corresponding to the inlined input vector of integers
Pure Function Inline_Integers_1D( InpVar, Separator, Trimed, Invert, Fmt ) result(OutStr)
  integer       ,dimension(:)                           ,intent(in)     ::      InpVar                          !< Array of real variable to be converted to character and inlined
  character(*)  ,optional                               ,intent(in)     ::      Separator                       !< Separator used between elements
  logical       ,optional                               ,intent(in)     ::      Trimed                          !< Indicator whether the strings should be trimed
  logical       ,optional                               ,intent(in)     ::      Invert                          !< Indicator whether the strings elements should be inverted (NOT IMPLEMENTED)
  character(*)  ,optional                               ,intent(in)     ::      Fmt                             !< Format used to describe the number
  character(:)  ,allocatable                                            ::      OutStr                          !< Scalar character string corresponding to the inlined input string
  character(:)  ,dimension(:)   ,allocatable                            ::      InpStr                          !< Array of strings to be inlined
  integer                                                               ::      i                               ! Indexes of the input character string elements
  character(:)  ,allocatable                                            ::      String
  do i = 1,size(InpVar)
    String      =       Convert_To_String( InpVar(i), Fmt=Fmt )
    call Add_Line_To_String( InpStr, String )
  end do
  OutStr        =       Inline( InpStr, Separator, Trimed, Invert )
End Function



Function Get_SubString( InpStr, InBetween ) result(OutStr)
  character(*)                                          ,intent(in)     ::      InpStr
  character(*)  ,dimension(2)                           ,intent(in)     ::      InBetween
  character(:)  ,allocatable                                            ::      OutStr
  character(:)  ,allocatable                                            ::      Str1
  character(:)  ,allocatable                                            ::      Str2
  integer                                                               ::      iIni
  integer                                                               ::      iFin
  Str1          =       trim(InBetween(1))
  Str2          =       trim(InBetween(2))
  OutStr         =       ""
  iIni          =       index( InpStr, Str1 )
  if ( iIni == 0 ) return
  iIni          =       iIni + len(Str1)
  iFin          =       index( InpStr, Str2 )
  if ( iFin == 0 ) return
  iFin          =       iFin - 1
  OutStr        =       InpStr(iIni:iFin)
End Function

Function Get_SubString_Indexes( InpStr, InBetween ) result(Indexes)
  character(*)                                          ,intent(in)     ::      InpStr
  character(*)  ,dimension(2)                           ,intent(in)     ::      InBetween
  integer       ,dimension(2)                                           ::      Indexes
  character(:)  ,allocatable                                            ::      Str1
  character(:)  ,allocatable                                            ::      Str2
  integer                                                               ::      iIni
  integer                                                               ::      iFin
  Str1          =       trim(InBetween(1))
  Str2          =       trim(InBetween(2))
  iIni          =       index( InpStr, Str1 )
  if ( iIni == 0 ) return
  iIni          =       iIni + len(Str1)
  iFin          =       index( InpStr, Str2, back=.true. )
  if ( iFin == 0 ) return
  iFin          =       iFin - 1
  Indexes(1)    =       iIni
  Indexes(2)    =       iFin
End Function

Function Add_Prefix_0d( String, Prefix ) result(Prefixed_String)
  character(*)                                          ,intent(in)     ::      String
  character(*)                                          ,intent(in)     ::      Prefix
  character(len(Prefix)+len(String))                                    ::      Prefixed_String
  Prefixed_String       =       Prefix // String
End Function

Function Add_Prefix_1d( Strings, Prefix ) result(Prefixed_Strings)
  character(*)  ,dimension(:)                           ,intent(in)     ::      Strings
  character(*)                                          ,intent(in)     ::      Prefix
  character(len(Prefix)+len(Strings))   ,dimension(size(Strings))       ::      Prefixed_Strings
  integer                                                               ::      i
  do i = 1,size(Strings)
    Prefixed_Strings(i)   =       Prefix // Strings(i)
  end do
End Function

Function Add_Suffix_0d( String, Suffix ) result(Suffixed_String)
  character(*)                                          ,intent(in)     ::      String
  character(*)                                          ,intent(in)     ::      Suffix
  character(len(Suffix)+len(String))                                    ::      Suffixed_String
  Suffixed_String       =       trim(String) // Suffix
End Function

Function Add_Suffix_1d( Strings, Suffix ) result(Suffixed_Strings)
  character(*)  ,dimension(:)                           ,intent(in)     ::      Strings
  character(*)                                          ,intent(in)     ::      Suffix
  character(len(Suffix)+len(Strings))   ,dimension(size(Strings))       ::      Suffixed_Strings
  integer                                                               ::      i
  do i = 1,size(Strings)
    Suffixed_Strings(i)   =       trim(Strings(i)) // Suffix
  end do
End Function

Pure Function Convert_Ratio( Numerator, Denominator ) result(Ratio)
  integer                                               ,intent(in)     ::      Numerator
  integer                                               ,intent(in)     ::      Denominator
  character(:)  ,allocatable                                            ::      Ratio
  integer                                                               ::      NDigits
  character(:)  ,allocatable                                            ::      NumStr
  character(:)  ,allocatable                                            ::      DenStr
  NDigits       =       Get_Number_Of_Digits( Denominator )                                                     ! Getting the number of digits of the denominator
  NumStr        =       Convert_To_String( Numerator, Pos='R', Len=NDigits )                                    ! Converting the numerator to a string
  DenStr        =       Convert_To_String( Denominator )                                                        ! Converting the numerator to a string
  Ratio         =       NumStr // "/" // DenStr                                                                 ! Setting the string corresponding to the 2 input numbers
End Function


Pure Function Remove_Character( String_Inp, CharOld ) result(String_Out)
  character(*)                                          ,intent(in)     ::      String_Inp                         !< String variable to be modified
  character(*)                                          ,intent(in)     ::      CharOld                        !<
  character(:)  ,allocatable                                            ::      String_Out                         !<
  character(*)  ,parameter                                              ::      CharNew = ''
  String_Out    =       Replace_Character( String_Inp, CharOld, CharNew )                                        ! End of loop on string characters
End Function

Pure Function Remove_Characters( String_Inp, CharOld ) result(String_Out)
  character(*)                                          ,intent(in)     ::      String_Inp                         !< String variable to be modified
  character(*)  ,dimension(:)                           ,intent(in)     ::      CharOld                        !<
  character(:)  ,allocatable                                            ::      String_Out                         !<
  integer                                                               ::      i
  character(*)  ,parameter                                              ::      CharNew = ''
  String_Out    =     String_Inp
  do i = 1,size(CharOld)
    String_Out  =     Replace_Character( String_Out, CharOld(i), CharNew )
  end do
End Function

Subroutine Set_Frame( Text, Frame )

  character(*)                                          ,intent(in)     ::    Text                            !< Character string around which a frame is to be created
  character(:)  ,allocatable    ,dimension(:)           ,intent(out)    ::    Frame                           !< Character string representing the framed text

  integer                                                               ::    i
  integer                                                               ::    MaxLen                          ! Maximum length of the output character string
  integer                                                               ::    NLines                          ! Number of lines of the frame including the text
  integer       ,parameter                                              ::    NSpace=2                        ! Number of character used for spacing before and after the text in the frame
  integer       ,parameter                                              ::    NTab=2                      ! Number of character used for spacing before and after the text in the frame
  character(*)  ,parameter                                              ::    Char_Frame='#'                  ! Character used for the frame
  character(3)                                                          ::    Char_MaxLen                     ! String representing the number of characters in the output string
  character(3)                                                          ::    Char_NSpace                     ! String representing the number of blanc characters before and after the text message
  character(3)                                                          ::    Char_NTab
  character(:)  ,allocatable                                            ::    String                          ! Text to be framed (same tha input text but with no blanks at the begining and at the end)
  character(:)  ,allocatable                                            ::    Fmt_Ext                         ! Format for the external frame, ie. the upper and lower lines of the frame
  character(:)  ,allocatable                                            ::    Fmt_Int                         ! Format for the internal frame, ie. the text inside the frame

  String        =       trim(adjustl(Text))
  NLines        =       3                                                                                       ! Setting the number of lines to 3
  MaxLen        =       len(String) + 2 * (NSpace + 1)                                                         ! Number of characters in the output string including blanks and the framing characters
  allocate( character(MaxLen) :: Frame(NLines) )                                                                ! Allocation of dimension and length of the output string
  write(Char_MaxLen,"(i3)") MaxLen                                                                              ! Interger-to-string convertion of the length of the output string
  write(Char_NSpace,"(i3)") NSpace                                                                              ! Interger-to-string convertion of the number of blanks characters
  write(Char_NTab,"(i3)") NTab

!   Fmt_Ext       =       "(" // Char_NTab // "x" // Char_MaxLen // "('"// Char_Frame // "'))"                                        ! Setting the format for the external frame, ie. the upper and lower lines of the frame
!   Fmt_Int       =       "(" // Char_NTab // "x,'"//Char_Frame//"',"//Char_NSpace//"x,'"//String//"',"//Char_NSpace//"x,'"//Char_Frame//"')"  ! Setting the format for the internal frame, ie. the text inside the frame

  Fmt_Ext       =       "(" // Char_MaxLen // "('"// Char_Frame // "'))"                                        ! Setting the format for the external frame, ie. the upper and lower lines of the frame
  Fmt_Int       =       "('"//Char_Frame//"',"//Char_NSpace//"x,'"//String//"',"//Char_NSpace//"x,'"//Char_Frame//"')"  ! Setting the format for the internal frame, ie. the text inside the frame

  write(Frame(1),Fmt_Ext)                                                                                       ! Writing the upper line of the frame
  do i = 2,Nlines-1                                                                                             ! Loop on all internal lines
    write(Frame(i),Fmt_Int)                                                                                     ! Writing the middle line of the frame
  end do                                                                                                        ! End loop on
  write(Frame(Nlines),Fmt_Ext)                                                                                  ! Writing the lower line of the frame

End Subroutine

Pure Function Equal( S1, S2, Trimed, CaseSensitive )

  character(*)                                          ,intent(in)     ::      S1
  character(*)                                          ,intent(in)     ::      S2
  logical               ,optional                       ,intent(in)     ::      Trimed                          !< Indicator whether
  logical               ,optional                       ,intent(in)     ::      CaseSensitive                  !< Indicator whether the search should be case sensitive
  logical                                                               ::      Equal

  logical                                                               ::      CaseSensitive_
  logical                                                               ::      i_Trimed
  character(:)  ,allocatable                                            ::      String_1
  character(:)  ,allocatable                                            ::      String_2

  String_1    =     S1
  String_2    =     S2

  CaseSensitive_ =       .true.
  if ( present(CaseSensitive) ) CaseSensitive_ = CaseSensitive

  i_Trimed      =       .true.
  if ( present(Trimed) ) i_Trimed = Trimed

  if ( .Not.CaseSensitive_ ) then
    String_1    =     UpperCase(String_1)
    String_2    =     UpperCase(String_2)
  end if

  if ( i_Trimed ) then
    String_1    =     trim(String_1)
    String_2    =     trim(String_2)
  end if

  Equal         =       String_1 == String_2

End Function


! **************************************************************************************************************
! **************************************************************************************************************
!                             PRIVATE PROCEDURES
! **************************************************************************************************************
! **************************************************************************************************************

Subroutine Error_Message( ErrMsg, ProcName )
  character(*)  ,optional                               ,intent(in)     ::      ErrMsg
  character(*)  ,optional                               ,intent(in)     ::      ProcName
  write(*,"('Fatal error leading to code termination')")
  if ( present(ProcName) ) write(LogUnit,"(4x,'Calling procedure   :   ',g0)") trim(ProcName)                   ! Printing the calling procedure
  if ( present(ErrMsg)   ) write(LogUnit,"(4x,'Error message       :   ',g0)") trim(ErrMsg)                     ! Printing the error message
  stop                                                                                                          ! Stopping the code
End Subroutine

! **************************************************************************************************************
! **************************************************************************************************************
!                                       WORK-IN-PROGRESS
! **************************************************************************************************************
! **************************************************************************************************************

Subroutine Parse_Fixed_Size( String_Input, Separator, String_Output, EscRHS, Ignore_Between )
  character(*)                                          ,intent(in)     ::      String_Input                    !< Input scalar character string to be proceeded
  character(*)                                          ,intent(in)     ::      Separator                       !< Separation character string
  character(*)  ,dimension(:)                           ,intent(out)    ::      String_Output                   !< Output vector character string containing sub-strings
  character(1)  ,optional       ,dimension(:)           ,intent(in)     ::      EscRHS                          !< RHS escape character
  character(*)  ,optional       ,dimension(:)           ,intent(in)     ::      Ignore_Between
  character(len=len(String_Input))                                      ::      String                    ! Local copy of the input string
  character(len=len(String_Input))                                      ::      String_LHS
  character(len=len(String_Input))                                      ::      String_RHS
  integer                                                               ::      i
  String_Output =       ""
  String        =       Compact(String_Input)                                                                   ! Storing local value of the input character string and compacting it
  if ( len_trim(String) == 0 ) return                                                                           ! If empty input string, then exiting the procedure
  i          =       0
  do                                                                                                            ! Infinite loop on all characters of the input string
    i        =       i + 1
    if ( i > size(String_Output) ) exit
    if ( len_trim(String)==0 ) exit                                                                             ! If the current local character string is empty, then exiting the procedure
    call Split_Fixed_Size( String, String_LHS, Separator, String_RHS, EscRHS=EscRHS, Ignore_Between=Ignore_Between )       ! Spliting the string into a left and right string separated by the Separator character
    String              =       String_RHS                                                                      ! Setting the remaining string to be proceed to the RHS string                                                                                                     ! If allocated output String_Output variable, then sub-strings have already been proceeded
    String_Output(i)    =       String_LHS                                                                      !
  end do                                                                                                        ! End of loop
End Subroutine

Subroutine Split_Fixed_Size( Input_String, LHS, Separator, RHS, EscLHS, EscRHS, i_Length_Eq_Inp, Ignore_Between, i_Debug )

  character(*)                                          ,intent(in)     ::      Input_String                          !< Input character string to be splitted into a LHS and a RHS parts
  character(*)                                          ,intent(out)    ::      LHS                             !< Left-hand-side part of the splitted string
  character(*)                                          ,intent(in)     ::      Separator                       !< Separation character string
  character(*)                                          ,intent(out)    ::      RHS                             !< Right-hand-side part of the splitted string
  character(1)  ,optional                               ,intent(in)     ::      EscLHS                          !< Left-hand-side escape character string
  character(1)  ,optional       ,dimension(:)           ,intent(in)     ::      EscRHS                          !< Right-hand-side escape character string
  logical       ,optional                               ,intent(in)     ::      i_Length_Eq_Inp                 ! Indicator to force the length of the output string to equal the one of the input string (required because otherwise the error:: fort: (4): Variable STRLHS has substring ending point 1 which is greater than the variable length of 0) Mayu be possible to
  character(*)  ,optional       ,dimension(:)           ,intent(in)     ::      Ignore_Between
  logical       ,optional                               ,intent(in)     ::      i_Debug

  logical                                                               ::      i_Debug_Loc
  logical                                                               ::      i_EscLHS                        ! Escape character indicator
  logical                                                               ::      i_EscRHS                        ! Escape character indicator
  logical                                                               ::      i_protect                       ! Character protection indicator
  logical                                                               ::      i_Length_Eq_Inp_Local
  logical                                                               ::      i_Ignore_Between
!   character(:)  ,allocatable    ,dimension(:)                           ::      Char_Ignore_Between
  integer                                                               ::      iStrInp                         ! Index of the input string
  integer                                                               ::      iRHS                            ! Index of the RHS string
  integer                                                               ::      iLHS                            ! Index of the LHS string
  integer                                                               ::      iSepIni                         ! Initial index of the separation character in the input string
  integer                                                               ::      iSepFin                         ! Final index of the separation character in the input string
  character(1)                                                          ::      Char1                           ! String single character
  character(1)                                                          ::      Char1_Next                      ! Next character
  character(len=len_trim(Input_String))                                 ::      StrInp                          ! Local copy of the input string
  integer                                                               ::      iEscRHS                         ! Index of element in the array of RHS escape character

  character(:)  ,allocatable                                            ::      String
  character(:)  ,allocatable                                            ::      Ignore_CharIni                  ! Initial character from which the separation character should be ignored
  character(:)  ,allocatable                                            ::      Ignore_CharFin                  ! final character   to   which the separation character should be ignored

  integer                                                               ::      Position_Ignore_CharIni
  integer                                                               ::      Position_Ignore_CharFin
  integer                                                               ::      i

  i_Debug_Loc = .false.; if ( present(i_Debug) ) i_Debug_Loc = i_Debug

  RHS   =       Input_String

! ==============================================================================================================
!       PROCESSING OPTIONAL INPUT ARGUMENTS
! ==============================================================================================================
  i_EscRHS              =       .false.                                                                         ! Initialization of the RHS escape character indicator
  i_EscLHS              =       .false.                                                                         ! Initialization of the LHS escape character indicator
  i_Length_Eq_Inp_Local =       .false.
  i_Ignore_Between      =       .false.
  if ( present(EscRHS) )                i_EscRHS = .true.                                                       ! If a RHS escape character is provided, the activating the RHS escape indicator
  if ( present(EscLHS) )                i_EscLHS = .true.                                                       ! If a LHS escape character is provided, the activating the LHS escape indicator
  if ( present(i_Length_Eq_Inp) )       i_Length_Eq_Inp_Local = i_Length_Eq_Inp
  if ( present(Ignore_Between) )        i_Ignore_Between = .true.
! ==============================================================================================================


  if ( i_EscRHS ) RHS = RemoveSpace( adjustl(RHS) )                                                          ! Removing blanks and adjusting the to LHS
  LHS       =       ''                                                                                   ! Initialization of the LHS character string

  iSepIni   =       index( RHS, Separator )                                                                 ! Index of initial position of the separation character in the input character
  iSepFin   =       iSepIni + len_trim(Separator) - 1                                                      ! Index of final   position of the separation character in the input character


! ==============================================================================================================
  if ( i_Ignore_Between ) then
    if ( size(Ignore_Between) == 1 ) then
      Ignore_CharIni    =       Ignore_Between(1)
      Ignore_CharFin    =       Ignore_Between(1)
    else if ( size(Ignore_Between) >= 2 ) then
      Ignore_CharIni    =       Ignore_Between(1)
      Ignore_CharFin    =       Ignore_Between(2)
    end if
    Position_Ignore_CharIni     =       0   !     This only works for ignore characters of 1 character
    Position_Ignore_CharFin     =       0
    do i = 1,len_trim(RHS)
      Char1     =       RHS(i:i)
      if ( (Position_Ignore_CharIni == 0) .and.  (Char1 == Ignore_CharIni) ) Position_Ignore_CharIni = i
      if ( (Position_Ignore_CharFin == 0) .and.  (Char1 == Ignore_CharFin) .and. (Position_Ignore_CharIni/=i) ) Position_Ignore_CharFin = i
    end do
    if ( ( Position_Ignore_CharIni < iSepIni ) .and. ( iSepIni < Position_Ignore_CharFin ) ) then
      String            =       RHS(Position_Ignore_CharFin+1:)
      iSepIni   =       index( String, Separator )                                                                  ! Index of initial position of the separation character in the input character
      iSepFin   =       iSepIni + len_trim(Separator) - 1                                                      ! Index of final   position of the separation character in the input character
      if ( iSepIni /= 0 ) then                                        ! if iSepIni == 0 , then the separation character is not present in the string (outside the ignore character). So the zero value is kept so that the procedure is exited below.
        iSepIni   =       iSepIni + (Position_Ignore_CharFin)
        iSepFin   =       iSepIni + len_trim(Separator) - 1                                                      ! Index of final   position of the separation character in the input character
      end if
    end if
  end if
! ==============================================================================================================


! ==============================================================================================================
!       TREATING CASES FOR WHICH THE INPUT STRING DOES NOT NEED ANY SPLITTING
! ==============================================================================================================
! This section deals with the cases when no splitting is required. There are 3 situations for which the input
! string does not need any splitting. For all these cases the input string is stored in the output LHS string
! variable 'LHS' and the output RHS string variable 'RHS' is set to an empty string.
! This ensure the "Parse" calling procedure to stop its splitting iteration, the iteration being stopped when
! RHS string variable 'RHS' has a zero length. The 3 situations are:
!  - When the input character correspond to an empty string (obvious)
!  - When the separation character is absent from the input string
!  - When the separation character corresponds to the last character of the input string and that the
!    separation character is an RHS escape character (yes, it's a tricky one !)
! ==============================================================================================================
  if ( len_trim(RHS) == 0 ) return                                                                                ! If the input string is empty, then exiting the procedure without changing the input string
! --------------------------------------------------------------------------------------------------------------
  if ( iSepIni == 0 ) then                                                                                      ! If the separation character is absent from the input string, then
    LHS      =       RHS                                                                                        ! Storing the input string in the LHS string
    RHS      =       ''                                                                                         ! Setting the RHS string to an empty string
    return                                                                                                      ! Exiting the procedure
  end if                                                                                                        ! End of if case on initial separation character index
! --------------------------------------------------------------------------------------------------------------
  if ( (i_EscRHS) .and. (iSepIni == len_trim(RHS)) )  then                                                      ! If a RHS escape character exists and if the last character of the input string corresponds to the separation character
  Char1 =       RHS(iSepIni:iSepIni)                                                                            ! Storing the last character
    do iEscRHS = 1,size(EscRHS)                                                                                 ! Loop on all RHS escape character
      if        (Char1 /= EscRHS(iEscRHS))      cycle                                                           ! Going to the next RHS escape character if the last character from the input string does not correspond to the current RHS escape character
        LHS  =       RHS                                                                                        ! Storing the input string in the LHS string
        RHS  =       ''                                                                                         ! Setting the RHS string to an empty string
        return                                                                                                  ! Exiting the procedure
    end do                                                                                                      ! End loop on RHS escape character
  end if                                                                                                        ! End of if case on initial separation character index
! ==============================================================================================================





  StrInp        =       RHS                                                                                     ! Storing the input string in a new variable
  RHS(:)        =       ''                                                                                      ! Initialization of RHS string
  iLHS          =       0                                                                                       ! Initialization of LHS string index
  iRHS          =       0                                                                                       ! Initialization of RHS string index
  iStrInp       =       0                                                                                       ! Initialization of input string index
  i_protect     =       .false.                                                                                 ! Initialization of the character protection indicator
  do                                                                                                            ! Infinit loop on all characters of the input string
    iStrInp     =       iStrInp + 1                                                                             ! Incrementing the input string index
    if  (iStrInp > len_trim(StrInp))    exit                                                                    ! Exiting the loop if the the index is greater than the string length
    Char1       =       StrInp(iStrInp:iStrInp)                                                                 ! Getting the current character
!   Case of a LHS escape character
!   ------------------------------
    if (i_EscLHS) then                                                                                          ! If an escape character exists, then
      if (i_protect) then                                                                                       ! If the current character is protected, then
        i_protect       =       .false.                                                                         ! Setting off the character protection mode
        iLHS            =       iLHS + 1                                                                        ! Index incrementation
        LHS(iLHS:iLHS)  =       Char1                                                                           ! Storing the current character in the LHS string
        cycle                                                                                                   ! Going to the next character
      end if                                                                                                    ! End of if case on character protection
      if ( Char1 == EscLHS ) then                                                                               ! If the current character corresponds to the escape character, then
        iLHS            =       iLHS + 1                                                                        ! Index incrementation
        LHS(iLHS:iLHS)  =       Char1                                                                           ! Storing the current character in the LHS string
        i_protect       =       .true.                                                                          ! Activating the character protection mode
        cycle                                                                                                   ! Going to the next character
      end if                                                                                                    ! End of if case on escape character
    end if                                                                                                      ! End of if case on existence of an escape character

!   Case of a RHS escape character
!   ------------------------------
    if (i_EscRHS) then                                                                                          ! If an RHS escape character exists, then
      if ( iStrInp == iSepIni ) then                                                                            ! If the current input string index correspond to the initial separation index
        Char1_Next      =       StrInp(iStrInp+1:iStrInp+1)                                                     ! Storing the value of the next character
        do iEscRHS = 1,size(EscRHS)                                                                             ! Loop on all RHS escape characters
          if    (Char1_Next == EscRHS(iEscRHS)) then                                                            ! If the next character corresponds to the RHS escape character, then
             iSepIni        =       iSepIni + index( StrInp(iStrInp+1:), Separator )                            ! Modyfiying the index of initial position of the separation character in the input character: Moving to the next separation character index
             iSepFin        =       iSepIni + len_trim(Separator) - 1                                           ! Modyfiying the index of final   position of the separation character in the input character
          end if                                                                                                ! End of if case on escape character
        end do                                                                                                  ! End loop on RHS escape characters
      end if                                                                                                    ! End if case on input string index
    end if                                                                                                      ! End of if case on RHS escape indicator

    if ( iStrInp < iSepIni ) then                                                                               ! If the current input string index is lower than the initial separation index
      iLHS              =       iLHS + 1                                                                        ! Index incrementation
      LHS(iLHS:iLHS)    =       Char1                                                                           ! Storing the current character in the LHS string
      cycle                                                                                                     ! Going to the next character
    end if                                                                                                      ! End if case on input string index

    if (iStrInp > iSepFin) then                                                                                 ! If the current input string index is greater than the final separation index
      iRHS              =       iRHS + 1                                                                        ! Index incrementation
      RHS(iRHS:iRHS)    =       Char1                                                                           ! Storing the current character in the LHS string
      cycle                                                                                                     ! Going to the next character
    end if                                                                                                      ! End if case on input string index

  end do                                                                                                        ! End loop on character of the input string

  RHS           =       adjustl(RHS)                                                                         ! Removing initial spaces
  RHS           =       trim(RHS)
  LHS           =       trim(LHS)

End Subroutine

End Module