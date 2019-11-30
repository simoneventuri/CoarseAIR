
  character(:)  ,allocatable                            ,intent(out)    ::  String
  character(*)                                ,optional ,intent(in)     ::  VarFmt
  character(*)                                ,optional ,intent(in)     ::  TypFormat
  character(*)                                ,optional ,intent(in)     ::  ComFmt
  integer                                     ,optional ,intent(out)    ::  Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)

  character(10000)                                                      ::  Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::  Local_Format
  integer                                                               ::  ios
  integer                                                               ::  i, iIni, iFin, NElements, NTotal
  character(:)  ,allocatable                                            ::  Separator
  integer                                                               ::  NumberOfRanges
  integer                                                               ::  Factor
  integer                                                               ::  Iter
  integer                                                   ,parameter  ::  NumberOfRangesMax  = 32
  integer                                                   ,parameter  ::  IterMax            = 10

  Local_Format    =       "(" // Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt ) // ")"

! ==============================================================================================================
!   FIRST CONVERTION ATTEMPT
! ==============================================================================================================
! First, lets try to convert directly the string.
! The 'write' wilol fail if 'Variable' contains too many elements (of the order of 1000)>
! If Variable contains
! ==============================================================================================================
  write( Long_String, Local_Format, iostat=ios ) Variable
  if ( ios == 0 ) then
    String = trim(Long_String)
    if ( present(Status) ) Status = ios
    return
  end if
! ==============================================================================================================


! ==============================================================================================================
!   SECOND CONVERTION ATTEMPT
! ==============================================================================================================
! If the conversion has failed, then we try to convert groups of elements one at the time.
!   if ( ios /= 0 ) then    ! output statement overflows record, unit -5, file Internal Formatted Write
! ==============================================================================================================
  String            =   ""
  NumberOfRanges    =   1
  Factor            =   2
  Iter              =   0
  NTotal            =   size(Variable)
  MainLoop: do
    Iter            =   Iter + 1
    NumberOfRanges  =   NumberOfRanges * Factor
    if ( NumberOfRanges >= NumberOfRangesMax ) exit      !       if ( Iter > IterMax ) exit
    iFin        =   0
    Separator   =   ""
    NElements   =   floor( real(NTotal) / NumberOfRanges ) + mod( NTotal, NumberOfRanges )
    do i = 1,NumberOfRanges
      iIni      =   1    + iFin
      iFin      =   iIni + NElements - 1
      if ( i == NumberOfRanges ) iFin = min(iFin,NTotal)
      write( Long_String, Local_Format, iostat=ios ) Variable(iIni:iFin)
      if ( ios /= 0 ) cycle MainLoop
      String    =   String // Separator // trim(Long_String)
      Separator =   '   '
      if ( iFin == size(Variable) ) exit MainLoop
    end do
  end do  MainLoop
  if ( ios == 0 ) then
    if ( present(Status) ) Status = ios
    return
  end if

! ==============================================================================================================
!   THIRD CONVERTION ATTEMPT
! ==============================================================================================================
  Local_Format    =   "(*(g0,1x))"
  write( Long_String, Local_Format, iostat=ios ) Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
! ==============================================================================================================
