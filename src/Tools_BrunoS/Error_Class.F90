! This file is intended to be replace by the 'true' Logger_Class module
Module Error_Class

  use Logger_Class    ,only:  Logger

  implicit none

  private
  public    ::    Error
  public    ::    CheckVariable

  integer ,parameter     ::    rkp = 8

  Interface           CheckVariable
    Module Procedure  CheckVariable_0d
    Module Procedure  CheckVariable_1d
    Module Procedure  CheckVariable_2d
  End Interface

  contains

Subroutine Error( Message )
  character(*)          ,intent(in)     ::    Message
  write(Logger%Unit,"('***** ERROR *****')")
  write(Logger%Unit,"(g0)") Message
  write(Logger%Unit,"('***** ERROR *****')")
  stop
End Subroutine

Subroutine CheckVariable_0d( Variable, ProcName, VarName )
  use ,intrinsic :: IEEE_Arithmetic ,only:  ieee_is_normal
  real(rkp)                 ,intent(in)   ::    Variable
  character(*)    ,optional ,intent(in)   ::    ProcName
  character(*)    ,optional ,intent(in)   ::    VarName
  if ( ieee_is_normal( Variable ) ) return
  write(Logger%Unit,"('NaN found in scalar variable')")
  if ( present(ProcName) ) write(Logger%Unit,"(' * Procedure: ',g0)") ProcName
  if ( present(VarName ) ) write(Logger%Unit,"(' * Variable:  ',g0)") VarName
  error stop
End Subroutine

Subroutine CheckVariable_1d( Variable, ProcName, VarName )
  use ,intrinsic :: IEEE_Arithmetic ,only:  ieee_is_normal
  real(rkp) ,dimension(:)   ,intent(in)   ::    Variable
  character(*)    ,optional ,intent(in)   ::    ProcName
  character(*)    ,optional ,intent(in)   ::    VarName
  integer                                 ::    i
  do i = 1,size(Variable,1)
    if ( ieee_is_normal( Variable(i) ) ) cycle
    write(Logger%Unit,"('NaN found in rank-1 variable')")
    write(Logger%Unit,"(' * Element: A(i) with i = ',g0)") i
  if ( present(ProcName) ) write(Logger%Unit,"(' * Procedure: ',g0)") ProcName
  if ( present(VarName ) ) write(Logger%Unit,"(' * Variable:  ',g0)") VarName
    error stop
  end do
End Subroutine

Subroutine CheckVariable_2d( Variable, ProcName, VarName )
  use ,intrinsic :: IEEE_Arithmetic ,only:  ieee_is_normal
  real(rkp) ,dimension(:,:) ,intent(in)   ::    Variable
  character(*)    ,optional ,intent(in)   ::    ProcName
  character(*)    ,optional ,intent(in)   ::    VarName
  integer                                 ::    i,j
  do j = 1,size(Variable,2)
  do i = 1,size(Variable,1)
    if ( ieee_is_normal( Variable(i,j) ) ) cycle
    write(Logger%Unit,"('NaN found in rank-2 variable')")
    write(Logger%Unit,"(' * Element: A(i,j) with i = ',g0,' j = ',g0)") i, j
  if ( present(ProcName) ) write(Logger%Unit,"(' * Procedure: ',g0)") ProcName
  if ( present(VarName ) ) write(Logger%Unit,"(' * Variable:  ',g0)") VarName
    error stop
  end do
  end do
End Subroutine

End Module
