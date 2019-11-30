Module Output_Module

  use Parameters_Module     ,only:  rkp
  use Logger_Class          ,only:  Logger

  implicit none

  private
  public    ::    OutputFile_Type

  Type      ::    OutputFile_Type
    character(:)  ,allocatable                ::    Name
    integer                                   ::    Unit    =   6
    integer                                   ::    irecl   =   0
    integer                                   ::    NBuffer =   0
    real(rkp)     ,allocatable  ,dimension(:) ::    Buffer
  contains
    procedure::   Open    =>    OpenFile
    procedure::   Write   =>    WriteVariable
  End Type

  logical   ,parameter    ::    i_Debug_Global = .False.
  logical   ,parameter    ::    Formatted = .False.
  integer   ,parameter    ::    isz      =     50000000

!   real(rkp) ,dimension(isz)    ::  GlobalStorage

  contains



Subroutine OpenFile( This, FileName, NEqtTot, Unit, i_Debug )

  class(OutputFile_Type)                          ,intent(out)    ::    This
  character(*)                              ,intent(in)     ::    FileName     ! Name of the file to be opened
  integer                                   ,intent(in)     ::    NEqtTot            ! Integer variable used to compute the record length
  integer                         ,optional ,intent(in)     ::    Unit
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    RecordLength
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "OpenFile")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()


  This%irecl      =     NEqtTot+2
  RecordLength    =     rkp * This%irecl !   * 4
  if (Formatted) RecordLength = RecordLength * 4
  This%Name       =     trim(FileName)
  if ( present(Unit) ) then
    This%Unit     =     Unit
    if (Formatted) then
      open( File=FileName, Unit=Unit,         Form='formatted', Access='direct', Status='unknown', Recl=RecordLength )
    else
      open( File=FileName, Unit=Unit,         Form='unformatted', Access='direct', Status='unknown', Recl=RecordLength )
    end if
  else
    if (Formatted) then
      open( File=FileName, NewUnit=This%Unit, Form='formatted', Access='direct', Status='unknown', Recl=RecordLength )
    else
      open( File=FileName, NewUnit=This%Unit, Form='unformatted', Access='direct', Status='unknown', Recl=RecordLength )
    end if
  end if

  This%NBuffer    =     isz
  allocate( This%Buffer(This%NBuffer) )

  if (i_Debug_Loc) then
    write(Logger%Unit,"(2x,'[OpenFile]: This%Name    = ',g0)") This%Name
    write(Logger%Unit,"(2x,'[OpenFile]: This%Unit    = ',g0)") This%Unit
    write(Logger%Unit,"(2x,'[OpenFile]: This%irecl   = ',g0)") This%irecl
    write(Logger%Unit,"(2x,'[OpenFile]: This%NBuffer = ',g0)") This%NBuffer
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine


! save unit 7 direct access records in memory to avoid problems with buffer flushing since record length is small compared to buffer size
! isz must be gt 2*irecl*NTraj+1 with irecl = 6 * natom - 6
! irec=0 is signal to write to disk.
Subroutine WriteVariable( This, irec, Var, i_Debug )

  use Error_Class           ,only:  Error
  use Parameters_Module     ,only:  Zero

  class(OutputFile_Type)                          ,intent(inout)  ::    This
  integer                                   ,intent(in)     ::    irec
  real(rkp) ,dimension(:)                   ,intent(in)     ::    Var         ! Variable to be written to memory or to a file
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc

  integer                                                   ::    NTraj
!   integer                                                   ::    Nneed
  integer                                                   ::    Nwrt
  integer                                                   ::    i, ig

! ! check for sufficient memory
!   if ( irec == 1 ) then
!     NTraj   =   int( Var(2) )           ! Extracting the number of trajectories
!     Nneed   =   NTraj * 2 * irecl + 1   ! Setting the required dimension
!     write(Logger%Unit,"('[wu7]: Required size: Nneed = ',g0)") Nneed
!     write(Logger%Unit,"('[wu7]: Current  size: isz   = ',g0)") isz
!     if ( Nneed > isz ) call Error( "Note enough space in isz" )
!   end if
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "WriteVariable")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
!   i_Debug_Loc = .True.
!   if (i_Debug_Loc) write(Logger%Unit,"(8x,'[WriteVariable]: Entering')")
! !   if (i_Debug_Loc) write(Logger%Unit,"(8x,'[WriteVariable]: Inputs:')")
!   if (i_Debug_Loc) write(Logger%Unit,"(8x,'[WriteVariable]:  -> irec = ',g0)") irec
! !   if (i_Debug_Loc) write(Logger%Unit,"(8x,'[WriteVariable]:  -> size(Var) = ',g0)") size(Var)

! Storing the input variable in memory (In This%Buffer)
  if ( irec > 0 ) then
!     if (i_Debug_Loc) write(Logger%Unit,"(8x,'[WriteVariable]: In memory')")
!     if (i_Debug_Loc) write(Logger%Unit,"(8x,'[WriteVariable]: -> This%irecl = ',g0)") This%irecl
    ig      =   (irec-1)*This%irecl
!     if (i_Debug_Loc) write(Logger%Unit,"(8x,'[WriteVariable]: -> ig = ',g0)") ig
    do i = 1,This%irecl
!     if (i_Debug_Loc) write(Logger%Unit,"(8x,'[WriteVariable]: -> i = ',g0,3x,'ig+i = ',g0)") i, ig+i
      if ( i <= size(Var) ) then
        This%Buffer(ig+i)   =   Var(i)
      else
        This%Buffer(ig+i)   =   Zero
      end if
    end do

! Write input variable to disk
  else
!     if (i_Debug_Loc) write(Logger%Unit,"(8x,'[WriteVariable]: In file')")
    NTraj   =   int( This%Buffer(2) )
    Nwrt    =   2 * NTraj + 1
!     if (i_Debug_Loc) write(Logger%Unit,"(8x,'[WriteVariable]: -> NTraj = ',g0)") NTraj
!     if (i_Debug_Loc) write(Logger%Unit,"(8x,'[WriteVariable]: -> Nwrt  = ',g0)") Nwrt
    do i = 1,Nwrt
      ig    =   (i-1)*This%irecl+1
!       if (i_Debug_Loc) write(Logger%Unit,"(8x,'[WriteVariable]: -> Nwrt  = ',g0)") Nwrt
      if (Formatted) then
        write(This%Unit,"(*(es15.8,1x),/)",rec=i) This%Buffer(ig:ig+This%irecl-1)
      else
        write(This%Unit,rec=i) This%Buffer(ig:ig+This%irecl-1)
      end if
    end do
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine


End Module
