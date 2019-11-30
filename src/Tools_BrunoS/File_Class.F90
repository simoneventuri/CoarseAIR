Module File_Class

!   use Error_Class               ,only:  Error

  implicit none

  private
  public    ::    File_Type

  Type                          ::      File_Type
    logical                     ::      Active = .False.
    integer                     ::      Unit
    character(:)  ,allocatable  ::      Name
    character(:)  ,allocatable  ::      Format
    integer                     ::      Status
    character(:)  ,allocatable  ::      ErrMsg
  contains
    procedure   ,public   ::    Open      =>    Open_File
    procedure   ,public   ::    Rewind    =>    Rewind_File
    procedure   ,public   ::    Close     =>    Close_File
    generic     ,public   ::    Write     =>    Write_C0D
    procedure   ,public   ::    Write_C0D
  End Type

  contains

Subroutine Open_File( This, FileName, iostat )

  class(File_Type)                                      ,intent(out)    ::      This                            !< Passed-object dummy argument
  character(*)                                          ,intent(in)     ::      FileName                        !< File name including relative path and extension
  integer   ,optional                                   ,intent(out)    ::      iostat

  allocate( This%Name , source = trim(adjustl(FileName)) )                                                      ! Setting the file name

!   inquire( File=This%Name, Exist=i_Exist )                                                                      ! Setting the file existence indicator
!   if ( .not.i_Exist ) call Error%Raise( "The file '"//This%Name//"' does not exists" )                                ! If the current FILE DOES NOT EXIST, then error

  open( NewUnit=This%Unit, File=This%Name, iostat=This%Status )                            ! Opening file for "read only mode"
!   if ( ios /= 0 ) call Error%Open( Unit=This%Unit, File=This%Name )                                             ! Checking for errors during file opening

  if ( present(iostat) ) iostat = This%Status

End Subroutine

Subroutine Close_File( This, Status, iostat )
  class(File_Type)                                      ,intent(inout)  ::      This                            !< Passed-object dummy argument
  character ,optional                                   ,intent(in)     ::      Status
  integer   ,optional                                   ,intent(out)    ::      iostat
  if ( present(Status) ) then
    close( This%Unit, iostat=This%Status, Status=Status )
  else
    close( This%Unit, iostat=This%Status )
  end if
  if ( present(iostat) ) iostat = This%Status
End Subroutine

Subroutine Rewind_File( This, iostat )
  class(File_Type)                                      ,intent(inout)  ::      This                            !< Passed-object dummy argument
  integer   ,optional                                   ,intent(out)    ::      iostat
  rewind( This%Unit, iostat=This%Status )
  if ( present(iostat) ) iostat = This%Status
End Subroutine

Subroutine Write_C0D( This, Variable, Format )
  class(File_Type)                                      ,intent(inout)  ::      This                            !< Passed-object dummy argument
  character(*)                                          ,intent(in)     ::      Variable
  character(*)                                ,optional ,intent(in)     ::      Format
  if ( present(Format) ) then
    write(This%Unit,Format,iostat=This%Status) Variable
  else
    write(This%Unit,*,     iostat=This%Status) Variable
  end if
End Subroutine

End Module
