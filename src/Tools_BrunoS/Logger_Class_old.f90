Module Logger_Class

  use ,intrinsic :: iso_fortran_env      ,only:  Output_Unit
  use LoggerItem_Class   ,only:  LoggerItem_Type

  implicit none

  private
  public        ::      Logger
  public        ::      Logger_Type
  public        ::      LogLevel_NOLOGS
  public        ::      LogLevel_ERROR
  public        ::      LogLevel_WARNING
  public        ::      LogLevel_INFO
  public        ::      LogLevel_DEBUG
  public        ::      LogLevel_HEAVYDEBUG
  public        ::      LogLevel_DEFAULT

  integer       ,parameter                      ::      LogLevel_NOLOGS      =   0
  integer       ,parameter                      ::      LogLevel_ERROR       =   1
  integer       ,parameter                      ::      LogLevel_WARNING     =   2
  integer       ,parameter                      ::      LogLevel_INFO        =   3
  integer       ,parameter                      ::      LogLevel_DEBUG       =   4
  integer       ,parameter                      ::      LogLevel_HEAVYDEBUG  =   5
  integer       ,parameter                      ::      LogLevel_DEFAULT     =   LogLevel_INFO

  integer       ,parameter                      ::      Indentation_Step          =     2
  character(*)  ,parameter                      ::      Default_Name_LogFile      =     "logfile.log"
  character(*)  ,parameter                      ::      Default_Logical_Format    =     "l1"
  character(*)  ,parameter                      ::      Default_Integer_Format    =     "i0"
  character(*)  ,parameter                      ::      Default_Real_Format       =     "g0"
  character(*)  ,parameter                      ::      Default_Character_Format  =     "a"
  character(*)  ,parameter                      ::      Default_Spacing2_Format   =     "3x"

  character(*)  ,parameter                      ::      Default_OpenStatus   = "REPLACE"
  character(*)  ,parameter                      ::      Default_OpenPosition = "REWIND"
  character(*)  ,parameter      ,dimension(2)   ::      Valid_OpenStatus     = ["REPLACE","OLD    "]
  character(*)  ,parameter      ,dimension(2)   ::      Valid_OpenPosition   = ["REWIND","APPEND"]

  logical ,parameter :: LocalDebug=.False.


  Type                                                  ::    Logger_Type
!     private
    character(:)        ,allocatable                    ::    FileName                                        !< Name of the log file
    integer                                             ::    Unit            =       Output_Unit             !< Unit of the log file
    logical                                             ::    Advancing       =       .True.
    integer                                             ::    iLogLev = 0
    integer                                             ::    NLevels = 0
    character(:)        ,allocatable ,dimension(:)      ::    ProcToLog                                         !< Names of the procedures for which the activation of logs should be enforced
    character(:)        ,allocatable ,dimension(:)      ::    ProcNotToLog                                      !< Names of the procedures for which the desactivation of logs should be enforced
    type(LoggerItem_Type)   ,dimension(:) ,allocatable  ::    Item                                            !< List of 'LoggerItem' objects
    type(LoggerItem_Type)   ,pointer                    ::    CurrentItem  => null()                                    !<
  contains
!   Procedures to initialize a Logger object and to perform file-related operations
    procedure   ,public                         ::      Initialize  =>  InitializeLogger                        !< Initializes a Logger object
    procedure   ,public                         ::      Reopen      =>  ReopenLogger                            !< Close and deleted current logfile, and then reopen it in 'append' mode
    procedure   ,public                         ::      Reset       =>  ResetLogger                             !< Reset the component of the Logger to their default values
    procedure   ,public                         ::      Close       =>  CloseLogger                             !< Close the logfile associated to a Logger object
    procedure   ,public                         ::      Backspace   =>  BackspaceLogger                         !< Backspaces the logfile associated to a Logger object
    procedure   ,public                         ::      Rewind      =>  RewindLogger                            !< Rewind the logfile associated to a Logger object
!     procedure   ,public                         ::      Deactivate



    procedure   ,public                         ::      Entering
    procedure   ,public                         ::      Exiting
    procedure   ,public                        ::      GetLogLevel
    procedure   ,public                        ::      StartDebugMode
    procedure   ,public                        ::      StopDebugMode



!   Procedures to enforce the activation of logs for a set of procedures
    generic     ,public                         ::      AddProcedureToLog => AddProcToLog0d, AddProcToLog1d
    procedure   ,public                         ::      GetProcedureToLog
    procedure   ,private                        ::      AddProcToLog0d
    procedure   ,private                        ::      AddProcToLog1d
!   Procedures to enforce the desactivation of logs for a set of procedures
    generic     ,public                         ::      AddProcedureNotToLog => AddProcNotToLog0d, AddProcNotToLog1d
    procedure   ,public                         ::      GetProcedureNotToLog
    procedure   ,private                        ::      AddProcNotToLog0d
    procedure   ,private                        ::      AddProcNotToLog1d
!   Procedures to access the properties of the Logger object
    generic     ,public                         ::      On => Active
!     procedure   ,public                         ::      Off
    procedure   ,public                         ::      Active                                                  !< Return a logical variable which indicates whether of not the Logger is active for the current log-level
    procedure   ,public                         ::      GetPrefix                                               !< Returns a character string corresponding to the prefix of the format used to write logs
    procedure   ,public                         ::      GetPath                                                 !< Get the path of the Logger in term of procedure names
    procedure   ,public                         ::      GetIndentation
!   Private procedures used during initialization or writing
    procedure   ,private                        ::      SetFileName                                             !< Sets the name of the file associated to the logger
    procedure   ,private                        ::      AddItem                                             !< Add an element to the list of 'LoggerItem' sub-objects and update the current index of log-level
    procedure   ,private                        ::      RemoveItem                                          !< Remove the last element from the list of 'LoggerItem' sub-objects and update the current index of log-level
    procedure   ,private                        ::      Set_Prefix_Indentation
    procedure   ,private                        ::      GetPrefixProcedure
!     procedure   ,private                        ::      Set_Local_Prefix
    procedure   ,private                        ::      Error_Open
    procedure   ,private                        ::      Write_NewLine
    procedure   ,private                        ::      Set_Advancing
    procedure   ,private                        ::      Set_Backspace
!   Procedures for writing logs
    procedure   ,private                        ::      Write_Blank_Line
    procedure   ,private                        ::      Write_1xV0
    procedure   ,private                        ::      Write_2xV0
    procedure   ,private                        ::      Write_3xV0
    procedure   ,private                        ::      Write_4xV0
    procedure   ,private                        ::      Write_5xV0
    procedure   ,private                        ::      Write_6xV0
    procedure   ,private                        ::      Write_7xV0
    procedure   ,private                        ::      Write_8xV0
    procedure   ,private                        ::      Write_9xV0
    procedure   ,private                        ::      Write_10xV0
    procedure   ,private                        ::      Write_11xV0
    procedure   ,private                        ::      Write_12xV0
    procedure   ,private                        ::      Write_13xV0
    procedure   ,private                        ::      Write_14xV0
    procedure   ,private                        ::      Write_15xV0
    procedure   ,private                        ::      Write_16xV0
    procedure   ,private                        ::      Write_17xV0
    procedure   ,private                        ::      Write_18xV0
    procedure   ,private                        ::      Write_20xV0
    procedure   ,private                        ::      Write_21xV0
    procedure   ,private                        ::      Write_22xV0
    procedure   ,private                        ::      Write_23xV0
    procedure   ,private                        ::      Write_24xV0
    procedure   ,private                        ::      Write_25xV0
    procedure   ,private                        ::      Write_26xV0
    procedure   ,private                        ::      Write_27xV0
    procedure   ,private                        ::      Write_28xV0
    procedure   ,private                        ::      Write_29xV0
    procedure   ,private                        ::      Write_30xV0
    procedure   ,private                        ::      Write_1xV0_1xV1
    procedure   ,private                        ::      Write_2xV0_1xV1
    procedure   ,private                        ::      Write_3xV0_1xV1
    procedure   ,private                        ::      Write_5xV0_1xV1
    procedure   ,private                        ::      Write_7xV0_1xV1
    procedure   ,private                        ::      Write_9xV0_1xV1
    procedure   ,private                        ::      Write_1xV0_2xV1
    procedure   ,private                        ::      Write_1xV0_2xV0V1
    procedure   ,private                        ::      Write_1xV0_3xV0V1
    procedure   ,private                        ::      Write_1xV0_1xV2
    procedure   ,private                        ::      Write_2xV0V1
    procedure   ,private                        ::      Write_3xV0V1
    procedure   ,private                        ::      Write_4xV0V1
    generic     ,public                         ::      Write   =>  Write_Blank_Line    , &
                                                                    Write_1xV0          , &
                                                                    Write_2xV0          , &
                                                                    Write_3xV0          , &
                                                                    Write_4xV0          , &
                                                                    Write_5xV0          , &
                                                                    Write_6xV0          , &
                                                                    Write_7xV0          , &
                                                                    Write_8xV0          , &
                                                                    Write_9xV0          , &
                                                                    Write_10xV0         , &
                                                                    Write_11xV0         , &
                                                                    Write_12xV0         , &
                                                                    Write_13xV0         , &
                                                                    Write_14xV0         , &
                                                                    Write_15xV0         , &
                                                                    Write_16xV0         , &
                                                                    Write_17xV0         , &
                                                                    Write_18xV0         , &
                                                                    Write_20xV0         , &
                                                                    Write_21xV0         , &
                                                                    Write_22xV0         , &
                                                                    Write_23xV0         , &
                                                                    Write_24xV0         , &
                                                                    Write_25xV0         , &
                                                                    Write_26xV0         , &
                                                                    Write_27xV0         , &
                                                                    Write_28xV0         , &
                                                                    Write_29xV0         , &
                                                                    Write_30xV0         , &
                                                                    Write_1xV0_1xV1     , &
                                                                    Write_2xV0_1xV1     , &
                                                                    Write_3xV0_1xV1     , &
                                                                    Write_5xV0_1xV1     , &
                                                                    Write_7xV0_1xV1     , &
                                                                    Write_9xV0_1xV1     , &
                                                                    Write_1xV0_2xV1     , &
                                                                    Write_1xV0_2xV0V1   , &
                                                                    Write_1xV0_3xV0V1   , &
                                                                    Write_1xV0_1xV2     , &
                                                                    Write_2xV0V1        , &
                                                                    Write_3xV0V1        , &
                                                                    Write_4xV0V1
  End Type

  type(Logger_Type) ,target                             ::      Logger

  Interface             Convert_Variable_To_String
    Module Procedure    Convert_Var0d_To_Str0d
    Module Procedure    Convert_Var1d_To_Str0d
    Module Procedure    Convert_Var1d_To_Str1d
    Module Procedure    Convert_Var2d_To_Str1d
  End Interface

  Interface             AddElementToArray
    Module Procedure    AddElementToArray_C0
    Module Procedure    AddElementToArray_C1
  End Interface

  Interface             Convert_To_String
    Module Procedure    Convert_Logical_To_String_0D , Convert_Logical_To_String_1D , Convert_Logical_To_Strings_1D
    Module Procedure    Convert_Int8_To_String_0D    , Convert_Int8_To_String_1D    , Convert_Int8_To_Strings_1D
    Module Procedure    Convert_Int16_To_String_0D   , Convert_Int16_To_String_1D    , Convert_Int16_To_Strings_1D
    Module Procedure    Convert_Int32_To_String_0D   , Convert_Int32_To_String_1D    , Convert_Int32_To_Strings_1D
    Module Procedure    Convert_Int64_To_String_0D   , Convert_Int64_To_String_1D    , Convert_Int64_To_Strings_1D
    Module Procedure    Convert_Real4_To_String_0D   , Convert_Real4_To_String_1D   , Convert_Real4_To_Strings_1D
    Module Procedure    Convert_Real8_To_String_0D   , Convert_Real8_To_String_1D   , Convert_Real8_To_Strings_1D
    Module Procedure    Convert_Real16_To_String_0D  , Convert_Real16_To_String_1D  , Convert_Real16_To_Strings_1D
    Module Procedure    Convert_String_To_String_0D  , Convert_String_To_String_1D  , Convert_String_To_Strings_1D
  End Interface

  Interface             Set_Optional_Argument
    Module Procedure    Set_Optional_Argument_Logical
    Module Procedure    Set_Optional_Argument_Integer
    Module Procedure    Set_Optional_Argument_Real
    Module Procedure    Set_Optional_Argument_Character
  End Interface

  contains

    ! **************************************************************************************************************
    !         PROCEDURES TO INITIALIZE A LOGGER OBJECT AND TO PERFORM FILE-RELATED OPERATIONS
    ! **************************************************************************************************************
!     Pure
    Subroutine InitializeLogger( This, FileName, Status, Position, Procedure, Indentation, i_Force_FileName )
      class(Logger_Type)                                    ,intent(inout)  ::      This                            !< Passed-object dummy argument corresponding
      character(*)                                          ,intent(in)     ::      FileName                        !< Name of the log file
      character(*)                                ,optional ,intent(in)     ::      Status
      character(*)                                ,optional ,intent(in)     ::      Position
      character(*)                                ,optional ,intent(in)     ::      Procedure                       !< Name of the calling procedure
      integer                                     ,optional ,intent(in)     ::      Indentation                     !< Indentation level
      logical                                     ,optional ,intent(in)     ::      i_Force_FileName
  character(:)  ,allocatable                                            ::      OpenStatus
  character(:)  ,allocatable                                            ::      OpenPosition
  integer                                                               ::      ios
  type(LoggerItem_Type)                                                 ::      LoggerItem
  call This%SetFileName( FileName, i_Force_FileName )                                                           ! Setting the name of the file associated to the Logger object
  call LoggerItem%Initialize( Procedure, Indentation )                                                          ! Initializing a 'LoggerItem' object with the input properties for the name of the procedure and the indentation level
  call This%AddItem( LoggerItem )                                                                               ! Adding the new 'LoggerItem' object to the list of 'LoggerItem' sub-objects
  OpenStatus   =       Get_OptOrDef_Value( Default_OpenStatus,   Valid_OpenStatus,   Status   )                 ! Getting the open status
  OpenPosition =       Get_OptOrDef_Value( Default_OpenPosition, Valid_OpenPosition, Position )                 ! Getting the open position
!   write(*,"('[InitializeLogger]: This%FileName = ',g0)") This%FileName
  open( NewUnit=This%Unit, File=This%FileName, Status=OpenStatus, Position=OpenPosition, IOStat=ios )           ! Opening the log file
  if ( ios /= 0 ) call This%Error_Open()                                                                        ! If opening error, then print error message and stop the code
End Subroutine
! Proce
    Pure Subroutine ResetLogger( This )
      class(Logger_Type)                                    ,intent(inout)    ::      This                            !< Passed-object dummy argument corresponding to the Logger object
  This%FileName       =   ''
  This%Unit           =   Output_Unit
  This%Advancing      =   .True.
  This%iLogLev        =   0
  This%NLevels        =   0
  if ( allocated(This%ProcToLog) )    deallocate( This%ProcToLog )
  if ( allocated(This%ProcNotToLog) ) deallocate( This%ProcNotToLog )
  if ( allocated(This%Item) )         deallocate( This%Item )
  if ( associated(This%CurrentItem) ) deallocate( This%CurrentItem )
    End Subroutine
!     Pure
    Subroutine ReopenLogger( This )
      class(Logger_Type)                                    ,intent(inout)  ::      This                            !< Passed-object dummy argument corresponding to the Logger object
  character(*)  ,parameter                                              ::      OpenStatus   = 'OLD'
  character(*)  ,parameter                                              ::      OpenPosition = 'APPEND'
  integer                                                               ::      ios
  integer :: i
  call This%Close( Status='DELETE')                                                                             ! Deleting the current Logger (the log file is deleted)
  open( NewUnit=This%Unit, File=This%FileName, Status=OpenStatus, Position=OpenPosition, IOStat=ios ) ! Opening the log file
  if ( ios /= 0 ) call This%Error_Open()
    End Subroutine
!     Pure
    Subroutine CloseLogger( This, Status )
      class(Logger_Type)                                    ,intent(inout)  ::      This                            !< Passed-object dummy argument corresponding to the Logger object
      character(*)                                ,optional ,intent(in)     ::      Status                          !< Status of the close instruction
  character(:)  ,allocatable                                            ::      Close_Status                    !< Local value of the close status
  Close_Status  =       ""
!   TODO: Add a checking that the "Status" string is a valid close status
  if ( present(Status) ) then
    Close_Status        =       Status
    close( This%Unit, Status=Close_Status )
  else
    close( This%Unit )
  end if
    End Subroutine
!     Pure
    Subroutine BackspaceLogger( This, Status )
      class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument corresponding to the Logger object
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  integer                                                               ::      ios
  backspace( This%Unit, iostat=ios )
  if ( present(Status) ) Status = ios
    End Subroutine
    Subroutine RewindLogger( This, Status )
      class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument corresponding to the Logger object
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  integer                                                               ::      ios
  rewind( This%Unit, iostat=ios )
  if ( present(Status) ) Status = ios
    End Subroutine


    Subroutine Entering( This, ProcedureName, LogLevel, DefLogLevel, MsgLogLevel, Writing, Active )
      class(Logger_Type)                                    ,intent(inout)  ::      This                          !< Passed-object dummy argument corresponding the Logger object
      character(*)                                          ,intent(in)     ::      ProcedureName                 !< Name of the calling procedure
      integer                                     ,optional ,intent(in)     ::      LogLevel
      integer                                     ,optional ,intent(in)     ::      DefLogLevel
      integer                                     ,optional ,intent(in)     ::      MsgLogLevel
      logical                                     ,optional ,intent(in)     ::      Writing                        !< Indicator whether the message should be written or not
      logical                                     ,optional ,intent(in)     ::      Active
  integer                                                               ::      Indentation
  type(LoggerItem_Type)                                                 ::      LoggerItem
  Indentation   =   This%GetIndentation()                                                                       ! Getting the current indentation levels
  call LoggerItem%Initialize( ProcedureName, Indentation, LogLevel, DefLogLevel, MsgLogLevel )                  ! Initializing a 'LoggerItem' object at current indentation level and with the input properties
  if ( LoggerItem%IsActive(ProcToLog=This%ProcToLog) ) call LoggerItem%Indent()                                 ! If the new 'LoggerItem' object is active, then increase its indentation level
  call This%AddItem( LoggerItem )                                                                               ! Adding the new 'LoggerItem' object to the list of 'LoggerItem' sub-objects

  if ( present(Active) ) then
    if (Active) then
      This%CurrentItem%MsgLogLevel = 0
    else
      This%CurrentItem%MsgLogLevel = 1000
    end if
  end if

  if ( Set_Optional_Argument(.True.,Writing) .and. This%Active() ) call This%Write("Entering")              ! Writing the "Entering" message


    End Subroutine
    Subroutine Exiting( This, Writing )
      class(Logger_Type)                                    ,intent(inout)  ::      This                          !< Passed-object dummy argument corresponding the Logger object
      logical                                     ,optional ,intent(in)     ::      Writing                       !< Indicator whether the message should be written or not                                                                                  ! Updating the prefix of the log message with the last indentation level and calling procedure
  if ( Set_Optional_Argument(.True.,Writing) .and. This%Active() ) call This%Write("Exiting")               ! Writing the "Exiting" message
  call This%RemoveItem()                                                                                  ! Removing the last procedure from the list of nested procedures
    End Subroutine

    ! **************************************************************************************************************
    !         PROCEDURES FOR ENFORCING THE ACTIVATION OF LOGS FOR A SET OF PROCEDURES
    ! **************************************************************************************************************
    Pure Subroutine GetProcedureToLog( This, ProcName )
      class(Logger_Type)                                    ,intent(in)     ::      This                        !< Passed-object dummy argument corresponding
      character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      ProcName                    !< Names of the procedures for which the activation of logs should be enforced
      integer ::  Length, i
!   allocate( ProcName, source = This%ProcToLog)
  Length  =   0
  do i = 1,size(This%ProcToLog)
    Length = max( Length , len_trim(This%ProcToLog(i)) )
  end do
  allocate( character(Length) :: ProcName(size(This%ProcToLog)) )
  ProcName = This%ProcToLog
    End Subroutine
    Pure Subroutine AddProcToLog0d( This, ProcName )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding
      character(*)                                          ,intent(in)     ::      ProcName                    !< Name of a procedure to add to the list of procedures for which the activation of logs should be enforced
  call AddElementToArray( ProcName, This%ProcToLog )
    End Subroutine
    Pure Subroutine AddProcToLog1d( This, ProcName )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding
      character(*)  ,dimension(:)                           ,intent(in)     ::      ProcName                    !< Name of the procedures to add to the list of procedures for which the activation of logs should be enforced
  call AddElementToArray( ProcName, This%ProcToLog )
    End Subroutine

    ! **************************************************************************************************************
    !         PROCEDURES FOR ENFORCING THE DESACTIVATION OF LOGS FOR A SET OF PROCEDURES
    ! **************************************************************************************************************
    Pure Subroutine GetProcedureNotToLog( This, ProcName )
      class(Logger_Type)                                    ,intent(in)     ::      This                        !< Passed-object dummy argument corresponding
      character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      ProcName                    !< Names of the procedures for which the desactivation of logs should be enforced
      integer ::  Length, i
!   allocate( ProcName, source = This%ProcNotToLog)
  Length  =   0
  do i = 1,size(This%ProcToLog)
    Length = max( Length , len_trim(This%ProcToLog(i)) )
  end do
  allocate( character(Length) :: ProcName(size(This%ProcNotToLog)) )
  ProcName = This%ProcNotToLog
    End Subroutine
    Pure Subroutine AddProcNotToLog0d( This, ProcName )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding
      character(*)                                          ,intent(in)     ::      ProcName                    !< Name of a procedure to add to the list of procedures for which the desactivation of logs should be enforced
  call AddElementToArray( ProcName, This%ProcNotToLog )
    End Subroutine
    Pure Subroutine AddProcNotToLog1d( This, ProcName )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding
      character(*)  ,dimension(:)                           ,intent(in)     ::      ProcName                    !< Name of the procedures to add to the list of procedures for which the desactivation of logs should be enforced
  call AddElementToArray( ProcName, This%ProcNotToLog )
    End Subroutine

    ! **************************************************************************************************************
    !                           PROCEDURES TO ACCESS THE PROPERTIES OF THE LOGGER OBJECT
    ! **************************************************************************************************************
!     Pure
    Function Active( This, LogLevel ) result(IsActive)
      class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument
      integer                                     ,optional ,intent(in)     ::      LogLevel
      logical                                                               ::      IsActive
  IsActive  =   .True.
  if ( associated(This%CurrentItem) ) IsActive = This%CurrentItem%IsActive( LogLevel, This%ProcToLog  )
    End Function

    Subroutine Deactivate( This )
      class(Logger_Type)                                    ,intent(inout)     ::      This                            !< Passed-object dummy argument
      if ( associated(This%CurrentItem) ) This%CurrentItem%MsgLogLevel = 1000
    End Subroutine

    Pure Function GetIndentation( This ) result(Indentation)
      class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument corresponding to the Logger object
      integer                                                               ::      Indentation
  Indentation   =   0
  if ( associated (This%CurrentItem) ) Indentation = This%CurrentItem%Indentation
    End Function
!     Pure
    Function GetPrefix( This, LogLevel, Error, Warning, Info, Debug, HeavyDebug ) result(Prefix)
      class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error
      logical                                     ,optional ,intent(in)     ::      Warning
      logical                                     ,optional ,intent(in)     ::      info
      logical                                     ,optional ,intent(in)     ::      Debug
      logical                                     ,optional ,intent(in)     ::      HeavyDebug
      character(:)  ,allocatable                                            ::      Prefix
!   Prefix    =   This%Set_Prefix_Indentation() // This%GetPrefixProcedure()
  Prefix    =   This%Set_Prefix_Indentation( LogLevel, Error, Warning, Info, Debug, HeavyDebug) // This%GetPrefixProcedure()
  if ( PresentAndTrue(Error) )       Prefix = Prefix // ",'<ERROR> '"
  if ( PresentAndTrue(Warning) )     Prefix = Prefix // ",'<WARNING> '"
  if ( PresentAndTrue(Info) )        Prefix = Prefix // ",'<INFO> '"
  if ( PresentAndTrue(Debug) )       Prefix = Prefix // ",'<DEBUG> '"
  if ( PresentAndTrue(HeavyDebug) )  Prefix = Prefix // ",'<HEAVYDEBUG> '"
  if ( len_trim(Prefix) > 1 ) then
    if ( This%Advancing ) Prefix = Prefix // ","  ! Only if in advancing mode, because when the cursor is not advancing, we do not want the prefix
  end if
  Prefix    =   "(" // Prefix
    End Function
    Pure Function GetPath( This ) result(Path)
      class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument
      character(:)  ,allocatable                                            ::      Path                            !< Path of the Logger in term of procedure names
  integer                                                               ::      i
  character(*)  ,parameter                                              ::      Separator = " > "
  Path        =       ""
  if ( .Not. allocated(This%Item) ) return
  do i = 1,size(This%Item)
  associate( LogLev => This%Item(i) )
    Path      =       Path // trim(LogLev%Name)
    if ( i /= size(This%Item) ) Path = Path // Separator
  end associate
  end do
    End Function
    Pure Function GetLogLevel( This ) result(LogLevel)
      class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument
      character(:)  ,allocatable                                            ::      LogLevel
  LogLevel = ""
  if ( associated (This%CurrentItem) ) LogLevel = LogLevelToString( This%CurrentItem%LogLevel )
    End Function

    Pure Subroutine StartDebugMode( This )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding
  if ( .Not. associated(This%CurrentItem) ) return
  This%CurrentItem%SavedMsgLogLevel  =   This%CurrentItem%MsgLogLevel
  This%CurrentItem%MsgLogLevel       =   LogLevel_DEBUG

!   if ( .Not. allocated(This%Item) ) return
!   associate( LogLev => This%Item(This%iLogLev) )
!     LogLev%SavedLogLevel  =   LogLev%LogLevel
!     LogLev%LogLevel       =   LogLevel_DEBUG
!   end associate
    End Subroutine
    Pure Subroutine StopDebugMode( This )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding
  if ( .Not. allocated(This%Item) ) return
  associate( LogLev => This%Item(This%iLogLev) )
    LogLev%MsgLogLevel       =   LogLev%SavedMsgLogLevel
  end associate
    End Subroutine



    ! **************************************************************************************************************
    !                           PRIVATE PROCEDURES USED DURING INITIALIZATION OR WRITING
    ! **************************************************************************************************************
!     Pure
    Subroutine SetFileName( This, FileName, i_Force_FileName )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding
      character(*)                                          ,intent(in)     ::      FileName                        !< FileName of the log file
      logical                                     ,optional ,intent(in)     ::      i_Force_FileName
#ifdef COARRAY
  character(10)                                                         ::      String                          ! Character string required to store the current image index (Required because write cannot make a length allocation)
#endif
  logical                                                               ::      i_Add_Image_Index
  i_Add_Image_Index     =       .True.
!   write(*,"('[SetFileName]: FileName = ',g0)") FileName
  This%FileName         =       trim(adjustl(FileName))                                                         ! Setting the FileName of the log file
!   write(*,"('[SetFileName]: This%FileName = ',g0)") This%FileName
  if ( present(i_Force_FileName) ) i_Add_Image_Index = .not. i_Force_FileName
#ifdef COARRAY
  if (i_Add_Image_Index) then
    if ( This_Image() == 1 ) then                                                                                 ! If the 1st image is considered
      This%FileName     =       trim(FileName)                                                                  ! Setting the FileName of the log file
    else                                                                                                          ! If an image other than the 1st one is considered, then adding the image index to the log file
      write(String,"(i0)") This_Image()                                                                           ! Converting the image index into a string
      This%FileName     =       trim(FileName) // "_" // trim(adjustl(String))                                  ! Setting the FileName of the log file with the image index
    end if                                                                                                        ! End if case on image index
  end if
#endif
    End Subroutine
    Pure Subroutine AddItem( This, LoggerItem )
      class(Logger_Type)                ,target                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding the Logger object
      type(LoggerItem_Type)                                 ,intent(in)     ::      LoggerItem                    !< LoggerItem object to be added to the list of LoggerItem sub-objects
  type(LoggerItem_Type)   ,dimension(:)   ,allocatable                  ::      LoggerItems
  if ( .not. allocated(This%Item) ) allocate( This%Item(0) )
  This%NLevels      =   size(This%Item)
  allocate( LoggerItems(This%NLevels+1) )
  LoggerItems(1:This%NLevels)    =       This%Item
  LoggerItems(This%NLevels+1)    =       LoggerItem
  call move_alloc( LoggerItems, This%Item )
  This%iLogLev      =   This%iLogLev + 1
  This%NLevels      =   size(This%Item)
  This%CurrentItem  =>  This%Item(This%iLogLev)
    End Subroutine
    Pure Subroutine RemoveItem( This )
      class(Logger_Type)                ,target                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding the Logger object
  type(LoggerItem_Type)  ,dimension(:)   ,allocatable                  ::      LoggerItems
  if ( .not. allocated(This%Item) ) allocate( This%Item(0) )
  This%NLevels      =   size(This%Item)
!   allocate( LoggerItems, source=This%Item(1:This%NLevels-1) )

  allocate( LoggerItems(1:This%NLevels-1) )
  LoggerItems  =  This%Item(1:This%NLevels-1)


  call move_alloc( LoggerItems, This%Item )
!   write(*,"('[RemoveItem]: Before: This%iLogLev = ',g0)") This%iLogLev
!   write(*,"('[RemoveItem]: Before: This%NLevels = ',g0)") This%NLevels
  This%iLogLev      =   This%iLogLev - 1
  This%NLevels      =   size(This%Item)
!   write(*,"('[RemoveItem]: After: This%iLogLev = ',g0)") This%iLogLev
!   write(*,"('[RemoveItem]: After: This%NLevels = ',g0)") This%NLevels
  if ( (This%iLogLev>0) .and. (This%iLogLev<=This%NLevels)  ) then
    This%CurrentItem  =>  This%Item(This%iLogLev)
  else
    This%CurrentItem  => null()
  end if
    End Subroutine
!     Pure
    Function Set_Prefix_Indentation( This, LogLevel, Error, Warning, Info, Debug, HeavyDebug ) result(String)
      class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error
      logical                                     ,optional ,intent(in)     ::      Warning
      logical                                     ,optional ,intent(in)     ::      info
      logical                                     ,optional ,intent(in)     ::      Debug
      logical                                     ,optional ,intent(in)     ::      HeavyDebug
      character(:)  ,allocatable                                            ::      String                          !< Character string corresponding to the indentation part of the format's prefix
  character(:)  ,allocatable                                            ::      Format, PrefixMode
  character(10)                                                         ::      Long_String
  integer                                                               ::      Indentation
  String       =       ""
  if ( This%iLogLev > 0 ) then
  associate( LogLev => This%Item(This%iLogLev) )
    if ( LogLev%Indentation > 0) then
!       write(*,"('[Set_Prefix_Indentation]:This%iLogLev = ',g0)") This%iLogLev
!       write(*,"('[Set_Prefix_Indentation]:associated(This%CurrentItem) = ',g0)") associated(This%CurrentItem)
!       write(*,"('[Set_Prefix_Indentation]:LogLev%Indentation = ',g0)") LogLev%Indentation
      Indentation   =   LogLev%Indentation
      Format        =   "(i0)"
      PrefixMode    =   ""
      if ( PresentAndTrue(Error) ) then
        Indentation =   Indentation - 1
        PrefixMode  =   "'E',"
      end if
      if ( PresentAndTrue(Warning) )then
        Indentation =   Indentation - 1
        PrefixMode  =   "'W',"
      end if
      if ( PresentAndTrue(Info) ) then
        Indentation =   Indentation - 1
        PrefixMode  =   "'I',"
      end if
      if ( PresentAndTrue(Debug) ) then
        Indentation =   Indentation - 1
        PrefixMode  =   "'D',"
      end if
      if ( PresentAndTrue(HeavyDebug) ) then
        Indentation =   Indentation - 1
        PrefixMode  =   "'H',"
      end if

!       write(*,"('[Set_Prefix_Indentation]:Indentation = ',g0)") Indentation
      write(Long_String,Format) Indentation
!       write(*,"('[Set_Prefix_Indentation]: Long_String = |',a,'|')") Long_String
!       write(Long_String,"(i0)") LogLev%Indentation
      String     =       PrefixMode // trim(adjustl(Long_String)) // "x,"
!       write(*,"('[Set_Prefix_Indentation]: String = |',a,'|')") String
    end if
  end associate
  end if
    End Function

    Function GetPrefixProcedure( This ) result(String)
      class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument
      character(:)  ,allocatable                                            ::      String                          !< Character string corresponding to the "procedure" part of the format's prefix
  String       =       ""
  if ( associated(This%CurrentItem) ) then
    if ( allocated(This%CurrentItem%Name) ) then
      String    =       "'[" // trim(This%CurrentItem%Name) // "]: '"
    end if
  end if
    End Function

!     Function Set_Local_Prefix( This, i_Prefix ) result(Local_Prefix)
!       class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument corresponding
!       logical                                     ,optional                 ::      i_Prefix                        !< Indicator whether the prefix has to be written
!       character(:)  ,allocatable                                            ::      Local_Prefix                    !< PRefix value
!     End Function


    Subroutine Set_Advancing( This, String, Advance )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      character(:)  ,allocatable                            ,intent(out)    ::      String
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
  String              =     "YES"
  This%Advancing      =     .True.
  if ( present(Advance) ) then
    if (.Not.Advance) then
      String          =     "NO"
      This%Advancing  =     .False.
    end if
  end if
    End Subroutine

    Subroutine Set_Backspace( This, Backspace )
      class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument corresponding to the Logger object
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
  if ( .Not. present(Backspace) ) return
  if ( .Not. Backspace ) return
  call This%Backspace()
    End Subroutine

    Subroutine Write_NewLine( This, NewLine )
      class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument corresponding to the Logger object
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
  if ( present(NewLine) ) then
    if (NewLine) write(This%Unit,*)
  end if
    End Subroutine

    Subroutine Error_Open( This )
      class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument corresponding
  write(*,"(4x,'[Initialize_Logger]: Error opening the Log file')")
  write(*,"(4x,'[Initialize_Logger]: FileName   = ',a)")  This%FileName
  write(*,"(4x,'[Initialize_Logger]: Unit   = ',i0)") This%Unit
  write(*,"(4x,'[Initialize_Logger]: Stopping the code')")
  error stop
    End Subroutine

    Subroutine Write_Blank_Line( This )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding
  write(This%Unit,*)
    End Subroutine

    Subroutine Write_1xV0(  This,                                                   &
                            V1,                                                     &
                            Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Underline, Fc, Fi, Fr, Fmt, &
                            F1                                                      )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error
      logical                                     ,optional ,intent(in)     ::      Warning
      logical                                     ,optional ,intent(in)     ::      info
      logical                                     ,optional ,intent(in)     ::      Debug
      logical                                     ,optional ,intent(in)     ::      HeavyDebug
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      Underline
      character(*)                                ,optional ,intent(in)     ::      F1
  character(:)  ,allocatable                                            ::      S1
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.

!   if ( .Not. This%Active(LogLevel) ) return

!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
!   LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_2xV0(  This,                                               &
                                  V1,  V2,                                            &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2                                             )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2
  character(:)  ,allocatable                                            ::      S1,  S2
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.

!   if ( .Not. This%Active(LogLevel) ) return

!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_3xV0(  This,                                               &
                                  V1,  V2,  V3,                                       &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3                                        )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3
  character(:)  ,allocatable                                            ::      S1,  S2,  S3
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_4xV0(   This,                                               &
                                  V1,  V2,  V3,  V4,                                  &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4                                   )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_5xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,                             &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5                              )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_6xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,                        &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6                         )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_7xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,                   &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7                    )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_8xV0(   This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,              &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8               )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_9xV0(   This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,         &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9          )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_10xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10    )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_11xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11,                                                &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11                                                 )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_12xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12,                                           &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12                                            )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_13xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13,                                      &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13                                       )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_14xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14,                                 &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14                                  )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_15xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15,                            &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15                             )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_16xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16,                       &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16                        )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_17xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17,                  &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17                   )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_18xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18,             &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18              )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_19xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19,        &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19         )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V19, S19, F19, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_20xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20    )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V19, S19, F19, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V20, S20, F20, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_21xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21,                                                &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21                                                 )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V19, S19, F19, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V20, S20, F20, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V21, S21, F21, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_22xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22,                                           &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22                                            )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V19, S19, F19, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V20, S20, F20, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V21, S21, F21, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V22, S22, F22, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_23xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23,                                      &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23                                       )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V19, S19, F19, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V20, S20, F20, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V21, S21, F21, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V22, S22, F22, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V23, S23, F23, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_24xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24,                                 &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24                                  )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V19, S19, F19, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V20, S20, F20, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V21, S21, F21, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V22, S22, F22, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V23, S23, F23, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V24, S24, F24, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_25xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24, V25,                            &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24, F25                             )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24, V25
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24, F25
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24, S25
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V19, S19, F19, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V20, S20, F20, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V21, S21, F21, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V22, S22, F22, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V23, S23, F23, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V24, S24, F24, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V25, S25, F25, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24, S25
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_26xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24, V25, V26,                       &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24, F25, F26                        )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24, V25, V26
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24, F25, F26
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24, S25, S26
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V19, S19, F19, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V20, S20, F20, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V21, S21, F21, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V22, S22, F22, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V23, S23, F23, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V24, S24, F24, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V25, S25, F25, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V26, S26, F26, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24, S25, S26
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_27xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24, V25, V26, V27,                  &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24, F25, F26, F27                   )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24, V25, V26, V27
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24, F25, F26, F27
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24, S25, S26, S27
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V19, S19, F19, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V20, S20, F20, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V21, S21, F21, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V22, S22, F22, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V23, S23, F23, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V24, S24, F24, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V25, S25, F25, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V26, S26, F26, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V27, S27, F27, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24, S25, S26, S27
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_28xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24, V25, V26, V27, V28,             &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24, F25, F26, F27, F28              )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24, V25, V26, V27, V28
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24, F25, F26, F27, F28
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24, S25, S26, S27, S28
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V19, S19, F19, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V20, S20, F20, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V21, S21, F21, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V22, S22, F22, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V23, S23, F23, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V24, S24, F24, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V25, S25, F25, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V26, S26, F26, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V27, S27, F27, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V28, S28, F28, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24, S25, S26, S27, S28
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_29xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24, V25, V26, V27, V28, V29,        &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24, F25, F26, F27, F28, F29         )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24, V25, V26, V27, V28, V29
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24, F25, F26, F27, F28, F29
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24, S25, S26, S27, S28, S29  ! Character strings corresponding to the input variable
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V19, S19, F19, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V20, S20, F20, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V21, S21, F21, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V22, S22, F22, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V23, S23, F23, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V24, S24, F24, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V25, S25, F25, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V26, S26, F26, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V27, S27, F27, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V28, S28, F28, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V29, S29, F29, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24, S25, S26, S27, S28, S29
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_30xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24, V25, V26, V27, V28, V29, V30,   &
                                  Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24, F25, F26, F27, F28, F29, F30    )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24, V25, V26, V27, V28, V29, V30
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24, F25, F26, F27, F28, F29, F30
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24, S25, S26, S27, S28, S29, S30  ! Character strings corresponding to the input variable
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V11, S11, F11, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V12, S12, F12, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V13, S13, F13, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V14, S14, F14, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V15, S15, F15, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V16, S16, F16, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V17, S17, F17, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V18, S18, F18, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V19, S19, F19, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V20, S20, F20, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V21, S21, F21, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V22, S22, F22, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V23, S23, F23, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V24, S24, F24, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V25, S25, F25, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V26, S26, F26, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V27, S27, F27, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V28, S28, F28, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V29, S29, F29, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V30, S30, F30, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24, S25, S26, S27, S28, S29, S30
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine



    Subroutine Write_1xV0_1xV1(  This,                                               &
                                    V1,  V2,                                            &
                                    Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F1,  F2                                             )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1
      class(*)      ,dimension(:)                           ,intent(in)     ::      V2
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2
  character(:)  ,allocatable                                            ::      S1,  S2
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
!   call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios


  select type (V2)
    type is (character(*))
      call Convert_To_String( V2, S2, VarFmt=F2, TypFormat=Fc, ComFmt=Fmt )
    class default;       ! Error
      call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  end select
  LineFormat    =       This%GetPrefix() // "*(a,a))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_2xV0_1xV1(  This,                                               &
                                    V1,  V2,  V3,                                       &
                                    Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F0, F1,  F2,  F3                                    )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2
      class(*)      ,dimension(:)                           ,intent(in)     ::      V3
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F0, F1,  F2,  F3
  character(:)  ,allocatable                                            ::      S0, S1,  S2, Fidx
  character(:)  ,allocatable  ,dimension(:)                             ::      S3
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Get_Integer_Format( size(V3), F0, Fidx, Status=ios );      if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
!                                             S1   i                                      S2  S3
  LineFormat    =       This%GetPrefix() // "a,"//Fidx//"," // Default_Spacing2_Format // ",a,a))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V3)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, i, S2, S3(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_3xV0_1xV1(  This,                                               &
                                    V1,  V2,  V3,  V4,                                  &
                                    Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F1,  F2,  F3,  F4                                   )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3
      class(*)      ,dimension(:)                           ,intent(in)     ::      V4
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_5xV0_1xV1(  This,                                               &
                                    V1,  V2,  V3,  V4,  V5,  V6,                        &
                                    Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F1,  F2,  F3,  F4,  F5,  F6                         )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5
      class(*)      ,dimension(:)                           ,intent(in)     ::      V6
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_7xV0_1xV1(  This,                                               &
                                    V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,              &
                                    Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8               )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7
      class(*)      ,dimension(:)                           ,intent(in)     ::      V8
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_9xV0_1xV1(  This,                                               &
                                    V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                    Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10    )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9
      class(*)      ,dimension(:)                           ,intent(in)     ::      V10
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V9 , S9 , F9 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V10, S10, F10, Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_1xV0_2xV1(  This,                                               &
                                      V1,  V2,  V3,                                       &
                                      Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                      F0, F1,  F2,  F3                                    )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1
      class(*)      ,dimension(:)                           ,intent(in)     ::      V2, V3
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F0,  F1,  F2,  F3
  character(:)  ,allocatable                                            ::      S0, S1, Fidx
  character(:)  ,allocatable  ,dimension(:)                             ::      S3,  S2
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Get_Integer_Format( size(V3), F0, Fidx, Status=ios );      if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
!                                             S1   i                                      S2  S3
  LineFormat    =       This%GetPrefix() // "a,"//Fidx//"," // Default_Spacing2_Format // ",a,a))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V3)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, i, S2(i), S3(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_1xV0_2xV0V1(  This,                                             &
                                      V1,  V2,  V3,  V4,  V5,                             &
                                      Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                      F0, F1,  F2,  F3,  F4,  F5                          )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1, V2, V4
      class(*)      ,dimension(:)                           ,intent(in)     ::      V3, V5
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F0, F1,  F2,  F3,  F4,  F5
  character(:)  ,allocatable                                            ::      S0, S1, S2, S4, Fidx
  character(:)  ,allocatable  ,dimension(:)                             ::      S3,  S5
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Get_Integer_Format( size(V3), F0, Fidx, Status=ios );      if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "a,"//Fidx//",*(" // Default_Spacing2_Format // ",a,a)))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V3)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, i, S2, S3(i), S4, S5(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_1xV0_3xV0V1(  This,                                             &
                                      V1,  V2,  V3,  V4,  V5,  V6,  V7,                   &
                                      Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                      F0, F1,  F2,  F3,  F4,  F5,  F6,  F7                )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1, V2, V4, V6
      class(*)      ,dimension(:)                           ,intent(in)     ::      V3, V5, V7
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F0, F1,  F2,  F3,  F4,  F5,  F6,  F7
  character(:)  ,allocatable                                            ::      S0, S1, S2, S4, S6, Fidx
  character(:)  ,allocatable  ,dimension(:)                             ::      S3,  S5,  S7
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Get_Integer_Format( size(V3), F0, Fidx, Status=ios );      if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "a,"//Fidx//",*(" // Default_Spacing2_Format // ",a,a)))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V3)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, i, S2, S3(i), S4, S5(i), S6, S7(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_1xV0_1xV2(  This,                                              &
                            V1, V2,                                                 &
                            Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                            F1, F2                                                  )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1
      class(*)      ,dimension(:,:)                         ,intent(in)     ::      V2
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1, F2
  character(:)  ,allocatable                                            ::      S1
  character(:)  ,allocatable  ,dimension(:)                             ::      S2
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  character(:)  ,allocatable                                            ::      Local_Prefix, FmtStr                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios, i, N
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status

  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1
!     N   =   This%GetPrefixLengthString()
  N   =   4 + 2
  call Convert_To_String(N,FmtStr)
  LineFormat  =  "("//FmtStr//"x,a)"
  do i = 1,size(S2,1)
    write( This%Unit,LineFormat, Advance=Adv, IOStat=ios ) S2(i)
  end do

  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_2xV0V1(  This,                                                   &
                              V1,  V2,  V3,  V4,                                      &
                              Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                              F1,  F2,  F3,  F4                                       )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1, V3
      class(*)      ,dimension(:)                           ,intent(in)     ::      V2, V4
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4
  character(:)  ,allocatable                                            ::      S1, S3
  character(:)  ,allocatable  ,dimension(:)                             ::      S2, S4
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V2)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, S2(i), S3, S4(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_3xV0V1(  This,                                                   &
                              V1,  V2,  V3,  V4,  V5,  V6,                            &
                              Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                              F1,  F2,  F3,  F4,  F5,  F6                             )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1, V3, V5
      class(*)      ,dimension(:)                           ,intent(in)     ::      V2, V4, V6
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6
  character(:)  ,allocatable                                            ::      S1, S3, S5
  character(:)  ,allocatable  ,dimension(:)                             ::      S2, S4, S6
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V2)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, S2(i), S3, S4(i), S5, S6(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine

    Subroutine Write_4xV0V1(  This,                                                   &
                              V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,                  &
                              Unused, i_Prefix, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                              F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8 )
      class(Logger_Type)                                    ,intent(inout)  ::      This                        !< Passed-object dummy argument corresponding to the Logger object
      class(*)                                              ,intent(in)     ::      V1, V3, V5, V7
      class(*)      ,dimension(:)                           ,intent(in)     ::      V2, V4, V6, V8
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                          !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                        !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                        !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      NewLine                         !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                         !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                       !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt                 !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8
  character(:)  ,allocatable                                            ::      S1, S3, S5, S7
  character(:)  ,allocatable  ,dimension(:)                             ::      S2, S4, S6, S8
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
!   character(:)  ,allocatable                                            ::      Local_Prefix                    ! Local prefix variable
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
!   Local_Prefix  =       This%Set_Local_Prefix( i_Prefix )                                                       ! Setting the local prefix value
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios

  LineFormat    =       This%GetPrefix() // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V2)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, S2(i), S3, S4(i), S5, S6(i), S7, S8(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
    End Subroutine




Function IsCurrentActive( This, ProcedureName, OptionalActive, DefaultActive ) result(WriteLogs_)
  class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument corresponding to the Logger object
  character(*)                                          ,intent(in)     ::      ProcedureName                       !< Name of the calling procedure
  logical                                     ,optional ,intent(in)     ::      OptionalActive
  logical                                     ,optional ,intent(in)     ::      DefaultActive
  logical                                                               ::      WriteLogs_

  logical                                                               ::      OptionalActive_
  logical                                                               ::      DefaultActive_
  integer                                                               ::      i

  DefaultActive_  = .False.
  if ( present(DefaultActive) ) DefaultActive_ = DefaultActive
  OptionalActive_ = DefaultActive_
  if ( present(OptionalActive) ) OptionalActive_ = OptionalActive

  WriteLogs_    =     OptionalActive_

  if ( allocated(This%ProcToLog) ) then
    do i = 1,size(This%ProcToLog)
      if ( .Not. ( LowerCase(ProcedureName) == LowerCase(This%ProcToLog(i)) ) ) cycle
      WriteLogs_  = .True.
      exit
    end do
  end if
!   if ( allocated(This%ProcToLog) ) then
!     if ( Procedure == This%ProcToLog(1) ) then ! @TODO
!       ActiveLoggerItem  = .True.
!     end if
!   end if
End Function


! **************************************************************************************************************
! **************************************************************************************************************
!                                           PRIVATE PROCEDURES
! **************************************************************************************************************
! **************************************************************************************************************


Subroutine Error_Unvalid_Value( Type_Value, Value, Valid_Values )
  character(*)                                          ,intent(in)     ::      Type_Value                      !< Type of value
  character(*)                                          ,intent(in)     ::      Value                           !< Erroneous value
  character(*)  ,dimension(:)                           ,intent(in)     ::      Valid_Values                    !< Valid values
  write(*,"(4x,'[Error_Unvalid_Value]: Error: Unvalid value ')")
  write(*,"(4x,'[Error_Unvalid_Value]: Type_Value   = ',a)") Type_Value
  write(*,"(4x,'[Error_Unvalid_Value]: Value        = ',a)") Value
  write(*,"(4x,'[Error_Unvalid_Value]: Valid_Values = ',*(a,3x))") Valid_Values(:)
  write(*,"(4x,'[Error_Unvalid_Value]: Stopping the code')")
  stop
End Subroutine


! This procedures set the format associated to a given variable which is returned in the character sting
! variable 'Format'. This variable is set to the following input variable using the following priority:
!  1) VarFmt: Used if present. It corresponds to a format specific to a given variable
!  2) TypFormat: Used if present and if 'VarFmt' is absent. It corresponds to a format specific to all variable of a given type.
!  2) ComFmt: Used if present and if 'VarFmt'/'TypeFormat' are absent. It corresponds to a format in common to all variables.
!  3) DefFormat: Used if all optional input variables are absent.
! TODO:
!   What should we do when a input format is not valid, that is, when the procedure 'Is_Valid_Format' returns a False value ?
!   Possible choices are: (1) an error, (2) something a bit less critical like taking the default format
!   For now, option (2) is conisdered

Pure Function Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt ) result(Format)

  character(*)                                          ,intent(in)     ::      DefFormat                     !< Default format:            Used if all optional format are absent
  character(*)                                ,optional ,intent(in)     ::      VarFmt                     !< Variable-specific format:  Used if present
  character(*)                                ,optional ,intent(in)     ::      TypFormat                     !< Type-specific format:      Used if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      ComFmt                     !< Common format:             Used if both the variable-specific and type-specific formats are absent
  character(:)  ,allocatable                                            ::      Format                        !< Selected format

  Format        =     DefFormat                                                                               ! Selecting the default format

  if ( present(VarFmt) ) then                                                                              ! If variable-specific format is present, then select it
    if ( Is_Valid_Format(VarFmt) ) Format = VarFmt                                                      ! Selecting the variable-specific format
    return                                                                                                    ! Exiting since the format has been selected
  end if                                                                                                      ! End if case on variable-specific format presence

  if ( present(TypFormat) ) then                                                                              ! If type-specific format is present, then select it
    if ( Is_Valid_Format(TypFormat) ) Format = TypFormat                                                      ! Selecting the type-specific format
    return                                                                                                    ! Exiting since the format has been selected
  end if                                                                                                      ! End if case on type-specific format presence

  if ( present(ComFmt) ) then                                                                              ! If common format is present, then select it
    if ( Is_Valid_Format(ComFmt) ) Format = ComFmt                                                      ! Selecting the common format
    return                                                                                                    ! Exiting since the format has been selected
  end if                                                                                                      ! End if case on common format presence

End Function

Pure Function Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt ) result(Format)

  character(*)                                          ,intent(in)     ::      DefFormat                     !< Default format:            Used if all optional format are absent
  character(*)                                ,optional ,intent(in)     ::      VarFmt                     !< Variable-specific format:  Used if present
  character(*)                                ,optional ,intent(in)     ::      TypFormat                     !< Type-specific format:      Used if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      ComFmt                     !< Common format:             Used if both the variable-specific and type-specific formats are absent
  character(:)  ,allocatable                                            ::      Format                        !< Selected format

  Format        =     "*(" // DefFormat // ",3x)"                                                             ! Selecting the default format

  if ( present(VarFmt) ) then                                                                              ! If variable-specific format is present, then select it
    if ( Is_Valid_Format(VarFmt) ) Format = "*(" // VarFmt // ",3x)"                                    ! Selecting the variable-specific format
    return                                                                                                    ! Exiting since the format has been selected
  end if                                                                                                      ! End if case on variable-specific format presence

  if ( present(TypFormat) ) then                                                                              ! If type-specific format is present, then select it
    if ( Is_Valid_Format(TypFormat) ) Format = "*(" // TypFormat // ",3x)"                                    ! Selecting the type-specific format
    return                                                                                                    ! Exiting since the format has been selected
  end if                                                                                                      ! End if case on type-specific format presence

  if ( present(ComFmt) ) then                                                                              ! If common format is present, then select it
    if ( Is_Valid_Format(ComFmt) ) Format = "*(" // ComFmt // ",3x)"                                    ! Selecting the common format
    return                                                                                                    ! Exiting since the format has been selected
  end if                                                                                                      ! End if case on common format presence

End Function


! TODO: Implement the algo telling if a given string is a valid fortran format.
Pure Elemental Function Is_Valid_Format( String ) result(Valid_Format)
  character(*)                                          ,intent(in)     ::      String
  logical                                                               ::      Valid_Format
  character(:)  ,allocatable                                            ::      String_Loc
  String_Loc = String
  Valid_Format  =       .True.
End Function


! **************************************************************************************************************
! **************************************************************************************************************
!                                           TOOLS
! **************************************************************************************************************
! **************************************************************************************************************

! This Function set and check a valid open status.
! If a valid optional open status is passed, then it is set other wise the default open status is taken
Function Get_OptOrDef_Value( Default_Value, Valid_Values, Optional_Value ) result( Output_Value )
  character(*)                                          ,intent(in)     ::      Default_Value                   !< Default value used if no optional value
  character(*)  ,dimension(:)                           ,intent(in)     ::      Valid_Values                    !< Valid values used to check validity of optional values if present
  character(*)                                ,optional ,intent(in)     ::      Optional_Value                  !< Optional values used if present and valid
  character(:)  ,allocatable                                            ::      Output_Value                    !< Output values
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
  character(*)                                  ,intent(in)     ::      Value                                   !< Value to be checked for validity
  character(*)          ,dimension( : )         ,intent(in)     ::      Valid_Values                            !< Valid values used for validity check
  logical                                                       ::      Valid                                   !< Indicator of input object validity
  integer                                                       ::      i                                       ! Index of valid strings
  Valid         =       .false.                                                                                 ! Initialization of the object validity indicator to false
  do i = 1,size(Valid_Values)                                                                                   ! Loop on all valid strings
    if ( trim(Value) == trim(Valid_Values(i)) ) Valid = .True.                                                  ! If the object if found in the list of valid strings, then setting validity indicator to True
  end do                                                                                                        ! End do loop on valid strings
End Function

! Subroutine Add_Element_To_Array( Element, Array )
!   implicit none
!   character(*)                                          ,intent(in)     ::      Element
!   character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      Array
!   integer                                                               ::      Length
!   character(:)  ,dimension(:)   ,allocatable                            ::      Array_tmp
! #ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
!   integer                                                               ::      i
! #endif
!   if ( .not. allocated(Array) ) allocate( character(0) :: Array(0) )
!   Length        =       max( len(Array), len(Element) )
!   allocate( character(Length) :: Array_tmp(size(Array)+1) )
! #ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
! !   ------------------------------
!   do i = 1,size(Array)
!     Array_tmp(i)    =       Array(i)
!   end do
! !   ------------------------------
! #else
!   Array_tmp(1:size(Array))    =       Array     ! COMPILER_BUG:GFORTRAN
! #endif
!   Array_tmp(size(Array)+1)    =       Element
!   call move_alloc( Array_tmp, Array )
! End Subroutine

Pure Subroutine AddElementToArray_C0( Element, Array )
  character(*)                                          ,intent(in)     ::      Element
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      Array
  character(:)  ,dimension(:)   ,allocatable                            ::      List_Elements
  integer                                                               ::      Length
  if ( .not. allocated(Array) ) allocate( character(0) :: Array(0) )
!   allocate( List_Elements, source = [Array,Element] )
  Length   = Max_Len_Trim( Array )
  allocate( character(Length) :: List_Elements(size(Array)+1) )
  List_Elements = [Array,Element]
  call move_alloc( List_Elements, Array )
End Subroutine

Pure Subroutine AddElementToArray_C1( Elements, Array )
  character(*)  ,dimension(:)                           ,intent(in)     ::      Elements
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      Array
  character(:)  ,dimension(:)   ,allocatable                            ::      List_Elements
  integer                                                               ::      Length
  if ( .not. allocated(Array) ) allocate( character(0) :: Array(0) )
!   allocate( List_Elements, source = [Array,Elements] )
  Length   = Max_Len_Trim( [Array,Elements] )
  allocate( character(Length) :: List_Elements(size(Array)+size(Elements)) )
  List_Elements = [Array,Elements]
  call move_alloc( List_Elements, Array )
End Subroutine

Pure Subroutine Remove_Element_From_Array( Array )
  implicit none
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      Array
  integer                                                               ::      NElements
  integer                                                               ::      Length
  character(:)  ,dimension(:)   ,allocatable                            ::      Array_tmp
#ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
  integer                                                       ::      n, i
#endif
  if ( .not. allocated(Array) ) allocate( character(0) :: Array(0) )
  NElements     =       size(Array)
  Length        =       len(Array(1:NElements-1))
#ifdef GFORTRAN_WORKAROUND_SOURCE_ALLOCATION
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
  character(*)  ,dimension(:)                           ,intent(in)     ::      Input_String
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

Pure Function Max_Len_Trim( Strings ) result( Length )
  character(*)  ,dimension(:)                   ,intent(in)             ::      Strings                         !< Array of character string
  integer                                                               ::      Length                          !< Maximum length without trailling blanks along all elements of the input string array
  integer                                                               ::      i                               ! Index of string' elements
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

! These procedures converts an input variable into a string using a given format.
! If all the optional input format are absent from the calling sequence, then the default format is used which
! depende on the type of variable. Otherwise, then format 'VarFmt' is used if present, or the format 'TypFormat'
! is used if present.

Subroutine Convert_Var0d_To_Str0d( Variable, String, VarFmt, ChaFmt, IntFmt, ReaFmt, ComFmt, Status )
  use iso_fortran_env ,only:  int8, int16, int32, int64
  class(*)                                              ,intent(in)     ::      Variable                        !< Variable to be converted into a character string
  character(:)  ,allocatable                            ,intent(out)    ::      String                          !< Output character string corresponding to the input variable
  character(*)                                ,optional ,intent(in)     ::      VarFmt                          !< Variable-specific format: Used if present
  character(*)                                ,optional ,intent(in)     ::      ChaFmt                          !< Character format: Used if variable is of type 'character(*)'and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      IntFmt                          !< Integer format: Used if variable is of type 'integer(*)'and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      ReaFmt                          !< Reals format: Used if variable is of type 'real(*)' and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      ComFmt                          !< Common format: Used if the variable-specific format 'Fv'  and the type-specific format 'Fc/Fi/Fr' associated to current variable's type are absent
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  integer                                                               ::      Local_Status
  select type (Variable)
    type is (logical);        call Convert_To_String( Variable, String, VarFmt=VarFmt,                   ComFmt=ComFmt, Status=Local_Status )
    type is (integer(int8));  call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(int16)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(int32)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(int64)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(4));        call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(8));        call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(16));       call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (character(*));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ChaFmt, ComFmt=ComFmt, Status=Local_Status )
    class default;       ! Error
      String        =       "?"
      Local_Status  =       -1
  end select
  if ( present(Status) ) Status = Local_Status
End Subroutine

Subroutine Convert_Var1d_To_Str0d( Variable, String, VarFmt, ChaFmt, IntFmt, ReaFmt, ComFmt, Status )
  class(*)      ,dimension(:)                           ,intent(in)     ::      Variable                        !< Variable to be converted into a character string
  character(:)  ,allocatable                            ,intent(out)    ::      String                          !< Output character string corresponding to the input variable
  character(*)                                ,optional ,intent(in)     ::      VarFmt                          !< Variable-specific format: Used if present
  character(*)                                ,optional ,intent(in)     ::      ChaFmt                          !< Character format: Used if variable is of type 'character(*)'and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      IntFmt                          !< Integer format: Used if variable is of type 'integer(*)'and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      ReaFmt                          !< Reals format: Used if variable is of type 'real(*)' and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      ComFmt                          !< Common format: Used if the variable-specific format 'Fv'  and the type-specific format 'Fc/Fi/Fr' associated to current variable's type are absent
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  integer                                                               ::      Local_Status
!
!   select type (Variable)
!     type is (character(:))
!       write(*,"('******************************* [Convert_Var1d_To_Str0d]: Variable = ',*(g0,3x))") Variable
!   end select


  select type (Variable)
    type is (logical);      call Convert_To_String( Variable, String, VarFmt=VarFmt,                   ComFmt=ComFmt, Status=Local_Status )
    type is (integer(1));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(2));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(4));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(8));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(4));      call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(8));      call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(16));     call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (character(*))
!       write(*,"('******************************* [Convert_Var1d_To_Str0d]: Variable = ',*(g0,3x))") Variable
      call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ChaFmt, ComFmt=ComFmt, Status=Local_Status )

    class default;       ! Error
      String        =       "?"
      Local_Status  =       -1
  end select
  if ( present(Status) ) Status = Local_Status
End Subroutine

Subroutine Convert_Var1d_To_Str1d( Variable, String, VarFmt, ChaFmt, IntFmt, ReaFmt, ComFmt, Status )
  class(*)      ,dimension(:)                           ,intent(in)     ::      Variable                        !< Variable to be converted into a character string
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String                          !< Output character string corresponding to the input variable
  character(*)                                ,optional ,intent(in)     ::      VarFmt                          !< Variable-specific format: Used if present
  character(*)                                ,optional ,intent(in)     ::      ChaFmt                          !< Character format: Used if variable is of type 'character(*)'and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      IntFmt                          !< Integer format: Used if variable is of type 'integer(*)'and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      ReaFmt                          !< Reals format: Used if variable is of type 'real(*)' and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      ComFmt                          !< Common format: Used if the variable-specific format 'Fv'  and the type-specific format 'Fc/Fi/Fr' associated to current variable's type are absent
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  integer                                                               ::      Local_Status
  select type (Variable)
    type is (logical);      call Convert_To_String( Variable, String, VarFmt=VarFmt,                   ComFmt=ComFmt, Status=Local_Status )
    type is (integer(1));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(2));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(4));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(8));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(4));      call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(8));      call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(16));     call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (character(*)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ChaFmt, ComFmt=ComFmt, Status=Local_Status )
    class default;       ! Error
      String        =       "?"
      Local_Status  =       -1
  end select
  if ( present(Status) ) Status = Local_Status
End Subroutine


Subroutine Convert_Var2d_To_Str1d( Variable, String, VarFmt, ChaFmt, IntFmt, ReaFmt, ComFmt, Status )
  class(*)      ,dimension(:,:)                         ,intent(in)     ::      Variable                        !< Variable to be converted into a character string
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String                          !< Output character string corresponding to the input variable
  character(*)                                ,optional ,intent(in)     ::      VarFmt                          !< Variable-specific format: Used if present
  character(*)                                ,optional ,intent(in)     ::      ChaFmt                          !< Character format: Used if variable is of type 'character(*)'and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      IntFmt                          !< Integer format: Used if variable is of type 'integer(*)'and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      ReaFmt                          !< Reals format: Used if variable is of type 'real(*)' and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      ComFmt                          !< Common format: Used if the variable-specific format 'Fv'  and the type-specific format 'Fc/Fi/Fr' associated to current variable's type are absent
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  integer                                                               ::      Local_Status
  character(:)  ,allocatable                                            ::      String_0d
  integer                                                               ::      i, Length
  allocate( character(10000) :: String(size(Variable,1)) )  ! Allocating to the number of rows
  select type (Variable)
    type is (logical)
      do i = 1,size(Variable,1)
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt,                   ComFmt=ComFmt, Status=Local_Status )
        String(i) = String_0d
      end do

    type is (integer(1))
      do i = 1,size(Variable,1)
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
        String(i) = String_0d
      end do

    type is (integer(2))
      do i = 1,size(Variable,1)
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
        String(i) = String_0d
      end do

    type is (integer(4))
      do i = 1,size(Variable,1)
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
        String(i) = String_0d
      end do

    type is (integer(8))
      do i = 1,size(Variable,1)
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
        String(i) = String_0d
      end do

    type is (real(4))
      do i = 1,size(Variable,1)
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
        String(i) = String_0d
      end do

    type is (real(8))
      do i = 1,size(Variable,1)
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
        String(i) = String_0d
      end do

    type is (real(16))
      do i = 1,size(Variable,1)
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
        String(i) = String_0d
      end do

    type is (character(*))
      do i = 1,size(Variable,1)
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=ChaFmt, ComFmt=ComFmt, Status=Local_Status )
        String(i) = String_0d
      end do

    class default;       ! Error
      deallocate( String )
!       allocate( String, source = ["?"] )

  Length = 1
  allocate( character(Length) :: String(1) )
  String = ["?"]


      Local_Status  =       -1
  end select
  String(:) = VecTrim_1D(String)
  if ( present(Status) ) Status = Local_Status
End Subroutine


Subroutine Convert_Logical_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  logical                                               ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Logical_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Int8_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  int8, int16, int32, int64
  integer(int8)                                         ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Int16_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  int8, int16, int32, int64
  integer(int16)                                        ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Int32_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  int8, int16, int32, int64
  integer(int32)                                        ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Int64_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  int8, int16, int32, int64
  integer(int64)                                        ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Real4_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  real(4)                                               ,intent(in)     ::      Variable                             !< Real number to be converted into a string
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Real_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Real8_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  real(8)                                               ,intent(in)     ::      Variable                             !< Real number to be converted into a string
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Real_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Real16_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  real(16)                                              ,intent(in)     ::      Variable                             !< Real number to be converted into a string
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Real_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

! Note: the 'present(ComFmt)' is removed from the if case because it will lead the a trim of the character
! when a common format is specified but is inteneted only for numeric variable (real, integer).

Subroutine Convert_String_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  character(*)                                          ,intent(in)     ::      Variable                             !< Real number to be converted into a string
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Character_Format
  if ( present(VarFmt) .or. present(TypFormat) .or. present(ComFmt) ) then
    Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
    Local_Format    =       "(" // Local_Format // ")"
    write( Long_String , Local_Format , iostat=ios ) Variable
    if ( ios /= 0 ) then
      String          =       Variable
    else
      String          =       trim(Long_String)
    end if
  else
    String          =       Variable
    ios             =       0
  end if
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Logical_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  logical       ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Logical_Format
  Local_Format    =       Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Int8_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  int8, int16, int32, int64
  integer(int8) ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
  Local_Format    =       Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Int16_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  int8, int16, int32, int64
  integer(int16) ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
  Local_Format    =       Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Int32_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  int8, int16, int32, int64
  integer(int32) ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
  Local_Format    =       Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Int64_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  int8, int16, int32, int64
  integer(int64),dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
  Local_Format    =       Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Real4_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  real(4)       ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Real_Format
  Local_Format    =       Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Real8_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  real(8)       ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Real_Format
  Local_Format    =       Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String, Local_Format, iostat=ios ) Variable
  if ( ios /= 0 ) then
    Local_Format    =   "(*(g0,1x))"
    write( Long_String, Local_Format, iostat=ios ) Variable
  end if
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Real16_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  real(16)      ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  character(*)  ,parameter                                              ::      DefFormat = Default_Real_Format
  Local_Format    =       Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_String_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  character(*)  ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,allocatable                            ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  integer                                                               ::      i
  character(*)  ,parameter                                              ::      DefFormat = Default_Character_Format
  if ( present(VarFmt) .or. present(TypFormat) .or. present(ComFmt) ) then
    Local_Format    =       Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt )
    Local_Format    =       "(" // Local_Format // ")"
    write( Long_String , Local_Format , iostat=ios ) Variable
    if ( ios == 0 ) then
      String          =       ""
      do i = 1,size(Variable)
        String        =       String // Variable(i) // "   "
      end do
      String          =       trim(String)
    else
      write(Long_String ,"(*(a,1x))") Variable
      String          =       trim(Long_String)
    end if
  else
    String          =       ""
!     write(*,"('******************************* [Convert_String_To_String_1D]: size(Variable) = ',g0)") size(Variable)
    do i = 1,size(Variable)
!       write(*,"('******************************* [Convert_String_To_String_1D]: i = ',g0)") i
!       write(*,"('******************************* [Convert_String_To_String_1D]: Variable(i) = ',g0)") Variable(i)
      String        =       String // Variable(i) // "   "
!       write(*,"('******************************* [Convert_String_To_String_1D]: String = ',g0)") String
!       write(*,"('******************************* [Convert_String_To_String_1D]: i = ',g0,3x,'i = ',g0,3x,'size(Variable) = ',g0,3x,'Variable(i) = ',g0,3x,'String = ',g0)") i, size(Variable), Variable(i), String
    end do
    String          =       trim(String)
    ios             =       0
  end if
  if ( present(Status) ) Status = ios
End Subroutine

Subroutine Convert_Logical_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  logical       ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, i, Length
  character(*)  ,parameter                                              ::      DefFormat = Default_Logical_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  Stat            =       0
  do i = 1,size(Variable)
    write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
    if ( ios /= 0 ) then
      write(Long_String(i) ,"(g0)") Variable(i)
      Stat        =       ios
    end if
  end do
!   allocate( String , source = VecTrim(Long_String), Stat=ios )
  Length = 0
  do i = 1,size(Long_String)
    Length = max( Length , len_trim(Long_String(i)) )
  end do
  allocate( character(Length) :: String(size(Long_String)) , Stat=ios )
  String(:) = VecTrim(Long_String)
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Subroutine Convert_Int8_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
 use iso_fortran_env ,only:  int8, int16, int32, int64
  integer(int8) ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, Length
!   integer(kind(Variable))                                               ::      i
  integer                                                               ::      i
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  Stat            =       0
  do i = 1,size(Variable)
    write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
    if ( ios /= 0 ) then
      write(Long_String(i) ,"(g0)") Variable(i)
      Stat        =       ios
    end if
  end do
!   allocate( String , source = VecTrim(Long_String), Stat=ios )
  Length = 0
  do i = 1,size(Long_String)
    Length = max( Length , len_trim(Long_String(i)) )
  end do
  allocate( character(Length) :: String(size(Long_String)) , Stat=ios )
  String(:) = VecTrim(Long_String)
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Subroutine Convert_Int16_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
 use iso_fortran_env ,only:  int8, int16, int32, int64
  integer(int16),dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, Length
!   integer(kind(Variable))                                               ::      i
  integer                                                               ::      i
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  Stat            =       0
  do i = 1,size(Variable)
    write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
    if ( ios /= 0 ) then
      write(Long_String(i) ,"(g0)") Variable(i)
      Stat        =       ios
    end if
  end do
!   allocate( String , source = VecTrim(Long_String), Stat=ios )
  Length = 0
  do i = 1,size(Long_String)
    Length = max( Length , len_trim(Long_String(i)) )
  end do
  allocate( character(Length) :: String(size(Long_String)) , Stat=ios )
  String(:) = VecTrim(Long_String)
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Subroutine Convert_Int32_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
 use iso_fortran_env ,only:  int8, int16, int32, int64
  integer(int32),dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, Length
!   integer(kind(Variable))                                               ::      i
  integer                                                               ::      i
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  Stat            =       0
  do i = 1,size(Variable)
    write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
    if ( ios /= 0 ) then
      write(Long_String(i) ,"(g0)") Variable(i)
      Stat        =       ios
    end if
  end do
!   allocate( String , source = VecTrim(Long_String), Stat=ios )
  Length = 0
  do i = 1,size(Long_String)
    Length = max( Length , len_trim(Long_String(i)) )
  end do
  allocate( character(Length) :: String(size(Long_String)) , Stat=ios )
  String(:) = VecTrim(Long_String)
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Subroutine Convert_Int64_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
 use iso_fortran_env ,only:  int8, int16, int32, int64
  integer(int64),dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, Length
!   integer(kind(Variable))                                               ::      i
  integer                                                               ::      i
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  Stat            =       0
  do i = 1,size(Variable)
    write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
    if ( ios /= 0 ) then
      write(Long_String(i) ,"(g0)") Variable(i)
      Stat        =       ios
    end if
  end do
!   allocate( String , source = VecTrim(Long_String), Stat=ios )
  Length = 0
  do i = 1,size(Long_String)
    Length = max( Length , len_trim(Long_String(i)) )
  end do
  allocate( character(Length) :: String(size(Long_String)) , Stat=ios )
  String(:) = VecTrim(Long_String)
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Subroutine Convert_Real4_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  real(4)       ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, i, Length
  character(*)  ,parameter                                              ::      DefFormat = Default_Real_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  Stat            =       0
  do i = 1,size(Variable)
    write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
    if ( ios /= 0 ) then
      write(Long_String(i) ,"(g0)") Variable(i)
      Stat        =       ios
    end if
  end do
!   allocate( String , source = VecTrim(Long_String), Stat=ios )
  Length = 0
  do i = 1,size(Long_String)
    Length = max( Length , len_trim(Long_String(i)) )
  end do
  allocate( character(Length) :: String(size(Long_String)) , Stat=ios )
  String(:) = VecTrim(Long_String)

  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Subroutine Convert_Real8_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  real(8)       ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, i, Length
  character(*)  ,parameter                                              ::      DefFormat = Default_Real_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  Stat            =       0
  do i = 1,size(Variable)
    write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
    if ( ios /= 0 ) then
      write(Long_String(i) ,"(g0)") Variable(i)
      Stat        =       ios
    end if
  end do
!   allocate( String , source = VecTrim(Long_String), Stat=ios )

  Length = 0
  do i = 1,size(Long_String)
    Length = max( Length , len_trim(Long_String(i)) )
  end do
  allocate( character(Length) :: String(size(Long_String)) , Stat=ios )
  String(:) = VecTrim(Long_String)

  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Subroutine Convert_Real16_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  real(16)      ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, i, Length
  character(*)  ,parameter                                              ::      DefFormat = Default_Real_Format
  Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  Stat            =       0
  do i = 1,size(Variable)
    write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
    if ( ios /= 0 ) then
      write(Long_String(i) ,"(g0)") Variable(i)
      Stat        =       ios
    end if
  end do
!   allocate( String , source = VecTrim(Long_String), Stat=ios )


  Length = 0
  do i = 1,size(Long_String)
    Length = max( Length , len_trim(Long_String(i)) )
  end do
  allocate( character(Length) :: String(size(Long_String)) , Stat=ios )
  String(:) = VecTrim(Long_String)

  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Subroutine Convert_String_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  character(*)  ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, i, Length
  character(*)  ,parameter                                              ::      DefFormat = Default_Character_Format
  if ( present(VarFmt) .or. present(TypFormat) .or. present(ComFmt) ) then
    Local_Format    =       Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt )
    Local_Format    =       "(" // Local_Format // ")"
    Stat            =       0
    do i = 1,size(Variable)
      write( Long_String(i) , Local_Format , iostat=ios ) Variable(i)
      if ( ios /= 0 ) then
        write(Long_String(i) ,"(g0)") Variable(i)
        Stat        =       ios
      end if
    end do
!     allocate( String , source = VecTrim(Long_String), Stat=ios )
    Length = 0
    do i = 1,size(Long_String)
      Length = max( Length , len_trim(Long_String(i)) )
    end do
    allocate( character(Length) :: String(size(Long_String)) , Stat=ios )
    String(:) = VecTrim(Long_String)
    if ( ios /= 0 ) Stat = ios
  else
!     allocate( String , source = VecTrim(Variable), Stat=Stat )
    Length = 0
    do i = 1,size(Variable)
      Length = max( Length , len_trim(Variable(i)) )
    end do
    allocate( character(Length) :: String(size(Variable)) , Stat=Stat )
    String(:) = VecTrim(Variable)
  end if
  if ( present(Status) ) Status = Stat
End Subroutine


! **************************************************************************************************************
! **************************************************************************************************************
!                                   PROCEDURES RELATED TO FORMATS
! **************************************************************************************************************
! **************************************************************************************************************

Subroutine Get_Integer_Format( Variable, Input_Format, Output_Format, Status )
  integer                                               ,intent(in)     ::      Variable
  character(*)                                ,optional ,intent(in)     ::      Input_Format
  character(:)  ,allocatable                            ,intent(out)    ::      Output_Format
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(:)  ,allocatable                                            ::      String
  integer                                                               ::      NDigits
  integer                                                               ::      ios
  if ( present(Input_Format) ) then
    Output_Format   =     Input_Format
    ios             =     0
  else
    NDigits         =       floor( log10( real( abs(Variable) ) ) ) + 1                                                     ! Getting the number of digits of the input integer number
    call Convert_To_String( NDigits, String, Status=ios )
    Output_Format   =       "(i" // trim(String) // ")"
  end if
  if ( present(Status) ) Status = ios
End Subroutine

Pure Function LowerCase(StrInp) result(StrOut)
  character(*)                  ,intent(in)     ::      StrInp                                                  !<
  character(:)  ,allocatable                    ::      StrOut
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
  StrOut    =   trim(adjustl(StrOut))
End Function
!
! Function Set_Scalar_Format_0D( Variable, Optional_Format ) result(Format)
!   class(*)                                              ,intent(in)     ::      Variable
!   character(*)                                ,optional ,intent(in)     ::      Optional_Format
!   character(:)  ,allocatable                                            ::      Format
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
!   class(*)      ,dimension(:)                           ,intent(in)     ::      Variable
!   character(*)                                ,optional ,intent(in)     ::      Optional_Format
!   character(:)  ,allocatable                                            ::      Format
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
!   class(*)                                              ,intent(in)     ::      Variable
!   character(*)                                ,optional ,intent(in)     ::      Optional_Format
!   character(:)  ,allocatable                                            ::      Format
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
!   class(*)      ,dimension(:)                           ,intent(in)     ::      Variable
!   character(*)                                ,optional ,intent(in)     ::      Optional_Format
!   character(:)  ,allocatable                                            ::      Format
!   select type (Variable)
!     type is (logical);          Format = Set_Generic_Vector_Format( Default_Logical_Format,    Optional_Format )
!     type is (integer);          Format = Set_Generic_Vector_Format( Default_Integer_Format,    Optional_Format )
!     type is (real(4));          Format = Set_Generic_Vector_Format( Default_Real_Format,       Optional_Format )
!     type is (real(8));          Format = Set_Generic_Vector_Format( Default_Real_Format,       Optional_Format )
!     type is (character(*));     Format = Set_Generic_Vector_Format( Default_Character_Format,  Optional_Format )
!     class default       ! Error
!   end select
! End Function

Pure Function Set_Optional_Argument_Logical( VarDef, VarOpt ) result(VarLoc)
  logical                                                       ,intent(in)     ::      VarDef                  !< Default value of the local variable if the optional variable is absent
  logical       ,optional                                       ,intent(in)     ::      VarOpt                  !< Optional varibale
  logical                                                                       ::      VarLoc                  !< Local variable to be set
  if ( present(VarOpt) ) then;  VarLoc = VarOpt                                                                 ! Setting the local variable to the optional variable if present ...
  else;                         VarLoc = VarDef; end if                                                         ! ... otherwise setting the local variable to the default variable
End Function

Pure Function Set_Optional_Argument_Integer( VarDef, VarOpt ) result(VarLoc)
  integer                                                       ,intent(in)     ::      VarDef                  !< Default value of the local variable if the optional variable is absent
  integer       ,optional                                       ,intent(in)     ::      VarOpt                  !< Optional varibale
  integer                                                                       ::      VarLoc                  !< Local variable to be set
  if ( present(VarOpt) ) then;  VarLoc = VarOpt                                                                 ! Setting the local variable to the optional variable if present ...
  else;                         VarLoc = VarDef; end if                                                         ! ... otherwise setting the local variable to the default variable
End Function

Pure Function Set_Optional_Argument_Real( VarDef, VarOpt ) result(VarLoc)
  real(8)                                                       ,intent(in)     ::      VarDef                  !< Default value of the local variable if the optional variable is absent
  real(8)       ,optional                                       ,intent(in)     ::      VarOpt                  !< Optional varibale
  real(8)                                                                       ::      VarLoc                  !< Local variable to be set
  if ( present(VarOpt) ) then;  VarLoc = VarOpt                                                                 ! Setting the local variable to the optional variable if present ...
  else;                         VarLoc = VarDef; end if                                                         ! ... otherwise setting the local variable to the default variable
End Function

Pure Function Set_Optional_Argument_Character( VarDef, VarOpt ) result(VarLoc)
  character(*)                                                  ,intent(in)     ::      VarDef                  !< Default value of the local variable if the optional variable is absent
  character(*)  ,optional                                       ,intent(in)     ::      VarOpt                  !< Optional varibale
  character(:)  ,allocatable                                                    ::      VarLoc                  !< Local variable to be set
  if ( present(VarOpt) ) then;  VarLoc = VarOpt                                                                 ! Setting the local variable to the optional variable if present ...
  else;                         VarLoc = VarDef; end if                                                         ! ... otherwise setting the local variable to the default variable
End Function


Pure Function LogLevelToString( LogLevel ) result(String)
  integer                                                       ,intent(in)     ::      LogLevel
  character(:)  ,allocatable                                                    ::      String
  select case (LogLevel)
    case( LogLevel_NOLOGS     ); String = "NOLOGS"
    case( LogLevel_ERROR      ); String = "ERROR"
    case( LogLevel_WARNING    ); String = "WARNING"
    case( LogLevel_INFO       ); String = "INFO"
    case( LogLevel_DEBUG      ); String = "DEBUG"
    case( LogLevel_HEAVYDEBUG:); String = "HEAVYDEBUG"
  end select
End Function

Pure Function IsPresent( Element, Array ) result(PresenceIndicator)
  character(*)                                          ,intent(in)     ::      Element
  character(:)  ,allocatable ,dimension(:)              ,intent(in)     ::      Array
  logical                                                               ::      PresenceIndicator
  integer                                                               ::      i
  PresenceIndicator   =   .False.
  if ( allocated(Array) ) then
    do i = 1,size(Array)
      if ( .Not. ( LowerCase(Element) == LowerCase(Array(i)) ) ) cycle
      PresenceIndicator  = .True.
      return
    end do
  end if
End Function

Pure Function PresentAndTrue( OptionalArgument ) result(Indicator)
  logical                                     ,optional ,intent(in)     ::      OptionalArgument
  logical                                                               ::      Indicator
  Indicator   =   .False.
  if ( Present(OptionalArgument) ) Indicator  = OptionalArgument
End Function

Pure Function VecTrim_1D( Input_String ) result(Output_String)
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

End Module
