Module Logger_Class

! #define Purity Pure
#define Purity

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


  Type                                                  ::    Logger_Type
    character(:)        ,allocatable                    ::    FileName                                        !< Name of the log file
    integer                                             ::    Unit            =       Output_Unit             !< Unit of the log file
    logical                                             ::    Advancing       =       .True.
    logical                                             ::    Activated       =       .True.
    integer                                             ::    iLogLev = 0
    integer                                             ::    NLevels = 0
    character(:)        ,allocatable ,dimension(:)      ::    ProcToLog                                         !< Names of the procedures for which the activation of logs should be enforced
    character(:)        ,allocatable ,dimension(:)      ::    ProcNotToLog                                      !< Names of the procedures for which the desactivation of logs should be enforced
    type(LoggerItem_Type)   ,dimension(:) ,allocatable  ::    Item                                            !< List of 'LoggerItem' objects
    type(LoggerItem_Type)   ,pointer                    ::    CurrentItem  => null()                                    !<
  contains
    Final                 ::    FinalizeLogger
!   Procedures to initialize a Logger object and to perform file-related operations
    procedure   ,public   ::    Initialize  =>  InitializeLogger                        !< Initializes a Logger object
    procedure   ,public   ::    Activate                                                !< Procedure to force activation of the Logger object
    procedure   ,public   ::    Deactivate                                              !< Procedure to force deactivation of the Logger object
    procedure   ,public   ::    Free        =>  FreeLogger                              !< Free a Logger object
    procedure   ,public   ::    Reopen      =>  ReopenLogger                            !< Close and deleted current logfile, and then reopen it in 'append' mode
    procedure   ,public   ::    Close       =>  CloseLogger                             !< Close the logfile associated to a Logger object
    procedure   ,public   ::    Backspace   =>  BackspaceLogger                         !< Backspaces the logfile associated to a Logger object
    procedure   ,public   ::    Rewind      =>  RewindLogger                            !< Rewind the logfile associated to a Logger object
    procedure   ,public   ::    Flush       =>  FlushLogger                             !< Flush the unit associated to a Logger object
    procedure   ,public   ::    Entering
    procedure   ,public   ::    Exiting
    procedure   ,public   ::    GetLogLevel
    procedure   ,public   ::    StartDebugMode
    procedure   ,public   ::    StopDebugMode
!   Procedures to enforce the activation of logs for a set of procedures
    generic     ,public   ::    AddProcedureToLog  => AddProcToLog0d, AddProcToLog1d
    procedure   ,public   ::    GetProcedureToLog
    procedure   ,private  ::    AddProcToLog0d
    procedure   ,private  ::    AddProcToLog1d
!   Procedures to enforce the desactivation of logs for a set of procedures
    generic     ,public   ::    AddProcedureNotToLog => AddProcNotToLog0d, AddProcNotToLog1d
    procedure   ,public   ::    GetProcedureNotToLog
    procedure   ,private  ::    AddProcNotToLog0d
    procedure   ,private  ::    AddProcNotToLog1d
!   Procedures to access the properties of the Logger object
    generic     ,public   ::    On => Active
    procedure   ,public   ::    Active                                                  !< Return a logical variable which indicates whether of not the Logger is active for the current log-level
    procedure   ,public   ::    GetPrefix                                               !< Returns a character string corresponding to the prefix of the format used to write logs
    procedure   ,public   ::    GetPath                                                 !< Get the path of the Logger in term of procedure names
    procedure   ,public   ::    GetIndentation
!   Private procedures used during initialization or writing
    procedure   ,private  ::    SetFileName                                             !< Sets the name of the file associated to the logger
    procedure   ,private  ::    AddItem                                             !< Add an element to the list of 'LoggerItem' sub-objects and update the current index of log-level
    procedure   ,private  ::    RemoveItem                                          !< Remove the last element from the list of 'LoggerItem' sub-objects and update the current index of log-level
    procedure   ,private  ::    Set_Prefix_Indentation
    procedure   ,private  ::    GetPrefixProcedure
!     procedure   ,private  ::    Set_Local_Prefix
    procedure   ,private  ::    Error_Open
    procedure   ,private  ::    Write_NewLine
    procedure   ,private  ::    Set_Advancing
    procedure   ,private  ::    Set_Backspace
!   Procedures for writing logs
    procedure   ,public   ::    Calling => CallingProcedure
    procedure   ,private  ::    Write_Blank_Line
    procedure   ,private  ::    Write_1xV0
    procedure   ,private  ::    Write_2xV0
    procedure   ,private  ::    Write_3xV0
    procedure   ,private  ::    Write_4xV0
    procedure   ,private  ::    Write_5xV0
    procedure   ,private  ::    Write_6xV0
    procedure   ,private  ::    Write_7xV0
    procedure   ,private  ::    Write_8xV0
    procedure   ,private  ::    Write_9xV0
    procedure   ,private  ::    Write_10xV0
    procedure   ,private  ::    Write_11xV0
    procedure   ,private  ::    Write_12xV0
    procedure   ,private  ::    Write_13xV0
    procedure   ,private  ::    Write_14xV0
    procedure   ,private  ::    Write_15xV0
    procedure   ,private  ::    Write_16xV0
    procedure   ,private  ::    Write_17xV0
    procedure   ,private  ::    Write_18xV0
    procedure   ,private  ::    Write_20xV0
    procedure   ,private  ::    Write_21xV0
    procedure   ,private  ::    Write_22xV0
    procedure   ,private  ::    Write_23xV0
    procedure   ,private  ::    Write_24xV0
    procedure   ,private  ::    Write_25xV0
    procedure   ,private  ::    Write_26xV0
    procedure   ,private  ::    Write_27xV0
    procedure   ,private  ::    Write_28xV0
    procedure   ,private  ::    Write_29xV0
    procedure   ,private  ::    Write_30xV0
    procedure   ,private  ::    Write_1xV1
    procedure   ,private  ::    Write_1xV0_1xV1
    procedure   ,private  ::    Write_2xV0_1xV1
    procedure   ,private  ::    Write_3xV0_1xV1
    procedure   ,private  ::    Write_5xV0_1xV1
    procedure   ,private  ::    Write_7xV0_1xV1
    procedure   ,private  ::    Write_9xV0_1xV1
    procedure   ,private  ::    Write_1xV0_2xV1
    procedure   ,private  ::    Write_1xV0_2xV0V1
    procedure   ,private  ::    Write_1xV0_3xV0V1
    procedure   ,private  ::    Write_1xV0_1xV2
    procedure   ,private  ::    Write_2xV0V1
    procedure   ,private  ::    Write_3xV0V1
    procedure   ,private  ::    Write_4xV0V1
    generic     ,public   ::    Write   =>  Write_Blank_Line    , &
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
                                            Write_1xV1          , &
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

  Interface

    ! **************************************************************************************************************
    !         PROCEDURES TO INITIALIZE A LOGGER OBJECT AND TO PERFORM FILE-RELATED OPERATIONS
    ! **************************************************************************************************************
    Module Subroutine InitializeLogger( This, FileName, Status, Position, Procedure, Indentation, i_Force_FileName )
      class(Logger_Type)                                    ,intent(inout)  ::  This                            !< Passed-object dummy argument
      character(*)                                          ,intent(in)     ::      FileName                    !< Name of the log file
      character(*)                                ,optional ,intent(in)     ::      Status
      character(*)                                ,optional ,intent(in)     ::      Position
      character(*)                                ,optional ,intent(in)     ::      Procedure                   !< Names of the calling procedures to reach the current procedure. Each procedure neames are separated by the '>' character.  Only used by the Logger object
      integer                                     ,optional ,intent(in)     ::      Indentation                 !< Indentation level
      logical                                     ,optional ,intent(in)     ::      i_Force_FileName
    End Subroutine
    Purity Module Subroutine Activate( This )
      class(Logger_Type)                                    ,intent(inout)  ::  This                            !< Passed-object dummy argument corresponding the Logger object
    End Subroutine
    Purity Module Subroutine Deactivate( This )
      class(Logger_Type)                                    ,intent(inout)  ::  This                            !< Passed-object dummy argument corresponding the Logger object
    End Subroutine
    Module Subroutine FreeLogger( This )
      class(Logger_Type)                                    ,intent(inout)  ::  This                            !< Passed-object dummy argument
    End Subroutine
    Module Subroutine FinalizeLogger( This )
      type(Logger_Type)                                     ,intent(inout)  ::  This                            !< Passed-object dummy argument
    End Subroutine
    Module Subroutine ReopenLogger( This )
      class(Logger_Type)                                    ,intent(inout)  ::  This                            !< Passed-object dummy argument
    End Subroutine
    Module Subroutine CloseLogger( This, Status )
      class(Logger_Type)                                    ,intent(inout)  ::  This                            !< Passed-object dummy argument
      character(*)                                ,optional ,intent(in)     ::      Status                      !< Status of the close instruction
    End Subroutine
    Module Subroutine BackspaceLogger( This, Status )
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
    End Subroutine
    Module Subroutine RewindLogger( This, Status )
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
    End Subroutine
    Module Subroutine FlushLogger( This )
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument
    End Subroutine


    Module Subroutine Entering( This, ProcedureName, LogLevel, DefLogLevel, MsgLogLevel, Writing )
      class(Logger_Type)                                    ,intent(inout)  ::  This                          !< Passed-object dummy argument corresponding the Logger object
      character(*)                                          ,intent(in)     ::      ProcedureName                 !< Name of the calling procedure
      integer                                     ,optional ,intent(in)     ::      LogLevel
      integer                                     ,optional ,intent(in)     ::      DefLogLevel
      integer                                     ,optional ,intent(in)     ::      MsgLogLevel
      logical                                     ,optional ,intent(in)     ::      Writing                        !< Indicator whether the message should be written or not
    End Subroutine
    Module Subroutine Exiting( This, Writing )
      class(Logger_Type)                                    ,intent(inout)  ::  This                          !< Passed-object dummy argument corresponding the Logger object
      logical                                     ,optional ,intent(in)     ::      Writing                       !< Indicator whether the message should be written or not                                                                                  ! Updating the prefix of the log message with the last indentation level and calling procedure
    End Subroutine

    ! **************************************************************************************************************
    !         PROCEDURES FOR ENFORCING THE ACTIVATION OF LOGS FOR A SET OF PROCEDURES
    ! **************************************************************************************************************
    Purity Module Subroutine GetProcedureToLog( This, ProcName )
      class(Logger_Type)                                    ,intent(in)     ::  This                        !< Passed-object dummy argument corresponding
      character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      ProcName                    !< Names of the procedures for which the activation of logs should be enforced
    End Subroutine
    Purity Module Subroutine AddProcToLog0d( This, ProcName )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument corresponding
      character(*)                                          ,intent(in)     ::      ProcName                    !< Name of a procedure to add to the list of procedures for which the activation of logs should be enforced
    End Subroutine
    Purity Module Subroutine AddProcToLog1d( This, ProcName )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument corresponding
      character(*)  ,dimension(:)                           ,intent(in)     ::      ProcName                    !< Name of the procedures to add to the list of procedures for which the activation of logs should be enforced
    End Subroutine

    ! **************************************************************************************************************
    !         PROCEDURES FOR ENFORCING THE DESACTIVATION OF LOGS FOR A SET OF PROCEDURES
    ! **************************************************************************************************************
    Purity Module Subroutine GetProcedureNotToLog( This, ProcName )
      class(Logger_Type)                                    ,intent(in)     ::  This                        !< Passed-object dummy argument corresponding
      character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      ProcName                    !< Names of the procedures for which the desactivation of logs should be enforced
    End Subroutine
    Purity Module Subroutine AddProcNotToLog0d( This, ProcName )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument corresponding
      character(*)                                          ,intent(in)     ::      ProcName                    !< Name of a procedure to add to the list of procedures for which the desactivation of logs should be enforced
    End Subroutine
    Purity Module Subroutine AddProcNotToLog1d( This, ProcName )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument corresponding
      character(*)  ,dimension(:)                           ,intent(in)     ::      ProcName                    !< Name of the procedures to add to the list of procedures for which the desactivation of logs should be enforced
    End Subroutine

    ! **************************************************************************************************************
    !                           PROCEDURES TO ACCESS THE PROPERTIES OF THE LOGGER OBJECT
    ! **************************************************************************************************************
    Purity Module Function Active( This, LogLevel ) result(IsActive)
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument
      integer                                     ,optional ,intent(in)     ::      LogLevel
      logical                                                               ::      IsActive
    End Function
    Purity Module Function GetIndentation( This ) result(Indentation)
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument
      integer                                                               ::      Indentation
    End Function
    Purity Module Function GetPrefix( This, LogLevel, Error, Warning, Info, Debug, HeavyDebug ) result(Prefix)
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      character(:)  ,allocatable                                            ::      Prefix
    End Function
    Purity Module Function GetPath( This ) result(Path)
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument
      character(:)  ,allocatable                                            ::      Path                            !< Path of the Logger in term of procedure names
    End Function
    Purity Module Function GetLogLevel( This ) result(LogLevel)
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument
      character(:)  ,allocatable                                            ::      LogLevel
    End Function

    Purity Module Subroutine StartDebugMode( This )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument corresponding
    End Subroutine
    Purity Module Subroutine StopDebugMode( This )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument corresponding
    End Subroutine



    ! **************************************************************************************************************
    !                           PRIVATE PROCEDURES USED DURING INITIALIZATION OR WRITING
    ! **************************************************************************************************************
    Purity Module Subroutine SetFileName( This, FileName, i_Force_FileName )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument corresponding
      character(*)                                          ,intent(in)     ::      FileName                        !< FileName of the log file
      logical                                     ,optional ,intent(in)     ::      i_Force_FileName
    End Subroutine
    Purity Module Subroutine AddItem( This, LoggerItem )
      class(Logger_Type)                ,target                    ,intent(inout)  ::  This                        !< Passed-object dummy argument corresponding the Logger object
      type(LoggerItem_Type)                                 ,intent(in)     ::      LoggerItem                    !< LoggerItem object to be added to the list of LoggerItem sub-objects
    End Subroutine
    Purity Module Subroutine RemoveItem( This )
      class(Logger_Type)                ,target                    ,intent(inout)  ::  This                        !< Passed-object dummy argument corresponding the Logger object
    End Subroutine
    Purity Module Function Set_Prefix_Indentation( This, LogLevel, Error, Warning, Info, Debug, HeavyDebug ) result(String)
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      character(:)  ,allocatable                                            ::      String                          !< Character string corresponding to the indentation part of the format's prefix
    End Function

    Purity Module Function GetPrefixProcedure( This ) result(String)
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument
      character(:)  ,allocatable                                            ::      String                          !< Character string corresponding to the "procedure" part of the format's prefix
    End Function

!     Module Function Set_Local_Prefix( This, i_Prefix ) result(Local_Prefix)
!       class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument corresponding
!       logical                                     ,optional                 ::      i_Prefix                    !< Indicator whether the prefix has to be written
!       character(:)  ,allocatable                                            ::      Local_Prefix                    !< PRefix value
!     End Function


    Module Subroutine Set_Advancing( This, String, Advance )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      character(:)  ,allocatable                            ,intent(out)    ::      String
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
    End Subroutine

    Module Subroutine Set_Backspace( This, Backspace )
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
    End Subroutine

    Module Subroutine Write_NewLine( This, NewLine )
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
    End Subroutine

    Module Subroutine Error_Open( This )
      class(Logger_Type)                                    ,intent(in)     ::  This                            !< Passed-object dummy argument corresponding
    End Subroutine

    Module Subroutine Write_Blank_Line( This )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument corresponding
    End Subroutine

    Module Subroutine CallingProcedure(  This,                                                   &
                            V1,                                                     &
                            Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Underline, Fc, Fi, Fr, Fmt, &
                            F1                                                      )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      character(*)                                          ,intent(in)     ::      V1
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      Underline
      character(*)                                ,optional ,intent(in)     ::      F1
    End Subroutine

    Module Subroutine Write_1xV0(  This,                                                   &
                            V1,                                                     &
                            Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Underline, Fc, Fi, Fr, Fmt, &
                            F1                                                      )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      Underline
      character(*)                                ,optional ,intent(in)     ::      F1
    End Subroutine

    Module Subroutine Write_2xV0(  This,                                               &
                                  V1,  V2,                                            &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2                                             )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2
    End Subroutine

    Module Subroutine Write_3xV0(  This,                                               &
                                  V1,  V2,  V3,                                       &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3                                        )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3
      character(:)  ,allocatable                                            ::      S1,  S2,  S3
    End Subroutine

    Module Subroutine Write_4xV0(   This,                                               &
                                  V1,  V2,  V3,  V4,                                  &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4                                   )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4
    End Subroutine

    Module Subroutine Write_5xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,                             &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5                              )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5
    End Subroutine

    Module Subroutine Write_6xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,                        &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6                         )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6
    End Subroutine

    Module Subroutine Write_7xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,                   &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7                    )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7
    End Subroutine

    Module Subroutine Write_8xV0(   This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,              &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8               )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8
    End Subroutine

    Module Subroutine Write_9xV0(   This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,         &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9          )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9
    End Subroutine

    Module Subroutine Write_10xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10    )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10
    End Subroutine

    Module Subroutine Write_11xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11,                                                &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11                                                 )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11
    End Subroutine

    Module Subroutine Write_12xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12,                                           &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12                                            )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12
    End Subroutine

    Module Subroutine Write_13xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13,                                      &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13                                       )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13
    End Subroutine

    Module Subroutine Write_14xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14,                                 &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14                                  )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14
    End Subroutine

    Module Subroutine Write_15xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15,                            &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15                             )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15
    End Subroutine

    Module Subroutine Write_16xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16,                       &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16                        )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16
    End Subroutine

    Module Subroutine Write_17xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17,                  &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17                   )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17
    End Subroutine

    Module Subroutine Write_18xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18,             &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18              )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18
    End Subroutine

    Module Subroutine Write_19xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19,        &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19         )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19
    End Subroutine

    Module Subroutine Write_20xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20    )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20
    End Subroutine

    Module Subroutine Write_21xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21,                                                &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21                                                 )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21
    End Subroutine

    Module Subroutine Write_22xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22,                                           &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22                                            )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22
    End Subroutine

    Module Subroutine Write_23xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23,                                      &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23                                       )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23
    End Subroutine

    Module Subroutine Write_24xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24,                                 &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24                                  )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24
    End Subroutine

    Module Subroutine Write_25xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24, V25,                            &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24, F25                             )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24, V25
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24, F25
    End Subroutine

    Module Subroutine Write_26xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24, V25, V26,                       &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24, F25, F26                        )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24, V25, V26
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24, F25, F26
    End Subroutine

    Module Subroutine Write_27xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24, V25, V26, V27,                  &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24, F25, F26, F27                   )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24, V25, V26, V27
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24, F25, F26, F27
    End Subroutine

    Module Subroutine Write_28xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24, V25, V26, V27, V28,             &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24, F25, F26, F27, F28              )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24, V25, V26, V27, V28
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24, F25, F26, F27, F28
    End Subroutine

    Module Subroutine Write_29xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24, V25, V26, V27, V28, V29,        &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24, F25, F26, F27, F28, F29         )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24, V25, V26, V27, V28, V29
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24, F25, F26, F27, F28, F29
    End Subroutine

    Module Subroutine Write_30xV0(  This,                                               &
                                  V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                  V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                  V21, V22, V23, V24, V25, V26, V27, V28, V29, V30,   &
                                  Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                  F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                  F21, F22, F23, F24, F25, F26, F27, F28, F29, F30    )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                                                                    V11, V12, V13, V14, V15, V16, V17, V18, V19, V20,   &
                                                                                    V21, V22, V23, V24, V25, V26, V27, V28, V29, V30
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10,   &
                                                                                    F11, F12, F13, F14, F15, F16, F17, F18, F19, F20,   &
                                                                                    F21, F22, F23, F24, F25, F26, F27, F28, F29, F30
    End Subroutine

    Module Subroutine Write_1xV1(   This,                                               &
                                    V1,                                                 &
                                    Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F1                                                  )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)      ,dimension(:)                           ,intent(in)     ::      V1
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1
    End Subroutine


    Module Subroutine Write_1xV0_1xV1(  This,                                               &
                                    V1,  V2,                                            &
                                    Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F1,  F2                                             )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1
      class(*)      ,dimension(:)                           ,intent(in)     ::      V2
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2
    End Subroutine

    Module Subroutine Write_2xV0_1xV1(  This,                                               &
                                    V1,  V2,  V3,                                       &
                                    Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F0, F1,  F2,  F3                                    )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2
      class(*)      ,dimension(:)                           ,intent(in)     ::      V3
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F0, F1,  F2,  F3
    End Subroutine

    Module Subroutine Write_3xV0_1xV1(  This,                                               &
                                    V1,  V2,  V3,  V4,                                  &
                                    Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F1,  F2,  F3,  F4                                   )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3
      class(*)      ,dimension(:)                           ,intent(in)     ::      V4
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4
    End Subroutine

    Module Subroutine Write_5xV0_1xV1(  This,                                               &
                                    V1,  V2,  V3,  V4,  V5,  V6,                        &
                                    Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F1,  F2,  F3,  F4,  F5,  F6                         )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5
      class(*)      ,dimension(:)                           ,intent(in)     ::      V6
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6
    End Subroutine

    Module Subroutine Write_7xV0_1xV1(  This,                                               &
                                    V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,              &
                                    Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8               )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7
      class(*)      ,dimension(:)                           ,intent(in)     ::      V8
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8
    End Subroutine

    Module Subroutine Write_9xV0_1xV1(  This,                                               &
                                    V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10,   &
                                    Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                    F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10    )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9
      class(*)      ,dimension(:)                           ,intent(in)     ::      V10
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9,  F10
    End Subroutine

    Module Subroutine Write_1xV0_2xV1(  This,                                               &
                                      V1,  V2,  V3,                                       &
                                      Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                      F0, F1,  F2,  F3                                    )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1
      class(*)      ,dimension(:)                           ,intent(in)     ::      V2, V3
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F0,  F1,  F2,  F3
    End Subroutine

    Module Subroutine Write_1xV0_2xV0V1(  This,                                             &
                                      V1,  V2,  V3,  V4,  V5,                             &
                                      Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                      F0, F1,  F2,  F3,  F4,  F5                          )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1, V2, V4
      class(*)      ,dimension(:)                           ,intent(in)     ::      V3, V5
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F0, F1,  F2,  F3,  F4,  F5
    End Subroutine

    Module Subroutine Write_1xV0_3xV0V1(  This,                                             &
                                      V1,  V2,  V3,  V4,  V5,  V6,  V7,                   &
                                      Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                                      F0, F1,  F2,  F3,  F4,  F5,  F6,  F7                )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1, V2, V4, V6
      class(*)      ,dimension(:)                           ,intent(in)     ::      V3, V5, V7
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F0, F1,  F2,  F3,  F4,  F5,  F6,  F7
    End Subroutine

    Module Subroutine Write_1xV0_1xV2(  This,                                              &
                            V1, V2,                                                 &
                            Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                            F1, F2                                                  )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1
      class(*)      ,dimension(:,:)                         ,intent(in)     ::      V2
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1, F2
    End Subroutine

    Module Subroutine Write_2xV0V1(  This,                                                   &
                              V1,  V2,  V3,  V4,                                      &
                              Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                              F1,  F2,  F3,  F4                                       )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1, V3
      class(*)      ,dimension(:)                           ,intent(in)     ::      V2, V4
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4
    End Subroutine

    Module Subroutine Write_3xV0V1(  This,                                                   &
                              V1,  V2,  V3,  V4,  V5,  V6,                            &
                              Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                              F1,  F2,  F3,  F4,  F5,  F6                             )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1, V3, V5
      class(*)      ,dimension(:)                           ,intent(in)     ::      V2, V4, V6
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6
    End Subroutine

    Module Subroutine Write_4xV0V1(  This,                                                   &
                              V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,                  &
                              Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, NewLine, Advance, Backspace, LogLevel, Status, Fc, Fi, Fr, Fmt, &
                              F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8 )
      class(Logger_Type)                                    ,intent(inout)  ::  This                        !< Passed-object dummy argument
      class(*)                                              ,intent(in)     ::      V1, V3, V5, V7
      class(*)      ,dimension(:)                           ,intent(in)     ::      V2, V4, V6, V8
      logical       ,dimension(:,:,:,:)           ,optional ,intent(in)     ::      Unused                      !< Unused optional variables required to avoid "The type/rank/keyword signature for this specific procedure matches another specific procedure that shares the same generic binding name."
      logical                                     ,optional ,intent(in)     ::      i_Prefix                    !< Indicator of the prefix presence
      integer                                     ,optional ,intent(in)     ::      LogLevel                    !< Indicator of the log level associated to current log message
      logical                                     ,optional ,intent(in)     ::      Error                       !< Indicator of an 'Error' log message
      logical                                     ,optional ,intent(in)     ::      Warning                     !< Indicator of an 'Warning' log message
      logical                                     ,optional ,intent(in)     ::      Info                        !< Indicator of an 'Info' log message
      logical                                     ,optional ,intent(in)     ::      Debug                       !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      HeavyDebug                  !< Indicator of an 'Debug' log message
      logical                                     ,optional ,intent(in)     ::      NewLine                     !< Indicator whether or not a black line must be written before writing the variables
      logical                                     ,optional ,intent(in)     ::      Advance                     !< Indicator whether or not the line should be advanced
      logical                                     ,optional ,intent(in)     ::      Backspace                   !< Indicator whether or not backpace should be used atfer the write statment
      integer                                     ,optional ,intent(out)    ::      Status                      !< Error status indicator (=0 if everthing is ok, /=0 if error)
      character(*)                                ,optional ,intent(in)     ::      Fc, Fi, Fr, Fmt             !< Type-specific format specificators for character (Fc), interger (Fi) and real (Fr) variables and common format specificator Fc
      character(*)                                ,optional ,intent(in)     ::      F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8
    End Subroutine

  End Interface

End Module