Module LoggerItem_Class

  implicit none

  private
  public        ::      LoggerItem_Type

  Type                                                  ::      LoggerItem_Type
    character(:)        ,allocatable                    ::      Name
    integer                                             ::      Index             =   0
    integer                                             ::      Indentation       =   0
    integer                                             ::      LogLevel          =   3   ! Log-level associated to current 'LoggerItem' object (Default=Info=3)
    integer                                             ::      SavedLogLevel     =   3   ! Copie of the Saved value of the 'Mode' (Used to restore)
    integer                                             ::      MsgLogLevel       =   3   ! Default log-level for all logging message performed by current current 'LoggerItem' object (Default=Info=3)
    integer                                             ::      SavedMsgLogLevel  =   3   ! Default log-level for all logging message performed by current current 'LoggerItem' object (Default=Info=3)
  contains
    private
    Final                 ::    FinalizeLoggerItem
    procedure ,public     ::    Initialize    =>  InitializeLoggerItem
    procedure ,public     ::    Free          =>  FreeLoggerItem
    procedure ,public     ::    Output        =>  OutputLoggerItem
    procedure ,public     ::    IsActive      =>  IsLoggerItemActive
    procedure ,public     ::    Indent        =>  IndentLoggerItem
    procedure ,public     ::    Deindent      =>  DeindentLoggerItem
    procedure             ::    SetLogLevel
    procedure             ::    SetMsgLogLevel
    procedure             ::    SetIndentation
  End Type

!   Interface             LoggerItem_Type
!     Module Procedure    ConstructLoggerItem
!   End Interface

  Interface
    Pure Module Subroutine InitializeLoggerItem( This, ProcedureName, Indentation, LogLevel, DefLogLevel, MsgLogLevel )
      class(LoggerItem_Type)                                ,intent(inout)  ::      This                          !< Passed-object dummy argument
      character(*)                                          ,intent(in)     ::      ProcedureName                 !< Name of the calling procedure
      integer                                     ,optional ,intent(in)     ::      Indentation
      integer                                     ,optional ,intent(in)     ::      LogLevel
      integer                                     ,optional ,intent(in)     ::      DefLogLevel
      integer                                     ,optional ,intent(in)     ::      MsgLogLevel
    End Subroutine
    Pure Module Subroutine FreeLoggerItem( This )
      class(LoggerItem_Type)                                ,intent(inout)  ::      This                          !< Passed-object dummy argument
    End Subroutine
    Pure Module Subroutine FinalizeLoggerItem( This )
      type(LoggerItem_Type)                                 ,intent(inout)  ::      This
    End Subroutine
    Pure Module Function OutputLoggerItem( This ) result(String)
      class(LoggerItem_Type)                                ,intent(in)     ::      This
      character(:)    ,allocatable                                          ::      String
    End Function
    Pure Module Subroutine IndentLoggerItem( This )
      class(LoggerItem_Type)                                ,intent(inout)  ::      This
    End Subroutine
    Pure Module Subroutine DeindentLoggerItem( This )
      class(LoggerItem_Type)                                ,intent(inout)  ::      This
    End Subroutine
    Pure Module Subroutine SetLogLevel( This, OptionalLogLevel, DefaultLogLevel )
      class(LoggerItem_Type)                                ,intent(inout)  ::      This
      integer                                     ,optional ,intent(in)     ::      OptionalLogLevel
      integer                                     ,optional ,intent(in)     ::      DefaultLogLevel
    End Subroutine
    Pure Module Subroutine SetMsgLogLevel( This, OptionalValue )
      class(LoggerItem_Type)                                ,intent(inout)  ::      This
      integer                                     ,optional ,intent(in)     ::      OptionalValue
    End Subroutine
    Pure Module Subroutine SetIndentation( This, Indentation )
      class(LoggerItem_Type)                                ,intent(inout)  ::      This                          !< Passed-object dummy argument
      integer                                     ,optional ,intent(in)     ::      Indentation
    End Subroutine
    Pure Module Function IsLoggerItemActive( This, LogLevel, ProcToLog ) result(Active)
      class(LoggerItem_Type)                                ,intent(in)     ::      This                            !< Passed-object dummy argument
      integer                                     ,optional ,intent(in)     ::      LogLevel
      character(:)  ,allocatable ,dimension(:)    ,optional ,intent(in)     ::      ProcToLog
      logical                                                               ::      Active
    End Function
  End Interface

End Module