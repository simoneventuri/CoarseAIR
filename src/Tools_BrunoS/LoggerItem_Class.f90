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
!     logical                                             ::      Active      =   .True.
  contains
    procedure ,public     ::    Initialize    =>  InitializeLoggerItem
    procedure ,public     ::    Output        =>  OutputLoggerItem
    procedure ,public     ::    IsActive      =>  IsLoggerItemActive
    procedure ,public     ::    Indent        =>  IndentLoggerItem
    procedure ,public     ::    Deindent      =>  DeindentLoggerItem
    procedure ,private    ::    SetLogLevel
    procedure ,private    ::    SetMsgLogLevel
    procedure ,private    ::    SetIndentation
  End Type

  integer       ,parameter                      ::      IndentationStep          =     2

!   Interface             LoggerItem_Type
!     Module Procedure    ConstructLoggerItem
!   End Interface

  contains

    Pure Subroutine InitializeLoggerItem( This, ProcedureName, Indentation, LogLevel, DefLogLevel, MsgLogLevel )
      class(LoggerItem_Type)                                ,intent(inout)    ::      This                          !< Passed-object dummy argument
      character(*)                                          ,intent(in)     ::      ProcedureName                 !< Name of the calling procedure
      integer                                     ,optional ,intent(in)     ::      Indentation
      integer                                     ,optional ,intent(in)     ::      LogLevel
      integer                                     ,optional ,intent(in)     ::      DefLogLevel
      integer                                     ,optional ,intent(in)     ::      MsgLogLevel
  This%Name   =   trim( ProcedureName )                                                                   ! Setting the name of the Item to the input procedure name
  call This%SetLogLevel( LogLevel, DefLogLevel )                                                          ! Setting the log level of current item to {OptionalLogLevel,DefaultLogLevel,LogLevel_DEFAULT} in this order
  call This%SetMsgLogLevel( MsgLogLevel )                                                                 ! Setting the default log level of the log message written from this object to {LogLevel_DEFAULT,MsgLogLevel} in this order
  call This%SetIndentation( Indentation )                                                                 ! Setting the default log level of the log message written from this object to {LogLevel_DEFAULT,MsgLogLevel} in this order
    End Subroutine
    Pure Function OutputLoggerItem( This ) result(String)
      class(LoggerItem_Type)                               ,intent(in)     ::      This
      character(:)    ,allocatable                                          ::      String
  character(1000) ::  LongString
  write(LongString,"('Index = ',i3,3x,'Name = ',a12,3x,'Indentation = ',i3,3x,'This%LogLevel = ',i1,3x,'MsgLogLevel = ',i1)") This%Index, This%Name, This%Indentation, This%LogLevel, This%MsgLogLevel
  String    =   trim(LongString)
    End Function


    Pure Subroutine IndentLoggerItem( This )
      class(LoggerItem_Type)                                ,intent(inout)  ::      This
  This%Indentation = This%Indentation + IndentationStep
    End Subroutine
    Pure Subroutine DeindentLoggerItem( This )
      class(LoggerItem_Type)                                ,intent(inout)  ::      This
  This%Indentation = This%Indentation - IndentationStep
    End Subroutine



    Pure Subroutine SetLogLevel( This, OptionalLogLevel, DefaultLogLevel )
      class(LoggerItem_Type)                                ,intent(inout)  ::      This
      integer                                     ,optional ,intent(in)     ::      OptionalLogLevel
      integer                                     ,optional ,intent(in)     ::      DefaultLogLevel
  This%LogLevel   =   3!LogLevel_DEFAULT
  if ( present(OptionalLogLevel) ) then
    This%LogLevel =   OptionalLogLevel
  else
    if ( present(DefaultLogLevel) ) This%LogLevel = DefaultLogLevel
  end if
  This%SavedLogLevel  =   This%LogLevel
    End Subroutine
    Pure Subroutine SetMsgLogLevel( This, OptionalValue )
      class(LoggerItem_Type)                                ,intent(inout)  ::      This
      integer                                     ,optional ,intent(in)     ::      OptionalValue
  This%MsgLogLevel   =   3  !LogLevel_DEFAULT                                                                   ! Setting the default message log level to the default value
  if ( present(OptionalValue) ) This%MsgLogLevel = OptionalValue                                                ! Setting the default message log level to the optional value if present
  This%SavedMsgLogLevel  =   This%MsgLogLevel                                                                   ! Saving the default message log level
    End Subroutine
    Pure Subroutine SetIndentation( This, Indentation )
      class(LoggerItem_Type)                                ,intent(inout)  ::      This                          !< Passed-object dummy argument
      integer                                     ,optional ,intent(in)     ::      Indentation
  This%Indentation      =   0
  if ( present(Indentation) ) This%Indentation = Indentation
    End Subroutine



    Function IsLoggerItemActive( This, LogLevel, ProcToLog ) result(Active)
      class(LoggerItem_Type)                                ,intent(in)     ::      This                            !< Passed-object dummy argument
      integer                                     ,optional ,intent(in)     ::      LogLevel
      character(:)  ,allocatable ,dimension(:)    ,optional ,intent(in)     ::      ProcToLog
      logical                                                               ::      Active
  integer                                                               ::      LogLevel_
  LogLevel_   =   This%MsgLogLevel    ! Setting the message log level to the default stored in the Logger Item
  if ( present(LogLevel) ) LogLevel_ = LogLevel
  Active    =   ( LogLevel_ <= This%LogLevel )
  if ( present(ProcToLog) ) then
    if ( IsPresent( This%Name, ProcToLog ) ) Active = .True.
  end if
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

End Module