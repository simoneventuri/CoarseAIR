SubModule(LoggerItem_Class) LoggerItem_SubClass

  implicit none

  integer       ,parameter                      ::      IndentationStep          =     2

  contains

! ! @TODO: I think the optional could be removed (safer)
! Module Procedure ConstructLoggerItem
!   LoggerItem%Name         =   ''
!   LoggerItem%Indentation  =   0
!   LoggerItem%Active       =   .True.
!   if ( present(Name) )        LoggerItem%Name         =   Name
!   if ( present(Indentation) ) LoggerItem%Indentation  =   Indentation
!   if ( present(Active) )      LoggerItem%Active       =   Active
!   if ( present(LogLev) )        LoggerItem%LogLevel         =   LogLev
! End Procedure

Module Procedure InitializeLoggerItem
  call This%Free()
  This%Name   =   trim( ProcedureName )                                                                   ! Setting the name of the Item to the input procedure name
  call This%SetLogLevel( LogLevel, DefLogLevel )                                                          ! Setting the log level of current item to {OptionalLogLevel,DefaultLogLevel,LogLevel_DEFAULT} in this order
  call This%SetMsgLogLevel( MsgLogLevel )                                                                 ! Setting the default log level of the log message written from this object to {LogLevel_DEFAULT,MsgLogLevel} in this order
  call This%SetIndentation( Indentation )                                                                 ! Setting the default log level of the log message written from this object to {LogLevel_DEFAULT,MsgLogLevel} in this order
End Procedure

Module Procedure FreeLoggerItem
  if ( allocated(This%Name) ) deallocate(This%Name)
  This%Index             =   0
  This%Indentation       =   0
  This%LogLevel          =   3
  This%SavedLogLevel     =   3
  This%MsgLogLevel       =   3
  This%SavedMsgLogLevel  =   3
End Procedure

Module Procedure FinalizeLoggerItem
  call This%Free()
End Procedure

! This procedure sets the indentation level.
Module Procedure SetIndentation
  This%Indentation      =   0
  if ( present(Indentation) ) This%Indentation = Indentation
!   if ( This%IsActive() ) This%Indentation = This%Indentation + IndentationStep
End Procedure

Module Procedure IndentLoggerItem
  This%Indentation = This%Indentation + IndentationStep
End Procedure
Module Procedure DeindentLoggerItem
  This%Indentation = This%Indentation - IndentationStep
End Procedure

! @TODO: I think the optional could be removed (safer)
Module Procedure OutputLoggerItem
  character(1000) ::  LongString
  write(LongString,"('Index = ',i3,3x,'Name = ',a12,3x,'Indentation = ',i3,3x,'This%LogLevel = ',i1,3x,'MsgLogLevel = ',i1)") This%Index, This%Name, This%Indentation, This%LogLevel, This%MsgLogLevel
  String    =   trim(LongString)
End Procedure

! This procedure sets the log levels associated to a given 'LoggerItem' object.
! The log level is set as followed:
! * to the 'OptionalLogLevel' if present or,
! * to the 'DefaultLogLevel_' value if present and if the 'OptionalLogLevel' argument is absent or,
! * to the 'LogLevel_DEFAULT' value if both the 'OptionalLogLevel' and 'DefaultLogLevel_' arguemnts are absent.
! The value of the 'LogLevel' component is then saved in the component 'SavedLogLevel'
Pure Module Subroutine SetLogLevel( This, OptionalLogLevel, DefaultLogLevel )
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

! This procedure sets the default log levels associated to the message written from a 'LoggerItem' object.
! The log level is set as followed:
! * to the 'OptionalValue' if present or,
! * to the 'LogLevel_DEFAULT' if the 'OptionalValue' argument is absent.
! The value of the 'MsgLogLevel' component is then saved in the component 'SavedMsgLogLevel'
Module Procedure SetMsgLogLevel
  This%MsgLogLevel   =   3  !LogLevel_DEFAULT                                                                   ! Setting the default message log level to the default value
  if ( present(OptionalValue) ) This%MsgLogLevel = OptionalValue                                                ! Setting the default message log level to the optional value if present
  This%SavedMsgLogLevel  =   This%MsgLogLevel                                                                   ! Saving the default message log level
End Procedure


Module Procedure IsLoggerItemActive
  integer                                                               ::      LogLevel_
  LogLevel_   =   This%MsgLogLevel    ! Setting the message log level to the default stored in the Logger Item
  if ( present(LogLevel) ) LogLevel_ = LogLevel
  Active    =   ( LogLevel_ <= This%LogLevel )
  if ( present(ProcToLog) ) then
    if ( IsPresent( This%Name, ProcToLog ) ) Active = .True.
  end if
End Procedure


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

End SubModule