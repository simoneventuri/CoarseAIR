SubModule(Logger_Class) Logger_SubClass

! #define Purity Pure
#define Purity

  implicit none

  integer       ,parameter                      ::      Indentation_Step          =   2
  character(*)  ,parameter                      ::      Default_Name_LogFile      =   "logfile.log"
  character(*)  ,parameter                      ::      Default_Logical_Format    =   "l1"
  character(*)  ,parameter                      ::      Default_Integer_Format    =   "i0"
  character(*)  ,parameter                      ::      Default_Real_Format       =   "g0"
  character(*)  ,parameter                      ::      Default_Character_Format  =   "a"
  character(*)  ,parameter                      ::      Default_Spacing2_Format   =   "3x"
  logical       ,parameter                      ::      Default_WriteInOut        =   .True.

  character(*)  ,parameter                      ::      Default_OpenStatus   = "REPLACE"
  character(*)  ,parameter                      ::      Default_OpenPosition = "REWIND"
  character(*)  ,parameter      ,dimension(2)   ::      Valid_OpenStatus     = ["REPLACE","OLD    "]
  character(*)  ,parameter      ,dimension(2)   ::      Valid_OpenPosition   = ["REWIND","APPEND"]

  logical ,parameter :: LocalDebug=.False.

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
    Module Procedure    Convert_Logical_To_String_0D, Convert_Logical_To_String_1D, Convert_Logical_To_Strings_1D
    Module Procedure    Convert_INT8_To_String_0D   , Convert_INT8_To_String_1D   , Convert_INT8_To_Strings_1D
    Module Procedure    Convert_INT16_To_String_0D  , Convert_INT16_To_String_1D  , Convert_INT16_To_Strings_1D
    Module Procedure    Convert_INT32_To_String_0D  , Convert_INT32_To_String_1D  , Convert_INT32_To_Strings_1D
    Module Procedure    Convert_INT64_To_String_0D  , Convert_INT64_To_String_1D  , Convert_INT64_To_Strings_1D
    Module Procedure    Convert_REAL32_To_String_0D , Convert_REAL32_To_String_1D , Convert_REAL32_To_Strings_1D
    Module Procedure    Convert_REAL64_To_String_0D , Convert_REAL64_To_String_1D , Convert_REAL64_To_Strings_1D
    Module Procedure    Convert_REAL128_To_String_0D, Convert_REAL128_To_String_1D, Convert_REAL128_To_Strings_1D
    Module Procedure    Convert_String_To_String_0D , Convert_String_To_String_1D , Convert_String_To_Strings_1D
  End Interface

  contains

! This procedure initializes a Logger object.
Module Procedure InitializeLogger
  character(:)  ,allocatable                                            ::      OpenStatus
  character(:)  ,allocatable                                            ::      OpenPosition
  integer                                                               ::      ios
  call This%Free()
  call This%SetFileName( FileName, i_Force_FileName )                                                           ! Setting the name of the file associated to the Logger object
  call This%Activate()
  call SetInitialLoggerItem( This, Procedure, Indentation )
  OpenStatus   =       Get_OptOrDef_Value( Default_OpenStatus,   Valid_OpenStatus,   Status   )                 ! Getting the open status
  OpenPosition =       Get_OptOrDef_Value( Default_OpenPosition, Valid_OpenPosition, Position )                 ! Getting the open position
  open( NewUnit=This%Unit, File=This%FileName, Status=OpenStatus, Position=OpenPosition, IOStat=ios )           ! Opening the log file
  if ( ios /= 0 ) call This%Error_Open()                                                                        ! If opening error, then print error message and stop the code
End Procedure

Module Procedure Activate
  This%Activated  =   .True.
End Procedure
Module Procedure Deactivate
  This%Activated  =   .False.
End Procedure

! This procedure free a Logger object.
Module Procedure FreeLogger
  call FinalizeLogger(This)
End Procedure

Module Procedure FinalizeLogger
  This%Unit       =   Output_Unit
  This%Advancing  =   .True.
  This%Activated  =   .True.
  This%iLogLev    =   0
  This%NLevels    =   0
  if ( allocated(This%FileName) )     deallocate(This%FileName)
  if ( allocated(This%ProcToLog) )    deallocate(This%ProcToLog)
  if ( allocated(This%ProcNotToLog) ) deallocate(This%ProcNotToLog)
  if ( allocated(This%Item) )         deallocate(This%Item)
!   if ( allocated(This%CurrentItem) )  deallocate(This%CurrentItem)
  if ( associated(This%CurrentItem) ) nullify(This%CurrentItem)
End Procedure

! This procedure reopens a Logger object.
! The current Logger is first closed and deleted. Then, it is initialized using the same filename than the
! one it currently has. This will open a new file (with the same name but eventually in a different directory)
! with the poition is set to 'APPEND'.
! This procedure is mainly used when one wnat to change the working directory of the Fortran application.
! Once the application has changed directory, this 'Reopen' procedure is called to delete the logfile in the
! previous working directory, and to open it in the new working directory. The 'APPEND' position ensures
! that the informaton previously written to the logfile are not overwritten.
Module Procedure ReopenLogger
  character(*)  ,parameter                                              ::      OpenStatus   = 'OLD'
  character(*)  ,parameter                                              ::      OpenPosition = 'APPEND'
  integer                                                               ::      ios
  call This%Close( Status='DELETE')                                                                             ! Deleting the current Logger (the log file is deleted)
  open( NewUnit=This%Unit, File=This%FileName, Status=OpenStatus, Position=OpenPosition, IOStat=ios ) ! Opening the log file
  if ( ios /= 0 ) call This%Error_Open()                                                                ! If opening error, then print error message and stop the code
End Procedure

! This procedure closes a Logger object.
Module Procedure CloseLogger
  character(:)  ,allocatable                                            ::      Close_Status                    !< Local value of the close status
  Close_Status  =       ""
!   TODO: Add a checking that the "Status" string is a valid close status
  if ( present(Status) ) then
    Close_Status        =       Status
    close( This%Unit, Status=Close_Status )
  else
    close( This%Unit )
  end if
End Procedure

! This procedure backspaces the logfile associated to a Logger object.
Module Procedure BackspaceLogger
  integer                                                               ::      ios
  backspace( This%Unit, iostat=ios )
  if ( present(Status) ) Status = ios
End Procedure

! This procedure rewinds the logfile associated to a Logger object
Module Procedure RewindLogger
  integer                                                               ::      ios
  rewind( This%Unit, iostat=ios )
  if ( present(Status) ) Status = ios
End Procedure

! This procedure flushes the unit associated to a Logger object
Module Procedure FlushLogger
  flush( This%Unit )
End Procedure

Module Procedure Entering
  logical                                                               ::      WriteInOut
  integer                                                               ::      Indentation
  type(LoggerItem_Type)                                                 ::      LoggerItem
  if ( .Not.This%Activated ) return
  Indentation   =   This%GetIndentation()                                                                       ! Getting the current indentation levels
  call LoggerItem%Initialize( ProcedureName, Indentation, LogLevel, DefLogLevel, MsgLogLevel )                  ! Initializing a 'LoggerItem' object at current indentation level and with the input properties
  if ( LoggerItem%IsActive(ProcToLog=This%ProcToLog) ) call LoggerItem%Indent()                                 ! If the new 'LoggerItem' object is active, then increase its indentation level
  call This%AddItem( LoggerItem )                                                                               ! Adding the new 'LoggerItem' object to the list of 'LoggerItem' sub-objects
  if (present(Writing)) then; WriteInOut = Writing; else; WriteInOut = Default_WriteInOut; end if
  if ( WriteInOut .and. This%Active() ) call This%Write("Entering")                                             ! Writing the "Entering" message
End Procedure

Module Procedure StartDebugMode
  if ( .Not. associated(This%CurrentItem) ) return
  This%CurrentItem%SavedMsgLogLevel  =   This%CurrentItem%MsgLogLevel
  This%CurrentItem%MsgLogLevel       =   LogLevel_DEBUG

!   if ( .Not. allocated(This%Item) ) return
!   associate( LogLev => This%Item(This%iLogLev) )
!     LogLev%SavedLogLevel  =   LogLev%LogLevel
!     LogLev%LogLevel       =   LogLevel_DEBUG
!   end associate
End Procedure

Module Procedure StopDebugMode
  if ( .Not. allocated(This%Item) ) return
  associate( LogLev => This%Item(This%iLogLev) )
    LogLev%MsgLogLevel       =   LogLev%SavedMsgLogLevel
  end associate
End Procedure

Module Procedure Exiting
  logical                                                               ::      WriteInOut
  if ( .Not.This%Activated ) return
  if (present(Writing)) then; WriteInOut = Writing; else; WriteInOut = Default_WriteInOut; end if
  if ( WriteInOut .and. This%Active() ) call This%Write("Exiting")                                              ! Writing the "Exiting" message
  call This%RemoveItem()                                                                                        ! Removing the last procedure from the list of nested procedures
End Procedure

Function IsCurrentActive( This, ProcedureName, OptionalActive, DefaultActive ) result(WriteLogs_)
  class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument corresponding to the Logger object
  character(*)                                          ,intent(in)     ::      ProcedureName                   !< Name of the calling procedure
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


!   if ( allocated(This%ProcToLog) ) then
!     do i = 1,size(This%ProcToLog)
!       if ( .Not. ( LowerCase(ProcedureName) == LowerCase(This%ProcToLog(i)) ) ) cycle
!       LogLevel  = .True.
!       exit
!     end do
!   end if

!   if ( allocated(This%ProcToLog) ) then
!     if ( Procedure == This%ProcToLog(1) ) then ! @TODO
!       ActiveLoggerItem  = .True.
!     end if
!   end if

! Function ActivateLogs( This, ProcedureName, OptionalActive, DefaultActive ) result(WriteLogs_)
!   class(Logger_Type)                                    ,intent(in)     ::      This                            !< Passed-object dummy argument corresponding to the Logger object
!   character(*)                                          ,intent(in)     ::      ProcedureName                       !< Name of the calling procedure
!   logical                                     ,optional ,intent(in)     ::      OptionalActive
!   logical                                     ,optional ,intent(in)     ::      DefaultActive
!   logical                                                               ::      WriteLogs_
!
!
!   logical                                                               ::      OptionalLogLevel_
!   logical                                                               ::      DefaultLogLevel_
!   integer                                                               ::      i
!
!   DefaultLogLevel_   = Set_Optional_Argument( LogLevel_DEFAULT, DefaultLogLevel  )
!   OptionalLogLevel_  = Set_Optional_Argument( DefaultLogLevel_, OptionalLogLevel )
!
!   if ( allocated(This%ProcToLog) ) then
!     do i = 1,size(This%ProcToLog)
!       if ( .Not. ( LowerCase(ProcedureName) == LowerCase(This%ProcToLog(i)) ) ) cycle
!       LogLevel  = .True.
!       exit
!     end do
!   end if
! !   if ( allocated(This%ProcToLog) ) then
! !     if ( Procedure == This%ProcToLog(1) ) then ! @TODO
! !       ActiveLoggerItem  = .True.
! !     end if
! !   end if
! End Function

! **************************************************************************************************************
!              PROCEDURES FOR ENFORCING THE ACTIVATION OF LOGS FOR A SET OF PROCEDURES
! **************************************************************************************************************
Module Procedure GetProcedureToLog
#ifdef WORKAROUND_GFORTRAN_SOURCE_ALLOCATION
  if ( allocated(This%ProcToLog) ) then
    allocate( ProcName(size(This%ProcToLog)), source = This%ProcToLog )
  end if
#else
  allocate( ProcName, source = This%ProcToLog )
#endif
End Procedure

Module Procedure AddProcToLog0d
  call AddElementToArray( ProcName, This%ProcToLog )
End Procedure

Module Procedure AddProcToLog1d
  call AddElementToArray( ProcName, This%ProcToLog )
End Procedure

! **************************************************************************************************************
!              PROCEDURES FOR ENFORCING THE DESACTIVATION OF LOGS FOR A SET OF PROCEDURES
! **************************************************************************************************************
Module Procedure GetProcedureNotToLog
#ifdef WORKAROUND_GFORTRAN_SOURCE_ALLOCATION
  if ( allocated(This%ProcNotToLog) ) then
    allocate( ProcName(size(This%ProcNotToLog)), source = This%ProcNotToLog )
  end if
#else
  allocate( ProcName, source = This%ProcNotToLog )
#endif
End Procedure

Module Procedure AddProcNotToLog0d
  call AddElementToArray( ProcName, This%ProcNotToLog )
End Procedure

Module Procedure AddProcNotToLog1d
  call AddElementToArray( ProcName, This%ProcNotToLog )
End Procedure

! **************************************************************************************************************
!                   PROCEDURES TO ACCESS THE PROPERTIES OF THE LOGGER OBJECT
! **************************************************************************************************************

! This procedure return a logical variable which indicates whether the 'Logger' is active (ie. it will should
! logs) for the default log-level of logging message of the current 'LoggerItem' level.
! If the 'Logger' object has no 'LoggerItem' sub-objects (ie. if the component 'CurrentItem' is not associated),
! then the Logger is always set as active. This behavior is required in order to write logs without having
! initialized the 'Logger' object through a call to 'Logger%Initialize'.
! A given log level can be passed as an optional input arguemnt. If so, the returned variable corresponds to a
! logical variable which indicates if the Logger is active for that particuler log level.
! @TODO: Comment: Note that ... override the Active indeicator if current procedure corresponds to a procedure
! which need to be logged
Module Procedure Active
  IsActive  =   .True.
  if ( .Not. This%Activated ) then
    IsActive  =   .False.
    return
  end if
  if ( associated(This%CurrentItem) ) IsActive = This%CurrentItem%IsActive( LogLevel, This%ProcToLog  )
End Procedure

! This procedure returns a character string corresponding to the prefix of the format used to write logs.
! This prefix is made of the indentation level and the name of the procedure.
Module Procedure GetPrefix
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
End Procedure

! This procedure returns a character string corresponding to the path of the Logger in term of procedure names.
! This path corresponds to the list of all procedures the Logger went through.
! Each name are separated by the character
! @TODO: Define a optional input arguement to select the Separator
Module Procedure GetPath
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
End Procedure

! This procedure returns a character string which corresponds to the name of the log level of the current
! LoggerItem object.
Module Procedure GetLogLevel
  LogLevel = ""
  if ( associated (This%CurrentItem) ) LogLevel = LogLevelToString( This%CurrentItem%LogLevel )
!   if ( .Not. allocated(This%Item) ) return
!   LogLevel     =   LogLevelToString( This%Item(This%iLogLev)%LogLevel )
End Procedure

! This procedure returns an integer which corresponds to the indentation level of the current LoggerItem object.
Module Procedure GetIndentation
  Indentation   =   0
  if ( associated (This%CurrentItem) ) Indentation = This%CurrentItem%Indentation
End Procedure


! **************************************************************************************************************
! **************************************************************************************************************
!                                   PRIVATE PROCEDURES USED TO SET THE LOGGER PARAMETERS
! **************************************************************************************************************
! **************************************************************************************************************

! This procedure sets the name of the file associated to the logger.
! If a coarray simulation is considered, then the index of the image is added as a suffix to the filename.
! This ensures that all images have different filenames. It is required since different images are not
! allowed to operate on the same file.
Module Procedure SetFileName
#ifdef COARRAY
  character(10)                                                         ::      String                          ! Character string required to store the current image index (Required because write cannot make a length allocation)
#endif
  logical                                                               ::      i_Add_Image_Index
  i_Add_Image_Index     =       .True.
  This%FileName         =       trim(adjustl(FileName))                                                         ! Setting the FileName of the log file
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
End Procedure

! This procedure adds an element to the list of 'LoggerItem' sub-objects and update the
! current index of log-level. Also, the number of 'LoggerItem' sub-objects is updated.
Module Procedure AddItem
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
End Procedure

! This procedure removes the last element from the list of 'LoggerItem' sub-objects and update
! the current index of log-level. Also, the number of 'LoggerItem' sub-objects is updated.
Module Procedure RemoveItem
  type(LoggerItem_Type)  ,dimension(:)   ,allocatable                  ::      LoggerItems
  if ( .not. allocated(This%Item) ) allocate( This%Item(0) )
  This%NLevels      =   size(This%Item)
  allocate( LoggerItems, source=This%Item(1:This%NLevels-1) )
  call move_alloc( LoggerItems, This%Item )
  This%iLogLev      =   This%iLogLev - 1
  This%NLevels      =   size(This%Item)
  if ( (This%iLogLev>0) .and. (This%iLogLev<=This%NLevels)  ) then
    This%CurrentItem  =>  This%Item(This%iLogLev)
  else
    This%CurrentItem  => null()
  end if
End Procedure



! Module Procedure SetFileName
! #ifdef COARRAY
!   character(10)                                                         ::      String                          ! Character string required to store the current image index (Required because write cannot make a length allocation)
! #endif
!   logical                                                               ::      i_Add_Image_Index
!   i_Add_Image_Index     =       .True.
!   This%FileName         =       trim(FileName)                                                                  ! Setting the FileName of the log file
!   if ( present(i_Force_FileName) ) i_Add_Image_Index = .not. i_Force_FileName
! #ifdef COARRAY
!   if (i_Add_Image_Index) then
!     if ( This_Image() == 1 ) then                                                                                 ! If the 1st image is considered
!       This%FileName     =       trim(FileName)                                                                  ! Setting the FileName of the log file
!     else                                                                                                          ! If an image other than the 1st one is considered, then adding the image index to the log file
!       write(String,"(i0)") This_Image()                                                                           ! Converting the image index into a string
!       This%FileName     =       trim(FileName) // "_" // trim(adjustl(String))                                  ! Setting the FileName of the log file with the image index
!     end if                                                                                                        ! End if case on image index
!   end if
! #endif
! End Procedure


! This procedure sets the "indentation" part of the format's prefix.
Module Procedure Set_Prefix_Indentation
  character(:)  ,allocatable                                            ::      Format, PrefixMode
  character(10)                                                         ::      Long_String
  integer                                                               ::      Indentation
  String       =       ""
  if ( This%iLogLev > 0 ) then
  associate( LogLev => This%Item(This%iLogLev) )
    if ( LogLev%Indentation > 0) then
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
      write(Long_String,Format) Indentation
      String     =       PrefixMode // trim(adjustl(Long_String)) // "x,"
    end if
  end associate
  end if
End Procedure

! This procedure sets the "procedure" part of the format's prefix.
Module Procedure GetPrefixProcedure
  String       =       ""
  if ( associated(This%CurrentItem) ) then
    if ( allocated(This%CurrentItem%Name) ) then
      String    =       "'[" // trim(This%CurrentItem%Name) // "]: '"
    end if
  end if
End Procedure

Module Procedure Set_Advancing
  String              =     "YES"
  This%Advancing      =     .True.
  if ( present(Advance) ) then
    if (.Not.Advance) then
      String          =     "NO"
      This%Advancing  =     .False.
    end if
  end if
End Procedure

Module Procedure Set_Backspace
  if ( .Not. present(Backspace) ) return
  if ( .Not. Backspace ) return
  call This%Backspace()
End Procedure

Module Procedure Write_NewLine
  if ( present(NewLine) ) then
    if (NewLine) write(This%Unit,"(a)") ''
  end if
End Procedure

Module Procedure Error_Open
  write(*,"(4x,'[Initialize_Logger]: Error opening the Log file')")
  write(*,"(4x,'[Initialize_Logger]: FileName   = ',a)")  This%FileName
  write(*,"(4x,'[Initialize_Logger]: Unit   = ',i0)") This%Unit
  write(*,"(4x,'[Initialize_Logger]: Stopping the code')")
  error stop
End Procedure




! **************************************************************************************************************
! **************************************************************************************************************
!                                   WRITING PROCEDURES
! **************************************************************************************************************
! **************************************************************************************************************

Module Procedure Write_Blank_Line
  write(This%Unit,*)
End Procedure

Module Procedure CallingProcedure
  character(:)  ,allocatable                                            ::      Message
  if ( .Not.This%Activated ) return
  Message =   "Calling " // V1
  call This%Write( Message,                                                   &
                   Unused, i_Prefix, Error, Warning, Info, Debug, HeavyDebug, &
                   NewLine, Advance, Backspace, LogLevel, Status, Underline,  &
                   Fc, Fi, Fr, Fmt, F1                                        )
End Procedure

Module Procedure Write_1xV0
  character(:)  ,allocatable                                            ::      S1
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
!   if ( .Not. This%Active(LogLevel) ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_2xV0
  character(:)  ,allocatable                                            ::      S1,  S2
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return

!   if ( .Not. This%Active(LogLevel) ) return

  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_3xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_4xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_5xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_6xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_7xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_8xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_9xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_10xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_11xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_12xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_13xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_14xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_15xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_16xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_17xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_18xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_19xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_20xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_21xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_22xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_23xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_24xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_25xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24, S25
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24, S25
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_26xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24, S25, S26
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24, S25, S26
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_27xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24, S25, S26, S27
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24, S25, S26, S27
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_28xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24, S25, S26, S27, S28
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24, S25, S26, S27, S28
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_29xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24, S25, S26, S27, S28, S29  ! Character strings corresponding to the input variable
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24, S25, S26, S27, S28, S29
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_30xV0
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                                                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                                                S21, S22, S23, S24, S25, S26, S27, S28, S29, S30  ! Character strings corresponding to the input variable
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10,   &
                                                   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,   &
                                                   S21, S22, S23, S24, S25, S26, S27, S28, S29, S30
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_1xV1
  character(:)  ,allocatable                                            ::      S1
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_1xV0_1xV1
  character(:)  ,allocatable                                            ::      S1,  S2
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_2xV0_1xV1
  character(:)  ,allocatable                                            ::      S1,  S2, Fidx
  character(:)  ,allocatable  ,dimension(:)                             ::      S3
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Get_Integer_Format( size(V3), F0, Fidx, Status=ios );      if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
!                                             S1   i                                      S2  S3
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "a,"//Fidx//"," // Default_Spacing2_Format // ",a,a))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V3)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, i, S2, S3(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_3xV0_1xV1
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_5xV0_1xV1
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_7xV0_1xV1
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_9xV0_1xV1
  character(:)  ,allocatable                                            ::      S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
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
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_1xV0_2xV1
  character(:)  ,allocatable                                            ::      S1, Fidx
  character(:)  ,allocatable  ,dimension(:)                             ::      S3,  S2
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Get_Integer_Format( size(V3), F0, Fidx, Status=ios );      if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
!                                             S1   i                                      S2  S3
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "a,"//Fidx//"," // Default_Spacing2_Format // ",a,a))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V3)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, i, S2(i), S3(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_1xV0_2xV0V1
  character(:)  ,allocatable                                            ::      S1, S2, S4, Fidx
  character(:)  ,allocatable  ,dimension(:)                             ::      S3,  S5
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Get_Integer_Format( size(V3), F0, Fidx, Status=ios );      if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "a,"//Fidx//",*(" // Default_Spacing2_Format // ",a,a)))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V3)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, i, S2, S3(i), S4, S5(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_1xV0_3xV0V1
  character(:)  ,allocatable                                            ::      S1, S2, S4, S6, Fidx
  character(:)  ,allocatable  ,dimension(:)                             ::      S3,  S5,  S7
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Get_Integer_Format( size(V3), F0, Fidx, Status=ios );      if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "a,"//Fidx//",*(" // Default_Spacing2_Format // ",a,a)))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V3)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, i, S2, S3(i), S4, S5(i), S6, S7(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_1xV0_1xV2
  character(:)  ,allocatable                                            ::      S1
  character(:)  ,allocatable  ,dimension(:)                             ::      S2
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  character(:)  ,allocatable                                            ::      FmtStr                    ! Local prefix variable
  character(:)  ,allocatable  ,dimension(:)                             ::      LinePrefix
  integer                                                               ::      Local_Status, ios, i
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  allocate( Character(len(S1)) :: LinePrefix(size(S2,1)) )
  LinePrefix(:)   =   ""
  LinePrefix(1)   =   S1
  do i = 1,size(S2,1)
    write( This%Unit,LineFormat, Advance=Adv, IOStat=ios ) LinePrefix(i) // S2(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_2xV0V1
  character(:)  ,allocatable                                            ::      S1, S3
  character(:)  ,allocatable  ,dimension(:)                             ::      S2, S4
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V2)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, S2(i), S3, S4(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_3xV0V1
  character(:)  ,allocatable                                            ::      S1, S3, S5
  character(:)  ,allocatable  ,dimension(:)                             ::      S2, S4, S6
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V2)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, S2(i), S3, S4(i), S5, S6(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

Module Procedure Write_4xV0V1
  character(:)  ,allocatable                                            ::      S1, S3, S5, S7
  character(:)  ,allocatable  ,dimension(:)                             ::      S2, S4, S6, S8
  character(:)  ,allocatable                                            ::      Adv                             ! Advance atatus: 'YES' (Default) or 'NO'
  character(:)  ,allocatable                                            ::      LineFormat                      ! Format for writing an entire line
  integer                                                               ::      Local_Status, ios
  integer                                                               ::      i                               ! Index of the elements
  if (Present(Unused)) then; end if                                                                             ! Just for the "Unused" variable to be used (yes, its tricky) in order to avoid the "This variable has not been used." warning message in the compilation phase.
  if ( .Not.This%Activated ) return
  Local_Status  =       0
  call Convert_Variable_To_String( V1 , S1 , F1 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V2 , S2 , F2 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V3 , S3 , F3 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V4 , S4 , F4 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V5 , S5 , F5 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V6 , S6 , F6 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V7 , S7 , F7 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios
  call Convert_Variable_To_String( V8 , S8 , F8 , Fc, Fi, Fr, Fmt, ios ); if ( ios /= 0 ) Local_Status = ios

  LineFormat    =       This%GetPrefix(LogLevel,Error,Warning,Info,Debug,HeavyDebug) // "*(a,a,:," // Default_Spacing2_Format // "))"
  call This%Write_NewLine( NewLine )                                                                            ! Writing a new line if required
  call This%Set_Advancing( Adv, Advance )                                                                       ! Setting the line advancement status
  do i = 1,size(V2)
    write( This%Unit, LineFormat, Advance=Adv, IOStat=ios ) S1, S2(i), S3, S4(i), S5, S6(i), S7, S8(i)
  end do
  call This%Set_Backspace( Backspace )                                                                          ! Setting the backspacing if required
  if ( ios /= 0 ) Local_Status = ios                                                                            ! Setting the local status indicator
  if ( present(Status) ) Status = Local_Status
End Procedure

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

Purity Function Set_Generic_Scalar_Format( DefFormat, VarFmt, TypFormat, ComFmt ) result(Format)

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

Purity Function Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt ) result(Format)

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
Purity Elemental Function Is_Valid_Format( String ) result(Valid_Format)
  character(*)                                          ,intent(in)     ::      String
  logical                                                               ::      Valid_Format
  character(:)  ,allocatable                                            ::      String_Loc
  String_Loc    =   String
  Valid_Format  =   .True.
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

Purity Function Is_Valid( Value, Valid_Values ) result(Valid)
  implicit none
  character(*)                                  ,intent(in)     ::      Value                                   !< Value to be checked for validity
  character(*)          ,dimension( : )         ,intent(in)     ::      Valid_Values                            !< Valid values used for validity check
  logical                                                       ::      Valid                                   !< Indicator of input object validity
  integer                                                       ::      i                                       ! Index of valid strings
  Valid         =       .False.                                                                                 ! Initialization of the object validity indicator to false
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
! #ifdef WORKAROUND_GFORTRAN_SOURCE_ALLOCATION
!   integer                                                               ::      i
! #endif
!   if ( .not. allocated(Array) ) allocate( character(0) :: Array(0) )
!   Length        =       max( len(Array), len(Element) )
!   allocate( character(Length) :: Array_tmp(size(Array)+1) )
! #ifdef WORKAROUND_GFORTRAN_SOURCE_ALLOCATION
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

Purity Subroutine AddElementToArray_C0( Element, Array )
  character(*)                                          ,intent(in)     ::      Element
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      Array
  character(:)  ,dimension(:)   ,allocatable                            ::      List_Elements
  if ( .not. allocated(Array) ) allocate( character(0) :: Array(0) )
  allocate( List_Elements, source = [Array,Element] )
  call move_alloc( List_Elements, Array )
End Subroutine

Purity Subroutine AddElementToArray_C1( Elements, Array )
  character(*)  ,dimension(:)                           ,intent(in)     ::      Elements
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      Array
  character(:)  ,dimension(:)   ,allocatable                            ::      List_Elements
  if ( .not. allocated(Array) ) allocate( character(0) :: Array(0) )
  allocate( List_Elements, source = [Array,Elements] )
  call move_alloc( List_Elements, Array )
End Subroutine

Purity Subroutine Remove_Element_From_Array( Array )
  implicit none
  character(:)  ,dimension(:)   ,allocatable            ,intent(inout)  ::      Array
  integer                                                               ::      NElements
  character(:)  ,dimension(:)   ,allocatable                            ::      Array_tmp
  if ( allocated(Array) ) then
#ifdef WORKAROUND_GFORTRAN_SOURCE_ALLOCATION
    allocate( Array_tmp(NElements-1), source=Array(1:NElements-1) )
#else
    allocate( Array_tmp, source=Array(1:NElements-1) )
#endif
    call move_alloc( Array_tmp, Array )
  end if
End Subroutine

!Purity Function VecTrim( Input_String ) result(Output_String)
!  character(*)  ,dimension(:)                           ,intent(in)     ::      Input_String
!  character(:)  ,dimension(:)           ,allocatable                    ::      Output_String
!  integer                                                               ::      i
!  allocate( character(Max_Len_Trim(Input_String)) :: Output_String(size(Input_String)) )
!  do i = 1,size(Input_String)
!    Output_String(i)    =       trim( Input_String(i) )
!  end do
!End Function

Purity Function Max_Len_Trim( Strings ) result( Length )
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

Purity Subroutine Convert_Var0d_To_Str0d( Variable, String, VarFmt, ChaFmt, IntFmt, ReaFmt, ComFmt, Status )
  use iso_fortran_env ,only:  INT8, INT16, INT32, INT64, REAL32, REAL64, REAL128
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
    type is (integer(INT8));  call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(INT16)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(INT32)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(INT64)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(REAL32));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(REAL64));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(REAL128));  call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (character(*));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ChaFmt, ComFmt=ComFmt, Status=Local_Status )
    class default;       ! Error
      String        =       "?"
      Local_Status  =       -1
  end select
  if ( present(Status) ) Status = Local_Status
End Subroutine

Purity Subroutine Convert_Var1d_To_Str0d( Variable, String, VarFmt, ChaFmt, IntFmt, ReaFmt, ComFmt, Status )
  use iso_fortran_env ,only:  INT8, INT16, INT32, INT64, REAL32, REAL64, REAL128
  class(*)      ,dimension(:)                           ,intent(in)     ::      Variable                        !< Variable to be converted into a character string
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
    type is (integer(INT8));  call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(INT16)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(INT32)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(INT64)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(REAL32));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(REAL64));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(REAL128));  call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (character(*));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ChaFmt, ComFmt=ComFmt, Status=Local_Status )
    class default;       ! Error
      String        =       "?"
      Local_Status  =       -1
  end select
  if ( present(Status) ) Status = Local_Status
End Subroutine

Purity Subroutine Convert_Var1d_To_Str1d( Variable, String, VarFmt, ChaFmt, IntFmt, ReaFmt, ComFmt, Status )
  use iso_fortran_env ,only:  INT8, INT16, INT32, INT64, REAL32, REAL64, REAL128
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
    type is (logical);        call Convert_To_String( Variable, String, VarFmt=VarFmt,                   ComFmt=ComFmt, Status=Local_Status )
    type is (integer(INT8));  call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(INT16)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(INT32)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (integer(INT64)); call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(REAL32));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(REAL64));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (real(REAL128));  call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
    type is (character(*));   call Convert_To_String( Variable, String, VarFmt=VarFmt, TypFormat=ChaFmt, ComFmt=ComFmt, Status=Local_Status )
    class default;       ! Error
      String        =       "?"
      Local_Status  =       -1
  end select
  if ( present(Status) ) Status = Local_Status
End Subroutine


Purity Subroutine Convert_Var2d_To_Str1d( Variable, String, VarFmt, ChaFmt, IntFmt, ReaFmt, ComFmt, Status )
  use iso_fortran_env ,only:  INT8, INT16, INT32, INT64, REAL32, REAL64, REAL128
  class(*)      ,dimension(:,:)                         ,intent(in)     ::      Variable                        !< Variable to be converted into a character string
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String                          !< Output character string corresponding to the input variable
  character(*)                                ,optional ,intent(in)     ::      VarFmt                          !< Variable-specific format: Used if present
  character(*)                                ,optional ,intent(in)     ::      ChaFmt                          !< Character format: Used if variable is of type 'character(*)'and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      IntFmt                          !< Integer format: Used if variable is of type 'integer(*)'and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      ReaFmt                          !< Reals format: Used if variable is of type 'real(*)' and if no variable-specific format
  character(*)                                ,optional ,intent(in)     ::      ComFmt                          !< Common format: Used if the variable-specific format 'Fv'  and the type-specific format 'Fc/Fi/Fr' associated to current variable's type are absent
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)

  integer                                                               ::      Local_Status
  integer                                                               ::      i, Length, NElements
  character(:)  ,allocatable                                            ::      String_0d
  character(:)  ,dimension(:) ,allocatable                              ::      StringTmp
  Length    =   10000
  NElements =   size(Variable,1)
  allocate( character(Length) :: String(NElements) )  ! Allocating to the number of rows
!   allocate( character(Length) :: StringTmp(NElements) )  ! Allocating to the number of rows

  select type (Variable)
    type is (logical)
      do i = 1,NElements
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt,                   ComFmt=ComFmt, Status=Local_Status )
        Length    =   len_trim(String_0d)
        if ( Length <= len(String) ) then
          String(i) = String_0d
        else
          allocate( Character(Length) :: StringTmp(NElements) )
          if (i>1) StringTmp(1:i-1) = String(1:i-1)
          StringTmp(i)  =   String_0d
          call move_alloc( StringTmp , String )
        end if
      end do

    type is (integer(INT8))
      do i = 1,NElements
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
        Length    =   len_trim(String_0d)
        if ( Length <= len(String) ) then
          String(i) = String_0d
        else
          allocate( Character(Length) :: StringTmp(NElements) )
          if (i>1) StringTmp(1:i-1) = String(1:i-1)
          StringTmp(i)  =   String_0d
          call move_alloc( StringTmp , String )
        end if
      end do

    type is (integer(INT16))
      do i = 1,NElements
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
        Length    =   len_trim(String_0d)
        if ( Length <= len(String) ) then
          String(i) = String_0d
        else
          allocate( Character(Length) :: StringTmp(NElements) )
          if (i>1) StringTmp(1:i-1) = String(1:i-1)
          StringTmp(i)  =   String_0d
          call move_alloc( StringTmp , String )
        end if
      end do

    type is (integer(INT32))
      do i = 1,NElements
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
        Length    =   len_trim(String_0d)
        if ( Length <= len(String) ) then
          String(i) = String_0d
        else
          allocate( Character(Length) :: StringTmp(NElements) )
          if (i>1) StringTmp(1:i-1) = String(1:i-1)
          StringTmp(i)  =   String_0d
          call move_alloc( StringTmp , String )
        end if
      end do

    type is (integer(INT64))
      do i = 1,NElements
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=IntFmt, ComFmt=ComFmt, Status=Local_Status )
        Length    =   len_trim(String_0d)
        if ( Length <= len(String) ) then
          String(i) = String_0d
        else
          allocate( Character(Length) :: StringTmp(NElements) )
          if (i>1) StringTmp(1:i-1) = String(1:i-1)
          StringTmp(i)  =   String_0d
          call move_alloc( StringTmp , String )
        end if
      end do

    type is (real(REAL32))
      do i = 1,NElements
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
        Length    =   len_trim(String_0d)
        if ( Length <= len(String) ) then
          String(i) = String_0d
        else
          allocate( Character(Length) :: StringTmp(NElements) )
          if (i>1) StringTmp(1:i-1) = String(1:i-1)
          StringTmp(i)  =   String_0d
          call move_alloc( StringTmp , String )
        end if
      end do

    type is (real(REAL64))
      do i = 1,NElements
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
        Length    =   len_trim(String_0d)
        if ( Length <= len(String) ) then
          String(i) = String_0d
        else
          allocate( Character(Length) :: StringTmp(NElements) )
          if (i>1) StringTmp(1:i-1) = String(1:i-1)
          StringTmp(i)  =   String_0d
          call move_alloc( StringTmp , String )
        end if
      end do

    type is (real(REAL128))
      do i = 1,NElements
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=ReaFmt, ComFmt=ComFmt, Status=Local_Status )
        Length    =   len_trim(String_0d)
        if ( Length <= len(String) ) then
          String(i) = String_0d
        else
          allocate( Character(Length) :: StringTmp(NElements) )
          if (i>1) StringTmp(1:i-1) = String(1:i-1)
          StringTmp(i)  =   String_0d
          call move_alloc( StringTmp , String )
        end if
      end do

    type is (character(*))
      do i = 1,NElements
        call Convert_To_String( Variable(i,:), String_0d, VarFmt=VarFmt, TypFormat=ChaFmt, ComFmt=ComFmt, Status=Local_Status )
        Length    =   len_trim(String_0d)
        if ( Length <= len(String) ) then
          String(i) = String_0d
        else
          allocate( Character(Length) :: StringTmp(NElements) )
          if (i>1) StringTmp(1:i-1) = String(1:i-1)
          StringTmp(i)  =   String_0d
          call move_alloc( StringTmp , String )
        end if
      end do

    class default;       ! Error
      deallocate( String )
      allocate( String, source = ["?"] )
      Local_Status  =       -1

  end select
   !String    =   VecTrim(String)

  if ( present(Status) ) Status = Local_Status
End Subroutine


Purity Subroutine Convert_Logical_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
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

Purity Subroutine Convert_INT8_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  INT8
  integer(INT8)                                         ,intent(in)     ::      Variable
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

Purity Subroutine Convert_INT16_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  INT16
  integer(INT16)                                        ,intent(in)     ::      Variable
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

Purity Subroutine Convert_INT32_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  INT32
  integer(INT32)                                        ,intent(in)     ::      Variable
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

Purity Subroutine Convert_INT64_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  INT64
  integer(INT64)                                        ,intent(in)     ::      Variable
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

Purity Subroutine Convert_REAL32_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  REAL32
  real(REAL32)                                          ,intent(in)     ::      Variable                             !< Real number to be converted into a string
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

Purity Subroutine Convert_REAL64_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  REAL64
  real(REAL64)                                          ,intent(in)     ::      Variable                             !< Real number to be converted into a string
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

Purity Subroutine Convert_REAL128_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  REAL128
  real(REAL128)                                         ,intent(in)     ::      Variable                             !< Real number to be converted into a string
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

Purity Subroutine Convert_String_To_String_0D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
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

Purity Subroutine Convert_Logical_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  logical       ,dimension(:)                           ,intent(in)     ::      Variable
  character(*)  ,parameter                                              ::      DefFormat = Default_Logical_Format
#include "Logger_NumberToString_Inline.F90"
End Subroutine

Purity Subroutine Convert_INT8_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  INT8
  integer(INT8) ,dimension(:)                           ,intent(in)     ::      Variable
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
#include "Logger_NumberToString_Inline.F90"
End Subroutine

Purity Subroutine Convert_INT16_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  INT16
  integer(INT16) ,dimension(:)                          ,intent(in)     ::      Variable
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
#include "Logger_NumberToString_Inline.F90"
End Subroutine

Purity Subroutine Convert_INT32_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  INT32
  integer(INT32) ,dimension(:)                          ,intent(in)     ::      Variable
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
#include "Logger_NumberToString_Inline.F90"
End Subroutine

Purity Subroutine Convert_INT64_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  INT64
  integer(INT64),dimension(:)                           ,intent(in)     ::      Variable
  character(*)  ,parameter                                              ::      DefFormat = Default_Integer_Format
#include "Logger_NumberToString_Inline.F90"
End Subroutine

Purity Subroutine Convert_REAL32_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  REAL32
  real(REAL32)  ,dimension(:)                           ,intent(in)     ::      Variable
  character(*)  ,parameter                                              ::      DefFormat = Default_Real_Format
#include "Logger_NumberToString_Inline.F90"
End Subroutine

Purity Subroutine Convert_REAL64_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  REAL64
  real(REAL64)  ,dimension(:)                           ,intent(in)     ::      Variable
  character(*)  ,parameter                                              ::      DefFormat = Default_Real_Format
#include "Logger_NumberToString_Inline.F90"
End Subroutine

Purity Subroutine Convert_REAL128_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  REAL128
  real(REAL128) ,dimension(:)                           ,intent(in)     ::      Variable
  character(*)  ,parameter                                              ::      DefFormat = Default_Real_Format
#include "Logger_NumberToString_Inline.F90"
#ifdef NODEF
  character(10000)                                                      ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios
  Local_Format    =       Set_Generic_Vector_Format( DefFormat, VarFmt, TypFormat, ComFmt )
  Local_Format    =       "(" // Local_Format // ")"
  write( Long_String , Local_Format , iostat=ios ) Variable
  if ( ios /= 0 ) write(Long_String ,"(g0)") Variable
  String          =       trim(Long_String)
  if ( present(Status) ) Status = ios
#endif
End Subroutine

Purity Subroutine Convert_String_To_String_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
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
    do i = 1,size(Variable)
      String        =       String // Variable(i) // "   "
    end do
    String          =       trim(String)
    ios             =       0
  end if
  if ( present(Status) ) Status = ios
End Subroutine

Purity Subroutine Convert_Logical_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  logical       ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, i
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
  !allocate( String , source = VecTrim(Long_String), Stat=ios )
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Purity Subroutine Convert_INT8_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  INT8
  integer(INT8) ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat
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
  !allocate( String , source = VecTrim(Long_String), Stat=ios )
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Purity Subroutine Convert_INT16_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
 use iso_fortran_env ,only:  INT16
  integer(INT16),dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat
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
  !allocate( String , source = VecTrim(Long_String), Stat=ios )
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Purity Subroutine Convert_INT32_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  INT32
  integer(INT32),dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat
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
  !allocate( String , source = VecTrim(Long_String), Stat=ios )
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Purity Subroutine Convert_INT64_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  INT64
  integer(INT64),dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat
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
  !allocate( String , source = VecTrim(Long_String), Stat=ios )
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Purity Subroutine Convert_REAL32_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  REAL32
  real(REAL32)  ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, i
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
  !allocate( String , source = VecTrim(Long_String), Stat=ios )
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Purity Subroutine Convert_REAL64_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  REAL64
  real(REAL64)  ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, i
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
  !allocate( String , source = VecTrim(Long_String), Stat=ios )
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Purity Subroutine Convert_REAL128_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  use iso_fortran_env ,only:  REAL128
  real(REAL128) ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, i
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
  !allocate( String , source = VecTrim(Long_String), Stat=ios )
  if ( ios /= 0 ) Stat = ios
  if ( present(Status) ) Status = Stat
End Subroutine

Purity Subroutine Convert_String_To_Strings_1D( Variable, String, VarFmt, TypFormat, ComFmt, Status )
  character(*)  ,dimension(:)                           ,intent(in)     ::      Variable
  character(:)  ,dimension(:) ,allocatable              ,intent(out)    ::      String
  character(*)                                ,optional ,intent(in)     ::      VarFmt
  character(*)                                ,optional ,intent(in)     ::      TypFormat
  character(*)                                ,optional ,intent(in)     ::      ComFmt
  integer                                     ,optional ,intent(out)    ::      Status                          !< Error status indicator (=0 if everthing is ok, /=0 if error)
  character(10000)  ,dimension( size(Variable) )                        ::      Long_String                         ! Local character string required to store the number of a very large string before allocation of the output variable
  character(:)  ,allocatable                                            ::      Local_Format
  integer                                                               ::      ios, Stat, i
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
    !allocate( String , source = VecTrim(Long_String), Stat=ios )
    if ( ios /= 0 ) Stat = ios
  else
    !allocate( String , source = VecTrim(Variable), Stat=Stat )
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

Purity Function LowerCase(StrInp) result(StrOut)
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
!     type is (real(REAL32));     Format = Set_Generic_Scalar_Format( Default_Real_Format,       Optional_Format )
!     type is (real(REAL64));     Format = Set_Generic_Scalar_Format( Default_Real_Format,       Optional_Format )
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
!     type is (real(REAL32));     Format = Set_Generic_Scalar_Format( Default_Real_Format,       Optional_Format )
!     type is (real(REAL64));     Format = Set_Generic_Scalar_Format( Default_Real_Format,       Optional_Format )
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
!     type is (real(REAL32));     Format = Set_Generic_Vector_Format( Default_Real_Format,       Optional_Format )
!     type is (real(REAL64));     Format = Set_Generic_Vector_Format( Default_Real_Format,       Optional_Format )
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
!     type is (real(REAL32));     Format = Set_Generic_Vector_Format( Default_Real_Format,       Optional_Format )
!     type is (real(REAL64));     Format = Set_Generic_Vector_Format( Default_Real_Format,       Optional_Format )
!     type is (character(*));     Format = Set_Generic_Vector_Format( Default_Character_Format,  Optional_Format )
!     class default       ! Error
!   end select
! End Function

Purity Function Set_Optional_Argument_Logical( VarDef, VarOpt ) result(VarLoc)
  logical                                                       ,intent(in)     ::      VarDef                  !< Default value of the local variable if the optional variable is absent
  logical       ,optional                                       ,intent(in)     ::      VarOpt                  !< Optional varibale
  logical                                                                       ::      VarLoc                  !< Local variable to be set
  if ( present(VarOpt) ) then;  VarLoc = VarOpt                                                                 ! Setting the local variable to the optional variable if present ...
  else;                         VarLoc = VarDef; end if                                                         ! ... otherwise setting the local variable to the default variable
End Function

Purity Function Set_Optional_Argument_Integer( VarDef, VarOpt ) result(VarLoc)
  integer                                                       ,intent(in)     ::      VarDef                  !< Default value of the local variable if the optional variable is absent
  integer       ,optional                                       ,intent(in)     ::      VarOpt                  !< Optional varibale
  integer                                                                       ::      VarLoc                  !< Local variable to be set
  if ( present(VarOpt) ) then;  VarLoc = VarOpt                                                                 ! Setting the local variable to the optional variable if present ...
  else;                         VarLoc = VarDef; end if                                                         ! ... otherwise setting the local variable to the default variable
End Function

Purity Function Set_Optional_Argument_Real( VarDef, VarOpt ) result(VarLoc)
  use iso_fortran_env ,only:  REAL64
  real(REAL64)                                                  ,intent(in)     ::      VarDef                  !< Default value of the local variable if the optional variable is absent
  real(REAL64)  ,optional                                       ,intent(in)     ::      VarOpt                  !< Optional varibale
  real(REAL64)                                                                  ::      VarLoc                  !< Local variable to be set
  if ( present(VarOpt) ) then;  VarLoc = VarOpt                                                                 ! Setting the local variable to the optional variable if present ...
  else;                         VarLoc = VarDef; end if                                                         ! ... otherwise setting the local variable to the default variable
End Function

Purity Function Set_Optional_Argument_Character( VarDef, VarOpt ) result(VarLoc)
  character(*)                                                  ,intent(in)     ::      VarDef                  !< Default value of the local variable if the optional variable is absent
  character(*)  ,optional                                       ,intent(in)     ::      VarOpt                  !< Optional varibale
  character(:)  ,allocatable                                                    ::      VarLoc                  !< Local variable to be set
  if ( present(VarOpt) ) then;  VarLoc = VarOpt                                                                 ! Setting the local variable to the optional variable if present ...
  else;                         VarLoc = VarDef; end if                                                         ! ... otherwise setting the local variable to the default variable
End Function


Purity Function LogLevelToString( LogLevel ) result(String)
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

Purity Function IsPresent( Element, Array ) result(PresenceIndicator)
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

Purity Function PresentAndTrue( OptionalArgument ) result(Indicator)
  logical                                     ,optional ,intent(in)     ::      OptionalArgument
  logical                                                               ::      Indicator
  Indicator   =   .False.
  if ( Present(OptionalArgument) ) Indicator  = OptionalArgument
End Function

Subroutine SetInitialLoggerItem( This, Procedure, Indentation )
  type(Logger_Type)                                     ,intent(inout)  ::      This                        !< Passed-object dummy argument
  character(*)                                ,optional ,intent(in)     ::      Procedure            !< Names of the calling procedures to reach the current procedure. Each procedure neames are separated by the '>' character.  Only used by the Logger object
  integer                                     ,optional ,intent(in)     ::      Indentation                 !< Indentation level
  character(:)  ,allocatable                                            ::      ProcedureName
  character(:)  ,allocatable  ,dimension(:)                             ::      ProcedureNames                  ! Vector containing the name of the calling procedures
  type(LoggerItem_Type)                                                 ::      LoggerItem
  integer                                                               ::      i                               ! Index of procedure names
  integer                                                               ::      ItemIndent
  integer                                                               ::      Indentation_
  call GetProcedureNames( Procedure, ProcedureNames )
  Indentation_ =    2
  if ( present(Indentation) ) Indentation_ = Indentation
  ItemIndent   =    0
  do i = 1,size(ProcedureNames)
    ProcedureName   =   trim( ProcedureNames(i) )                                                                 ! Setting the procedure name used to construct the current Logger-Item object
    ItemIndent      =   ItemIndent + Indentation_                                                                 ! Setting the indentation level used to construct the current Logger-Item object
    call LoggerItem%Initialize( ProcedureName, ItemIndent )                                                       ! Initializing a 'LoggerItem' object
    call This%AddItem( LoggerItem )                                                                               ! Adding the new 'LoggerItem' object to the list of 'LoggerItem' sub-objects
  end do
!   call LoggerItem%Initialize( Procedure, Indentation )                                                          ! Initializing a 'LoggerItem' object with the input properties for the name of the procedure and the indentation level
!   call This%AddItem( LoggerItem )                                                                               ! Adding the new 'LoggerItem' object to the list of 'LoggerItem' sub-objects
End Subroutine


Subroutine GetProcedureNames( Procedures, Names )
  character(*)                                ,optional ,intent(in)     ::      Procedures                      !< Names of the calling procedures to reach the current procedure. Each procedure neames are separated by the '>' character.  Only used by the Logger object
  character(:)  ,allocatable  ,dimension(:)             ,intent(out)    ::      Names                 !< Vector containing the name of the calling procedures

  character(*)  ,parameter                                              ::      Separator = ">"                 ! Character separating the names of procedures
  integer                                                               ::      iSep, jSep                      ! Position iof the separator character
  integer                                                               ::      NProc                           ! Number of procedure names
  integer                                                               ::      TotLength                       ! Length of the input character string
  integer                                                               ::      MaxLength                       ! Maximum length of the procedure names
  integer                                                               ::      i                               ! Index of procedure names
  character(:)  ,allocatable                                            ::      CurrentNames
  character(:)  ,allocatable                                            ::      LeftNames

  TotLength =   len_trim(Procedures)

  if ( .Not. present(Procedures) ) then
    allocate( character(0) :: Names(0) )
    return
  end if

  if ( TotLength == 0 ) then
    allocate( character(0) :: Names(0) )
    return
  end if

  LeftNames     =   Procedures                                                                                ! Copying the input character string
  NProc         =   0                                                                                         ! Initializing the number of procedures
  MaxLength     =   0                                                                                         ! Initializing the length of procedures
  do                                                                                                          ! Loop
    if ( len_trim(LeftNames) == 0 ) exit                                                                      ! If the string which still to be processed, then exiting the loop
    iSep        =       index(LeftNames,Separator)                                                            ! Getting the position of the separator character
    if ( iSep == 0  ) then                                                                                    ! If the separator character has not been found, then there is only one procedure
      CurrentNames  =   LeftNames                                                                             ! Setting the name of current procedure to the entire string
      LeftNames     =   ''                                                                                    ! Setting the string which still to be processed to an empty string
    else                                                                                                      ! If the separator character has not been found, then
      CurrentNames  =   trim( LeftNames(1:iSep-1) )                                                           ! Setting the name of current procedure to the string at the LHS of the separator
      if ( iSep == TotLength ) then; LeftNames   =   ''                                                       ! If the separator is the last character, then setting the string which still to be processed to an empty string
      else; LeftNames = LeftNames(iSep+1:); end if                                                            ! ... otherwise setting the string which still to be processed to the string at the RHS of the separator
    end if                                                                                                    ! End if case on the presence of the separator
    NProc       =   NProc + 1                                                                                 ! Incrementing the number of procedures
    MaxLength   =   max( MaxLength , len_trim(CurrentNames) )                                                 ! Updating the maximum length of the procedure names
  end do

  allocate( character(MaxLength) :: Names(NProc) )
  LeftNames     =   Procedures                                                                                ! Copying the input character string
  i             =   0                                                                                         ! Initializing the index of procedures name
  do                                                                                                          ! Loop
    if ( len_trim(LeftNames) == 0 ) exit                                                                      ! If the string which still to be processed, then exiting the loop
    iSep        =       index(LeftNames,Separator)                                                            ! Getting the position of the separator character
    if ( iSep == 0  ) then                                                                                    ! If the separator character has not been found, then there is only one procedure
      CurrentNames  =   LeftNames                                                                             ! Setting the name of current procedure to the entire string
      LeftNames     =   ''                                                                                    ! Setting the string which still to be processed to an empty string
    else                                                                                                      ! If the separator character has not been found, then
      CurrentNames  =   trim( LeftNames(1:iSep-1) )                                                           ! Setting the name of current procedure to the string at the LHS of the separator
      if ( iSep == TotLength ) then; LeftNames   =   ''                                                       ! If the separator is the last character, then setting the string which still to be processed to an empty string
      else; LeftNames = LeftNames(iSep+1:); end if                                                            ! ... otherwise setting the string which still to be processed to the string at the RHS of the separator
    end if                                                                                                    ! End if case on the presence of the separator
    i               =   i + 1                                                                                 ! Incrementing the names index
    Names(i)        =   CurrentNames                                                                          ! Setting the current procedure name
  end do
End Subroutine


End SubModule
