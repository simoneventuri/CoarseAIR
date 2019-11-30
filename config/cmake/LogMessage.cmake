#==========================================================================
# Utilities for writing cmake log messages.
#==========================================================================
# File Name:      LogMessage.cmake
# Functions:
#   LogMessage    Write a log message
#   setLogsOn     Activate logs
#   setLogsOff    Deactivate logs
#==========================================================================

function(LogMessage message)
  if (LOG_MESSAGE_STATUS)
    set ( LOG_MESSAGE_PREFIX "[${PROJECT_NAME}]:" )
    message ( STATUS "${LOG_MESSAGE_PREFIX} ${message}" )
  endif()
endfunction(LogMessage)

function( setLogsOn )
  set( LOG_MESSAGE_STATUS "ON" CACHE INTERNAL "Status for log messages" )
endfunction( setLogsOn )

function( setLogsOff )
  set( LOG_MESSAGE_STATUS "OFF" CACHE INTERNAL "Status for log messages" )
endfunction( setLogsOff )

# function( setLogLevel TARGET_LOG_LEVEL )
#   set( LOG_LEVEL "INFO" CACHE STRING "Select the logs level")
#   set_property( CACHE LOG_LEVEL PROPERTY STRINGS "ERROR" "WARNING" "INFO" "DEBUG" "HEAVYDEBUG" )
# #   TARGET_LOG_LEVEL
# endfunction( setLogsOn )
#
#
# function(InfoMessage message)
#   if ( LOG_LEVEL STREQUAL "Coarray" )
#     set ( LOG_MESSAGE_PREFIX "[${PROJECT_NAME}]:" )
#     message ( STATUS "${LOG_MESSAGE_PREFIX} ${message}" )
#   endif()
# endfunction(LogMessage)


