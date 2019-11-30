function(add_sources)
  get_property(is_defined GLOBAL PROPERTY SRCS_LIST DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY SRCS_LIST
      BRIEF_DOCS "List of source files"
      FULL_DOCS  "List of source files to be compiled in one library")
  endif()
  # make absolute paths
  set(SRCS)
  foreach(s IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${s}")
      get_filename_component(s "${s}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${s}")
  endforeach()
  # append to global list
  set_property(GLOBAL APPEND PROPERTY SRCS_LIST "${SRCS}")
endfunction(add_sources)

# This function gets the list of source files stored in the 'SRCS_LIST' variable
# and reset the variable to an empty string.
function( get_sources SRCS )
  get_property( SRCS_LOCAL GLOBAL PROPERTY SRCS_LIST )
  SET( SRCS ${SRCS_LOCAL} PARENT_SCOPE )
  set_property( GLOBAL PROPERTY SRCS_LIST "" )                                        # Clear the list of source files from the global variable
endfunction(get_sources)

