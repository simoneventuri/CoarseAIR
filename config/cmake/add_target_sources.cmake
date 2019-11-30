#- Add sources for a target
#
#  add_target_sources(<target> <source1> [<source2> ...])
#
function(add_target_sources target)
  # define the <target>_SRCS properties if necessary
  get_property(is_defined GLOBAL PROPERTY ${target}_SRCS DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY ${target}_SRCS
      BRIEF_DOCS "Sources for the ${target} target"
      FULL_DOCS "List of source files for the ${target} target")
  endif()
  # create list of sources (absolute paths)
  set(SRCS)
  foreach(s IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${s}")
      get_filename_component(s "${s}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${s}")
  endforeach()
  # append to global property
  set_property(GLOBAL APPEND PROPERTY "${target}_SRCS" "${SRCS}")
endfunction()
