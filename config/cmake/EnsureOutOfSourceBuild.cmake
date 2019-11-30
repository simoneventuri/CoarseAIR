# - ensure_out_of_source_build()
# Ensures that the build directory is different from source directory.
# If the build directory is the source directory then it will bump
# an error message and stop the compilation processus.

macro (ensure_out_of_source_build)
  string ( COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" In_Source)
#   string ( REGEX MATCH "${CMAKE_SOURCE_DIR}/" In_Source_Subdir "${CMAKE_BINARY_DIR}")

  get_filename_component ( PARENTDIR ${CMAKE_SOURCE_DIR} PATH )
  string ( COMPARE EQUAL "${CMAKE_SOURCE_DIR}"  "${PARENTDIR}" In_Source_Subdir)

  if ( In_Source OR In_Source_Subdir )
  message(FATAL_ERROR "
ERROR BUILDING ${PROJECT_NAME}:
In source builds of this project are not allowed.
A separate build directory is required.
Please create one and run cmake from the build directory.
Also note that cmake has just added files to your source code directory.
We suggest getting a new copy of the source code.
Otherwise, delete `CMakeCache.txt' and the directory `CMakeFiles'.
  ")
  endif ()

# message("
# CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}
# CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}
# PARENTDIR: ${PARENTDIR}
# PATH: ${PATH}
# ")

endmacro (ensure_out_of_source_build)
