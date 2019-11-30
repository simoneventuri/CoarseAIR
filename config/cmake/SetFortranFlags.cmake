

#==============================================================================
# LOADING OPENBLAS
#==============================================================================
  message ( STATUS "[OPENBLAS]: Configuring OpenBlas" )
  if ( NOT DEFINED USE_OPENBLAS )
    set( USE_OPENBLAS "YES" CACHE STRING "Using the OpenBlas library" FORCE )
  endif ()
  set_property ( CACHE USE_OPENBLAS PROPERTY STRINGS "NO" "YES" )
  message ( STATUS "[OPENBLAS]: USE_OPENBLAS = ${USE_OPENBLAS}" )
  if ( USE_OPENBLAS STREQUAL "YES" )
    message ( STATUS "[OPENBLAS]: Including 'FindOpenBLAS' module" )
    include(FindOpenBLAS)
    find_package( OpenBLAS REQUIRED )
    set( OpenBLAS_LIB_DIR "${OpenBLAS_INCLUDE_DIR}/../lib" )
    message ( STATUS "[OPENBLAS]: OpenBLAS             =${OpenBLAS}" )
    message ( STATUS "[OPENBLAS]: OpenBLAS_LIB         =${OpenBLAS_LIB}" )
    message ( STATUS "[OPENBLAS]: OpenBLAS_INCLUDE_DIR =${OpenBLAS_INCLUDE_DIR}" )
    message ( STATUS "[OPENBLAS]: OpenBLAS_LIB_DIR     =${OpenBLAS_LIB_DIR}" )
    set( OPENBLAS_FLAGS "-L${OpenBLAS_LIB_DIR} -lopenblas" )
    add_definitions( -DUSE_OPENBLAS )
  elseif ( USE_OPENBLAS STREQUAL "NO" )
  else()
    message(FATAL_ERROR "Wrong value for 'USE_OPENBLAS': ${USE_OPENBLAS}")
  endif()
#==============================================================================

#==============================================================================
# MKL
#==============================================================================
    message ( STATUS "[MKL]: Configuring MKL" )
    if ( NOT DEFINED USE_MKL )
      set( USE_MKL "YES" CACHE STRING "Using the intel MKL library" FORCE )
    endif ()
    set_property ( CACHE USE_MKL PROPERTY STRINGS "NO" "YES" )
#     if ( (USE_OPENBLAS STREQUAL "NO") AND (USE_LAPACK STREQUAL "NO") )
#       message ( STATUS "[MKL]: Neither OpenBlas nor Lapack libraries are used => Using MKL by default" )
#       set( USE_MKL "YES" )
#     else()
#       set( USE_MKL "NO" )
#     endif()
    message ( STATUS "[MKL]: USE_MKL = ${USE_MKL}" )
    if ( USE_MKL STREQUAL "YES" )
      set( MKL_FLAGS "-mkl -DUSE_MKL" )  # -DMKL_INLINE
      add_definitions( -DUSE_MKL )
#       set( USE_LAPACK "YES" )
    elseif ( USE_MKL STREQUAL "NO" )
    else()
      message(FATAL_ERROR "Wrong value for 'USE_MKL': ${USE_MKL}")
    endif()
#==============================================================================


#==============================================================================
# LAPACK
#==============================================================================
    message ( STATUS "[LAPACK]: Configuring Lapack" )
#    set( USE_LAPACK "YES" )
    if ( NOT DEFINED USE_LAPACK )
      set( USE_LAPACK "YES" CACHE STRING "Using the Lapack library" FORCE )
    endif ()
    message ( STATUS "[LAPACK]: USE_LAPACK = ${USE_LAPACK}" )
    set_property ( CACHE USE_LAPACK PROPERTY STRINGS "NO" "YES" )
    if ( USE_LAPACK STREQUAL "YES" )
      find_package( LAPACK REQUIRED )
      message ( STATUS "[LAPACK]: LAPACK_FOUND = " ${LAPACK_FOUND} )
endif()
#==============================================================================

#==============================================================================
# LOADING OPENMP
#==============================================================================
set( OMP_FLAGS "-DUSE_OPENMP -fopenmp")

#==============================================================================


  set( OMP_FLAGS "")
#   set( OPENBLAS_FLAGS "" )
#   set( MKL_FLAGS "" )

#==============================================================================
# Setup the Fortran compiler flags depending on the compiler used
#==============================================================================
# First, the name of the compiler is found. Then, this name is used to set the
# compiler flags
#------------------------------------------------------------------------------
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message("Setting the gfortran compiler flags... ")

  add_definitions( -DGFORTRAN_WORKAROUND_SOURCE_ALLOCATION )
  add_definitions( -DGFORTRAN_WORKAROUND_PUREPROC_POLYMORPHIC )    #    INTENT(OUT) argument 'lhs' of pure procedure ‘...’ at (1) may not be polymorphic
  add_definitions( -DGCC_COMPILER )
  add_definitions( -DWORKAROUND_GFORTRAN_SOURCE_ALLOCATION )        # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=44672
  add_definitions( -DWORKAROUND_GFORTRAN_RECURSIVE_ASSOCIATION )    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=64678
  add_definitions( -DWORKAROUND_GFORTRAN_ASSIGN_ALLOCATABLE_CHARACTER )
  add_definitions( -DWORKAROUND_GFORTRAN_INQUIRE_DIRECTORY )
  add_definitions( -DWORKAROUND_GFORTRAN_SELECT_TYPE )


  if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
    message( STATUS "Setting the gfortran Fortran compiler flags for Windows... ")
    message(FATAL_ERROR "... Windows builds are currently not supported")
  else ()
    set( COMMON_FLAG "-cpp -fno-unsafe-math-optimizations -fdefault-double-8 -fdefault-real-8 -ffree-line-length-none -frealloc-lhs")
#      -fno-realloc-lhs
    set( COMMON_FLAG "${COMMON_FLAG} ${OPENBLAS_FLAGS} ")
    if ( USE_LAPACK STREQUAL "YES" )
      set( COMMON_FLAG "${COMMON_FLAG}" )
#      set( COMMON_FLAG "${COMMON_FLAG} -llapack" )
    endif()


    if ( CMAKE_BUILD_TYPE STREQUAL "debug" )
      set( DEBUG_FLAG "-g  -fbounds-check -fbacktrace -fdump-core -ggdb -pg -Wall -Wno-line-truncation -finit-real=inf")
      set( DEBUG_FLAG "${DEBUG_FLAG} -Wno-unused-dummy-argument -Wno-unused-function -Wno-uninitialized -Wno-maybe-uninitialized -Wno-surprising -Wno-unused-variable" )
      set( OPTIM_FLAG "-O0")
#       -Wrealloc-lhs-all -Wrealloc-lhs
      set( FFLAGS "${COMMON_FLAG} ${OPTIM_FLAG} ${DEBUG_FLAG} " )
    endif ()
    
    if ( CMAKE_BUILD_TYPE STREQUAL "release" )
      #set( DEBUG_FLAG "-g  -fbounds-check -fbacktrace -fdump-core -ggdb -pg -Wall -Wno-line-truncation ")
      #set( DEBUG_FLAG "${DEBUG_FLAG} -Wno-unused-dummy-argument -Wno-unused-function -Wno-uninitialized -Wno-maybe-uninitialized -Wno-surprising -Wno-unused-variable -Wno-conversion" )
      set( OPTIM_FLAG "-O3 -fPIC -fbounds-check -fbacktrace ${OMP_FLAGS}")
      set( FFLAGS "${COMMON_FLAG} ${OPTIM_FLAG} ${DEBUG_FLAG}" )
    endif ()
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS}")
    set( CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")                                                    # Set the share library link flags: This is to get ride of the "-i_dynamic" option
  endif()

elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")

  if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
    message( STATUS "Setting the Intel Fortran compiler flags for Windows... ")
    message(FATAL_ERROR "... Windows builds are currently not supported")
  else ()
    message( STATUS "Setting the Intel Fortran compiler options for Linux/X-OS")


    set( COMMON_FLAG " -fpp -diag-disable 7712,6717,406")
    set( COMMON_FLAG "${COMMON_FLAG} ${MKL_FLAGS}  ")

    if ( CMAKE_BUILD_TYPE STREQUAL "debug" )
      set( DEBUG_FLAG "-fPIC -ftrapuv -auto -check all -traceback -debug-parameters all -g -pg -warn all -warn uninitialized -warn usage -warn uncalled -warn declarations -warn general -warn alignments -implicitnone -nozero -check noarg_temp_created")
      set( OPTIM_FLAG "-O0")
      set( FFLAGS "${COMMON_FLAG} ${OPTIM_FLAG} ${DEBUG_FLAG}" )
    endif ()
    if ( CMAKE_BUILD_TYPE STREQUAL "release" )
      set( DEBUG_FLAG "")
#      set( OPTIM_FLAG "-fast -fPIC -r8")
      set( OPTIM_FLAG "-ipo -O3 -no-prec-div -xHost -xAVX -fPIC -r8")
      set( FFLAGS "${COMMON_FLAG} ${OPTIM_FLAG} ${DEBUG_FLAG}" )
endif()

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS}")
    set( CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")                                                    # Set the share library link flags: This is to get ride of the "-i_dynamic" option


# #==============================================================================
# # LAPACK
# #==============================================================================
#     message ( STATUS "[LAPACK]: Configuring Lapack" )
#     if ( NOT DEFINED USE_LAPACK )
#       set( USE_LAPACK "NO" CACHE STRING "Using the Lapack library" FORCE )
#     endif ()
#     message ( STATUS "[LAPACK]: USE_LAPACK = ${USE_LAPACK}" )
#     set_property ( CACHE USE_LAPACK PROPERTY STRINGS "NO" "YES" )
#     if ( USE_LAPACK STREQUAL "YES" )
# #       set( FFLAGS "-llapack" )
# #       message ( STATUS "[LAPACK]: Adding Fortran flags: ${FFLAGS}" )
# #       set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS}")
#
# #       find_package( Lapack REQUIRED )
#
#       find_package( LAPACK REQUIRED )
#       message ( STATUS "[LAPACK]: LAPACK_FOUND = " ${LAPACK_FOUND} )
#
#     endif()
# #==============================================================================
#
#
  endif()

else ()
   message(FATAL_ERROR "Unknown compiler: ${Fortran_COMPILER_NAME}")
endif()

add_definitions( -DVVCT_VERSION=${PROJECT_VERSION} )
add_definitions( -DVVCT_VERSION_MAJOR=${PROJECT_VERSION_MAJOR} )
add_definitions( -DVVCT_VERSION_MINOR=${PROJECT_VERSION_MINOR} )
add_definitions( -DVVCT_VERSION_PATCH=${PROJECT_VERSION_PATCH} )

message ( STATUS "Debugging Fortran options:    DEBUG_FLAG: ${DEBUG_FLAG}")
message ( STATUS "Optimization Fortran options: OPTIM_FLAG: ${OPTIM_FLAG}")
message ( STATUS "Fortran options:     CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
#==============================================================================

