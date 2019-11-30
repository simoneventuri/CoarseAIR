cmake_minimum_required (VERSION 3.00)
project( CoarseAIR VERSION 1.1)
set( PROJECT_ACRONYM_DEFINITION "Quasi-Classical Trajectory Code" )

######################################################
# Setting the directory which contains CMake modules #
######################################################
# Appends the cmake path inside the MAKE_MODULE_PATH variable which stores the
# list of directories where to lookk for additional CMake modules.
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config/cmake ${CMAKE_MODULE_PATH} )


###########################
# Including cmake modules #
###########################
include( LogMessage )                   # Include utilities to raise log message
include( add_sources )
include( EnsureOutOfSourceBuild )
ensure_out_of_source_build()            # Ensures an out-of-source-build.
setLogsOn()
# setLogsOff()


#################################################
# Setting variables related to the project name #
#################################################
string( TOLOWER ${PROJECT_NAME} PROJECT_NAME_LC )
string( TOUPPER ${PROJECT_NAME} PROJECT_NAME_UC )
string( TOLOWER ${PROJECT_NAME} LIBRARY_NAME )
set( ${PROJECT_NAME_UC}_LIBRARY_NAME ${LIBRARY_NAME} )
set( STATIC_LIBRARY ${LIBRARY_NAME}-static )
LogMessage( "PROJECT_NAME:          ${PROJECT_NAME}" )
LogMessage( "PROJECT_VERSION:       ${PROJECT_VERSION}" )
LogMessage( "LIBRARY_NAME:          ${LIBRARY_NAME}" )
LogMessage( "STATIC_LIBRARY:        ${STATIC_LIBRARY}" )


#####################################
# Setting the component to be build #
#####################################
# Set the build options
option( BUILD_STATIC_LIBRARY  "Build static library"  ON  )
option( BUILD_APPLICATIONS    "Build applications"    ON  )
option( BUILD_SHARED_LIBRARY  "Build shared library"  OFF )
option( BUILD_EXAMPLES        "Build examples"        OFF )
option( BUILD_TESTS           "Build tests"           OFF )
LogMessage( "BUILD_STATIC_LIBRARY:    ${BUILD_STATIC_LIBRARY}" )
LogMessage( "BUILD_SHARED_LIBRARY:    ${BUILD_SHARED_LIBRARY}" )
LogMessage( "BUILD_EXAMPLES      :    ${BUILD_EXAMPLES}" )
LogMessage( "BUILD_TESTS         :    ${BUILD_TESTS}" )


##########################
# Setting the build type #
##########################
# Two buid types are available:
# - "debug":    Uses debug compiler options
# - "release":  Uses optimized compiler options
# The default build type is "debug".
if (NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE "debug" CACHE STRING "Select the type of build" FORCE )
endif ()
set_property ( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "debug" "release" )
LogMessage( "Build type: CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )


################################
# Setting the Fortran compiler #
################################
# The Fortran compiler is set here. Currently, only intel Fortran compiler is supported.
# If no Fortran compiler has been defined, then ifort is used by default.
# IMPORTRANT: The selection of the Fortran compiler must be done BEFORE 'enable_language(Fortran)'.
set( DEFAULT_Fortran_COMPILER "gfortran" )
if( NOT DEFINED CMAKE_Fortran_COMPILER)
    set( CMAKE_Fortran_COMPILER ${DEFAULT_Fortran_COMPILER} CACHE STRING "" FORCE)
else()
    set( CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER} CACHE STRING "" FORCE)
endif()
enable_language( Fortran )   # This line MUST go after the selection of the Fortran compiler
string( REPLACE "." ";" VERSION_LIST ${CMAKE_Fortran_COMPILER_VERSION} )
list(GET VERSION_LIST 0 Fortran_COMPILER_VERSION_MAJOR)
list(GET VERSION_LIST 1 Fortran_COMPILER_VERSION_MINOR)
list(GET VERSION_LIST 2 Fortran_COMPILER_VERSION_PATCH)
set( CMAKE_Fortran_COMPILER_VERSION "${Fortran_COMPILER_VERSION_MAJOR}.${Fortran_COMPILER_VERSION_MINOR}.${Fortran_COMPILER_VERSION_PATCH}" )
set( Fortran_COMPILER_CONFIG ${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION} )
string( TOLOWER ${Fortran_COMPILER_CONFIG} Fortran_COMPILER_CONFIG )
LogMessage( "Fortran compiler:          CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}" )
LogMessage( "Fortran compiler version:  CMAKE_Fortran_COMPILER_VERSION = ${CMAKE_Fortran_COMPILER_VERSION}" )
LogMessage( "Fortran compiler config:   Fortran_COMPILER_CONFIG = ${Fortran_COMPILER_CONFIG}" )


#############################################
# Setting the name of current configuration #
#############################################
set( PARALLEL_TYPE "" )
if ( PARALLEL_ENVIRONMENT STREQUAL "Coarray" )
  set( PARALLEL_TYPE "-${PARALLEL_ENVIRONMENT}" )
endif ()
set( PROJECT_BUILD_CONFIGURATION ${PROJECT_VERSION}${PARALLEL_TYPE}-${CMAKE_BUILD_TYPE}-${Fortran_COMPILER_CONFIG} )
string( TOLOWER ${PROJECT_BUILD_CONFIGURATION} PROJECT_BUILD_CONFIGURATION )
set( PROJECT_FULL_NAME ${PROJECT_NAME_LC}-${PROJECT_BUILD_CONFIGURATION} )
LogMessage( "PROJECT_BUILD_CONFIGURATION = ${PROJECT_BUILD_CONFIGURATION}" )
LogMessage( "PROJECT_FULL_NAME = ${PROJECT_FULL_NAME}" )


######################################
# Setting the Fortran compiler flags #
######################################
# The Fortran compiler flags are set in the module 'SetFortranFlags' which can be modified
# to meet you needs.
# The Fortran flags sould not be set more than once.  Therefore, we have to define a variable
# which indicates to the sub-project that the Fortran flags have alreay been defined.
include ( SetFortranFlags )


##############################################
# Setting the names of the build directories #
##############################################
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib )  # Set the name of the directory where to store libraries
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin )  # Set the name of the directory where to store executables
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod )  # Set the name of the directory where to store Fortran module files '*.mod'
include_directories( ${CMAKE_Fortran_MODULE_DIRECTORY} )


#####################################################
# Setting the names of the installation directories #
#####################################################
# Setup the directories where to install the code
# if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
#   LogMessage ( "Changing the default installation prefix is used")
#   LogMessage ( "-> from: CMAKE_INSTALL_PREFIX   = ${CMAKE_INSTALL_PREFIX}")
# #   set( CMAKE_INSTALL_PREFIX "/opt/${PROJECT_NAME_LC}/${PROJECT_BUILD_CONFIGURATION}" CACHE PATH "Installation prefix" FORCE  )
# #   set( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/../install/${PROJECT_BUILD_CONFIGURATION}" CACHE PATH "Installation prefix" FORCE  )
#   set( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/../install/Test" CACHE PATH "Installation prefix" FORCE  )
#   LogMessage ( "-> to:   CMAKE_INSTALL_PREFIX   = ${CMAKE_INSTALL_PREFIX}")
# endif()
set( PROJECT_INSTALL_DIR      ${CMAKE_INSTALL_PREFIX}/${PROJECT_FULL_NAME} ) # Installation directory for libraries (both static and shared)
get_filename_component( PROJECT_INSTALL_DIR "${PROJECT_INSTALL_DIR}" ABSOLUTE)
set( LIBRARY_INSTALL_DIR      ${PROJECT_INSTALL_DIR}/lib )          # Installation directory for libraries (both static and shared)
set( MODULES_INSTALL_DIR      ${PROJECT_INSTALL_DIR}/mod )          # Installation directory for Fortran modules
set( EXECUTABLE_INSTALL_DIR   ${PROJECT_INSTALL_DIR}/bin )          # Installation directory for executables
set( CMAKE_INSTALL_DIR        ${PROJECT_INSTALL_DIR}/cmake )        # Installation directory for cmake configuration files
set( DATABASE_INSTALL_DIR     ${PROJECT_INSTALL_DIR}/dtb )          # Installation directory for database files
set( INCLUDE_INSTALL_DIR      ${MODULES_INSTALL_DIR} )              # Installation directory for cmake configuration files
LogMessage ( "Installation directory prefix:                        CMAKE_INSTALL_PREFIX   = ${CMAKE_INSTALL_PREFIX}")
LogMessage ( "Installation directory for the projet:                PROJECT_INSTALL_DIR    = ${PROJECT_INSTALL_DIR}")
LogMessage ( "Installation directory for libraries:                 LIBRARY_INSTALL_DIR    = ${LIBRARY_INSTALL_DIR}")
LogMessage ( "Installation directory for Fortran modules:           MODULES_INSTALL_DIR    = ${MODULES_INSTALL_DIR}")
LogMessage ( "Installation directory for include files:             INCLUDE_INSTALL_DIR    = ${INCLUDE_INSTALL_DIR}")
LogMessage ( "Installation directory for executables:               EXECUTABLE_INSTALL_DIR = ${EXECUTABLE_INSTALL_DIR}")
LogMessage ( "Installation directory for cmake configuration files: CMAKE_INSTALL_DIR      = ${CMAKE_INSTALL_DIR}")
LogMessage ( "Installation directory for database files:            DATABASE_INSTALL_DIR   = ${DATABASE_INSTALL_DIR}")
if ( DEFINED MODULEFILE_INSTALL_DIR )
LogMessage ( "Installation directory for the modulefile:            MODULEFILE_INSTALL_DIR = ${MODULEFILE_INSTALL_DIR}")
endif()

LogMessage ( "GSL  Root Directory:                                   GSL_ROOT_DIR           = ${GSL_ROOT_DIR}")
LogMessage ( "GSL  Include:                                          GSL_INC                = ${GSL_INCLUDE_DIR}")
LogMessage ( "GSL  Library:                                          GSL_LIB                = ${GSL_LIBRARIES}")
#find_package( GSL REQUIRED )
#LogMessage ( "FGSL Include:                                          FGSL_INC               = ${FGSL_INC}")
#LogMessage ( "FGSL Library:                                          FGSL_LIB               = ${FGSL_LIB}")

#LogMessage ( "HDF5 Include:                                          HDF5_INCLUDE_DIRS      = ")
#LogMessage ( "HDF5 Library:                                          HDF5_LIBRARIES         = ${HDF5_LIBRARIES}")
#find_package( HDF5 REQUIRED )

###############################
# Building the static library #
###############################
add_subdirectory( src )                                                             # Adding the directory containing the source files
get_sources( SRCS )                                                                 # Get the list of source files and reseting the global variable
if ( BUILD_STATIC_LIBRARY )
  set( TARGET_NAME ${LIBRARY_NAME}-static )                                         # Set the name of the target to be build: 'spark-static'
  set( STATIC_LIBRARY ${TARGET_NAME} )                                              # Save the name of the target associated to the static library for re-use when building application codes)
  LogMessage( "Building the static library: ${LIBRARY_NAME}" )
  add_library( ${TARGET_NAME} STATIC ${SRCS} )                                      # Create a static library from the list of source files
  set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME ${LIBRARY_NAME})     # Change the name of the static library from 'spark-static' to 'spark'
  install(  TARGETS       ${STATIC_LIBRARY}                                         # Installing the static library
            DESTINATION   ${LIBRARY_INSTALL_DIR}
            EXPORT        ${LIBRARY_NAME}-targets )
endif()


###############################
# Building the shared library #
###############################
if ( BUILD_SHARED_LIBRARY )
  set( TARGET_NAME ${LIBRARY_NAME}-shared )                                         # Set the name of the target to be build: 'spark-shared'
  LogMessage( "Building the shared library: ${LIBRARY_NAME}" )
  add_library( ${TARGET_NAME} SHARED ${SRCS} )                                      # Create a shared library from the list of source files
  set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME ${LIBRARY_NAME})     # Change the name of the shared library from 'spark-shared' to 'spark'
  install(  TARGETS       ${TARGET_NAME}                                            # Install the shared library
            DESTINATION   ${LIBRARY_INSTALL_DIR}
            EXPORT        ${LIBRARY_NAME}-targets )
endif()


##################################
# Building the application codes #
##################################

if (BUILD_APPLICATIONS)
  set( APPLICATION_NAME ${LIBRARY_NAME}-initialize.x )                                # Set the name of the application (ie the executable) to be build
  set( TARGET_NAME ${APPLICATION_NAME} )                                              # Set the name of the target to be build
  LogMessage( "Building the ${APPLICATION_NAME} code" )                               # Write log message
  add_subdirectory( app/Initialize )                                                  # Adding the directory containing the source files
  get_sources( SRCS )                                                                 # Get the list of source files and reseting the global variable
  add_executable( ${TARGET_NAME} ${SRCS} )                                            # Create a application code from the list of source files
  set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME ${APPLICATION_NAME} )  # Change the name of the executable from 'spark-driver' to 'spark'
  target_link_libraries( ${TARGET_NAME} PUBLIC ${STATIC_LIBRARY} )                    # Linking the 'spark' executable to the 'spark-cfd' library, which will also link to the 'spark' and 'spark-ode' libraries due to the PUBLIC dependency
  install( TARGETS ${TARGET_NAME} DESTINATION ${EXECUTABLE_INSTALL_DIR} )             # Installing the exectuable
endif()

if (BUILD_APPLICATIONS)
  set( APPLICATION_NAME ${LIBRARY_NAME}-computelevels.x )                             # Set the name of the application (ie the executable) to be build
  set( TARGET_NAME ${APPLICATION_NAME} )                                              # Set the name of the target to be build
  LogMessage( "Building the ${APPLICATION_NAME} code" )                               # Write log message
  add_subdirectory( app/ComputeLevels )                                               # Adding the directory containing the source files
  get_sources( SRCS )                                                                 # Get the list of source files and reseting the global variable
  add_executable( ${TARGET_NAME} ${SRCS} )                                            # Create a application code from the list of source files
  set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME ${APPLICATION_NAME} )  # Change the name of the executable from 'spark-driver' to 'spark'
  target_link_libraries( ${TARGET_NAME} PUBLIC ${STATIC_LIBRARY} )                    # Linking the 'spark' executable to the 'spark-cfd' library, which will also link to the 'spark' and 'spark-ode' libraries due to the PUBLIC dependency
  install( TARGETS ${TARGET_NAME} DESTINATION ${EXECUTABLE_INSTALL_DIR} )             # Installing the exectuable
endif()

if (BUILD_APPLICATIONS)
  set( APPLICATION_NAME ${LIBRARY_NAME}-mergelevels.x )                               # Set the name of the application (ie the executable) to be build
  set( TARGET_NAME ${APPLICATION_NAME} )                                              # Set the name of the target to be build
  LogMessage( "Building the ${APPLICATION_NAME} code" )                               # Write log message
  add_subdirectory( app/MergeLevels )                                                 # Adding the directory containing the source files
  get_sources( SRCS )                                                                 # Get the list of source files and reseting the global variable
  add_executable( ${TARGET_NAME} ${SRCS} )                                            # Create a application code from the list of source files
  set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME ${APPLICATION_NAME} )  # Change the name of the executable from 'spark-driver' to 'spark'
  target_link_libraries( ${TARGET_NAME} PUBLIC ${STATIC_LIBRARY} )                    # Linking the 'spark' executable to the 'spark-cfd' library, which will also link to the 'spark' and 'spark-ode' libraries due to the PUBLIC dependency
  install( TARGETS ${TARGET_NAME} DESTINATION ${EXECUTABLE_INSTALL_DIR} )             # Installing the exectuable
endif()

if (BUILD_APPLICATIONS)
  set( APPLICATION_NAME ${LIBRARY_NAME}-preproclevels.x )                             # Set the name of the application (ie the executable) to be build
  set( TARGET_NAME ${APPLICATION_NAME} )                                              # Set the name of the target to be build
  LogMessage( "Building the ${APPLICATION_NAME} code" )                               # Write log message
  add_subdirectory( app/PreprocLevels )                                               # Adding the directory containing the source files
  get_sources( SRCS )                                                                 # Get the list of source files and reseting the global variable
  add_executable( ${TARGET_NAME} ${SRCS} )                                            # Create a application code from the list of source files
  set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME ${APPLICATION_NAME} )  # Change the name of the executable from 'spark-driver' to 'spark'
  target_link_libraries( ${TARGET_NAME} PUBLIC ${STATIC_LIBRARY} )                    # Linking the 'spark' executable to the 'spark-cfd' library, which will also link to the 'spark' and 'spark-ode' libraries due to the PUBLIC dependency
  install( TARGETS ${TARGET_NAME} DESTINATION ${EXECUTABLE_INSTALL_DIR} )             # Installing the exectuable
endif()

if (BUILD_APPLICATIONS)
  set( APPLICATION_NAME ${LIBRARY_NAME}-runtrajectories.x )                           # Set the name of the application (ie the executable) to be build
  set( TARGET_NAME ${APPLICATION_NAME} )                                              # Set the name of the target to be build
  LogMessage( "Building the ${APPLICATION_NAME} code" )                               # Write log message
  add_subdirectory( app/RunTrajectories )                                             # Adding the directory containing the source files
  get_sources( SRCS )                                                                 # Get the list of source files and reseting the global variable
  add_executable( ${TARGET_NAME} ${SRCS} )                                            # Create a application code from the list of source files
  set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME ${APPLICATION_NAME} )  # Change the name of the executable from 'spark-driver' to 'spark'
  target_link_libraries( ${TARGET_NAME} PUBLIC ${STATIC_LIBRARY} )                    # Linking the 'spark' executable to the 'spark-cfd' library, which will also link to the 'spark' and 'spark-ode' libraries due to the PUBLIC dependency
  install( TARGETS ${TARGET_NAME} DESTINATION ${EXECUTABLE_INSTALL_DIR} )             # Installing the exectuable
endif()

if (BUILD_APPLICATIONS)
  set( APPLICATION_NAME ${LIBRARY_NAME}-trajectoriesstats.x )                         # Set the name of the application (ie the executable) to be build
  set( TARGET_NAME ${APPLICATION_NAME} )                                              # Set the name of the target to be build
  LogMessage( "Building the ${APPLICATION_NAME} code" )                               # Write log message
  add_subdirectory( app/TrajectoriesStats )                                           # Adding the directory containing the source files
  get_sources( SRCS )                                                                 # Get the list of source files and reseting the global variable
  add_executable( ${TARGET_NAME} ${SRCS} )                                            # Create a application code from the list of source files
  set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME ${APPLICATION_NAME} )  # Change the name of the executable from 'spark-driver' to 'spark'
  target_link_libraries( ${TARGET_NAME} PUBLIC ${STATIC_LIBRARY} )                    # Linking the 'spark' executable to the 'spark-cfd' library, which will also link to the 'spark' and 'spark-ode' libraries due to the PUBLIC dependency
  install( TARGETS ${TARGET_NAME} DESTINATION ${EXECUTABLE_INSTALL_DIR} )             # Installing the exectuable
endif()

if (BUILD_APPLICATIONS)
  set( APPLICATION_NAME ${LIBRARY_NAME}-posttrajectories.x )                          # Set the name of the application (ie the executable) to be build
  set( TARGET_NAME ${APPLICATION_NAME} )                                              # Set the name of the target to be build
  LogMessage( "Building the ${APPLICATION_NAME} code" )                               # Write log message
  add_subdirectory( app/PostTrajectories )                                            # Adding the directory containing the source files
  get_sources( SRCS )                                                                 # Get the list of source files and reseting the global variable
  add_executable( ${TARGET_NAME} ${SRCS} )                                            # Create a application code from the list of source files
  set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME ${APPLICATION_NAME} )  # Change the name of the executable from 'spark-driver' to 'spark'
  target_link_libraries( ${TARGET_NAME} PUBLIC ${STATIC_LIBRARY} )                    # Linking the 'spark' executable to the 'spark-cfd' library, which will also link to the 'spark' and 'spark-ode' libraries due to the PUBLIC dependency
  install( TARGETS ${TARGET_NAME} DESTINATION ${EXECUTABLE_INSTALL_DIR} )             # Installing the exectuable
endif()

# if (BUILD_APPLICATIONS)
#  LogMessage ( "GSL  Include:           GSL_INC                = ${GSL_INCLUDE_DIR}")
#  LogMessage ( "GSL  Library:           GSL_LIB                = ${GSL_LIBRARIES}")
#  set( APPLICATION_NAME ${LIBRARY_NAME}-derivequantities.x )                          # Set the name of the application (ie the executable) to be build
#  set( TARGET_NAME ${APPLICATION_NAME} )                                              # Set the name of the target to be build
#  LogMessage( "Building the ${APPLICATION_NAME} code" )                               # Write log message
#  add_subdirectory( app/DeriveQuantities )                                            # Adding the directory containing the source files
#  get_sources( SRCS )                                                                 # Get the list of source files and reseting the global variable
#  add_executable( ${TARGET_NAME} ${SRCS} )                                            # Create a application code from the list of source files
#  set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME ${APPLICATION_NAME} )  # Change the name of the executable from 'spark-driver' to 'spark'
#  target_link_libraries( ${TARGET_NAME} PUBLIC ${STATIC_LIBRARY} ${GSL_LIBRARIES} )   # Linking the 'spark' executable to the 'spark-cfd' library, which will #also link to the 'spark' and 'spark-ode' libraries due to the PUBLIC dependency
#  target_include_directories( ${TARGET_NAME} PUBLIC ${GSL_INCLUDE_DIR}  )
#  install( TARGETS ${TARGET_NAME} DESTINATION ${EXECUTABLE_INSTALL_DIR} )             # Installing the exectuable
# endif()

if (BUILD_APPLICATIONS)
  set( APPLICATION_NAME ${LIBRARY_NAME}-plotpes.x )                                   # Set the name of the application (ie the executable) to be build
  set( TARGET_NAME ${APPLICATION_NAME} )                                              # Set the name of the target to be build
  LogMessage( "Building the ${APPLICATION_NAME} code" )                               # Write log message
  add_subdirectory( app/PlotPES )                                                     # Adding the directory containing the source files
  get_sources( SRCS )                                                                 # Get the list of source files and reseting the global variable
  add_executable( ${TARGET_NAME} ${SRCS} )                                            # Create a application code from the list of source files
  set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME ${APPLICATION_NAME} )  # Change the name of the executable from 'spark-driver' to 'spark'
  target_link_libraries( ${TARGET_NAME} PUBLIC ${STATIC_LIBRARY} ${GSL_LIBRARIES} )         # Linking the 'spark' executable to the 'spark-cfd' library, which will also link to the 'spark' and
  target_include_directories( ${TARGET_NAME} PUBLIC ${GSL_INCLUDE_DIR} )
  install( TARGETS ${TARGET_NAME} DESTINATION ${EXECUTABLE_INSTALL_DIR} )             # Installing the exectuable
endif()


##############################################################
# Installing the Fortran module files and the database files #
##############################################################
# if ( BUILD_STATIC_LIBRARY OR BUILD_SHARED_LIBRARY )
#   set( DATABASE_SOURCE_DIR ${PROJECT_SOURCE_DIR}/dtb )                                          # Directory containing the database files within the sources
#   install( DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/   DESTINATION ${MODULES_INSTALL_DIR}  ) # Installing the Fortran module files
#   install( DIRECTORY ${DATABASE_SOURCE_DIR}/               DESTINATION ${DATABASE_INSTALL_DIR} ) # Installing the kinetic files
# endif()


##################################################
# Creating the modulefile in the build directory #
##################################################
configure_file(
    ${PROJECT_SOURCE_DIR}/config/modulefile.in
    ${CMAKE_BINARY_DIR}/${PROJECT_FULL_NAME}.modulefile )
if ( DEFINED MODULEFILE_INSTALL_DIR )
  install(  FILES         ${CMAKE_BINARY_DIR}/${PROJECT_FULL_NAME}.modulefile
            DESTINATION   ${MODULEFILE_INSTALL_DIR}
            RENAME        ${PROJECT_BUILD_CONFIGURATION} )
endif()


#################################################################
# Creating and installing the Cmake Package Configuration Files #
#################################################################
# The package configuration and version files are required so that
# other projects can locate the current project by using "find_package".
include(CMakePackageConfigHelpers)

# We configure our template.
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config/pkg/config.cmake.in
  ${PROJECT_BINARY_DIR}/config/pkg/config.cmake
  INSTALL_DESTINATION     ${CMAKE_INSTALL_DIR}
  PATH_VARS               LIBRARY_INSTALL_DIR MODULES_INSTALL_DIR INCLUDE_INSTALL_DIR CMAKE_INSTALL_DIR
)

write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/config/pkg/config-version.cmake
  VERSION                 1.0.0
  COMPATIBILITY           AnyNewerVersion
)

install(  EXPORT        ${LIBRARY_NAME}-targets
          DESTINATION   ${CMAKE_INSTALL_DIR} )

install(  FILES         ${PROJECT_BINARY_DIR}/config/pkg/config.cmake
          DESTINATION   ${CMAKE_INSTALL_DIR}
          RENAME        ${LIBRARY_NAME}-config.cmake )

install(  FILES         ${PROJECT_BINARY_DIR}/config/pkg/config-version.cmake
          DESTINATION   ${CMAKE_INSTALL_DIR}
          RENAME        ${LIBRARY_NAME}-config-version.cmake )
