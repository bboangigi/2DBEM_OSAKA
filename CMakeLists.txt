## -------------------------------------------------------------------------- ##
#
#   Project: Example
#
#   Description:
#       CMake & Fortran project for graduate students
#
#   Author:
#       Young-Myung CHOI
#         Dept. of Naval Arch. and Ocean Eng.,
#         Pusan National University, Republic of Korea
#
## -------------------------------------------------------------------------- ##

## Default Settings --------------------------------------------------------- ##

##... Set the minimum version of cmake
cmake_minimum_required(VERSION 2.9 FATAL_ERROR)

##... CMake Project Name
project(FortranProjectBase)

##... Set Default Fortran Compiler & Eable the Fortran Compilation Environment
enable_language(Fortran)
if (NOT DEFINED CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER gfortran)
endif()

## CMake Configuration ------------------------------------------------------ ##

##... Load Default CMake Configuration
include( config/config.cmake )

##... Set Fortran Compiler Flags
include( config/fortranFlag.cmake )

##... Auto Compilation Functions
include( config/Func_AutoCompile.cmake )

## Print CMake Project Information ------------------------------------------ ##
message("")
message("${BoldYellow}Project: ${BoldWhite}${CMAKE_PROJECT_NAME}")
message("")
message("  - ${Green}Projct Path       ${White}: ${PROJECT_COMPILE_PATH} ${Green}")
message("  - ${Green}GIT_BRANCH        ${White}: ${GIT_BRANCH} ${Green}")
message("  - ${Green}GIT_COMMIT_HASH   ${White}: ${GIT_COMMIT_HASH} ${Green}")
message("")
## SRC 폴더의 CMakeLists.txt와 연관 있음.  
message("  - ${Green}Lib. Install Path ${White}: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${Green}")
message("  - ${Green}Exe. Install Path ${White}: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${Green}")
message("")

# set(FORTRAN_COMPILE_MODE "Release")
# set(FORTRAN_COMPILE_MODE "Debug")

if( ${FORTRAN_COMPILE_MODE} STREQUAL "Release" )
message("  - ${Green}Compilation Mode  ${White}: Release ${Green}")
elseif( ${FORTRAN_COMPILE_MODE} STREQUAL "Debug" )
message("  - ${Green}Compilation Mode  ${White}: Debug ${Green}")
endif()

message("${Green}")
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
message("  - ${Green}Fortran Compiler  ${White}: Intel Fortran ${Green}")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
message("  - ${Green}Fortran Compiler  ${White}: GNU Fortran ${Green}")
endif()
message("  - ${Green}Global Flag       ${White}: ${CMAKE_Fortran_FLAGS} ${Green}")
message("${ColourReset}")
## Compile Library ---------------------------------------------------------- ##

#... Library call folder name
set(LIB_NAME_Libraries_2DBEM_OSAKA lib2DBEM_OSAKA)
add_subdirectory(src)

## Compile Applications ----------------------------------------------------- ##

#... Application
set(APP_NAME_app2DBEM app2DBEM)
set(APP_EXTLIB_appTutorial ${LAPACK_FLAG})
set(APP_DEPEND_app2DBEM lib2DBEM_OSAKA)  # --- app2DBEM-library connection folder name / commented when library not use
add_subdirectory(applications/app2DBEM)