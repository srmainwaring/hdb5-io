#=============================================================================
#
# CMake configuration file for the HDB5_IO project
#
#=============================================================================


#-----------------------------------------------------------------------------
# Check backward compatibility
#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.11)

#-----------------------------------------------------------------------------
# Prevent in-source builds
#-----------------------------------------------------------------------------
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source build is not possible and not recommended. Choose an empty directory for build output.")
endif (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

#-----------------------------------------------------------------------------
# Define project name and specify version information
#-----------------------------------------------------------------------------
project(hdb5_io)

#-----------------------------------------------------------------------------
# Setting the configuration of the build (Debug, Release...)
#-----------------------------------------------------------------------------
#set(CMAKE_BUILD_TYPE Debug)

# OPTIONS
set(USE_VTK OFF CACHE BOOL "Enable VTK Support")
mark_as_advanced(USE_VTK)

option(HDB5IO_USE_VTK "Enable VTK Support" ${USE_VTK})
option(HDB5IO_BUILD_TESTS "Activate build tests" ON)

#-----------------------------------------------------------------------------
# C++ dialect
#-----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 14)


# Access to local cmake scripts
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#=============================================================================
# Retrieving the current Git revision
#=============================================================================
set(PRE_CONFIGURE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/src/hdb5_io/version.cpp.in CACHE FILEPATH "" FORCE)
set(POST_CONFIGURE_FILE ${CMAKE_CURRENT_BINARY_DIR}/src/hdb5_io/version.cpp CACHE FILEPATH "" FORCE)
set(GIT_WORKING_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE FILEPATH "" FORCE)
set(GIT_STATE_FILE ${CMAKE_CURRENT_BINARY_DIR}/git-state-hash CACHE FILEPATH "" FORCE)
include(cmake/git_watcher.cmake)

if (NOT EXISTS ${POST_CONFIGURE_FILE})
    configure_file(${PRE_CONFIGURE_FILE} ${POST_CONFIGURE_FILE})
endif()

# =============================================================================
# HDB5_IO DEPENDENCIES
# =============================================================================

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/URL.conf.cmake)

include(Add_HDF5)
include(Add_mathutils)
include(Add_meshoui)

message(STATUS "USE VTK : ${HDB5IO_USE_VTK}" )

if (HDB5IO_USE_VTK)
    find_package(VTK 6.0 REQUIRED NO_MODULE)
    include(${VTK_USE_FILE})
endif()

get_target_property(INC eigen INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Found eigen : ${INC}")

set(EIGEN3_INCLUDE_DIRS ${INC})
#set(EIGEN3_DIR ${INC})

include(Add_HighFive)

#=============================================================================
# Adding the HDB5_IO project sources
#=============================================================================
add_subdirectory(src/hdb5_io)

#=============================================================================
# Adding tests
#=============================================================================
if (${HDB5IO_BUILD_TESTS})
    enable_testing()
    include(Add_googletest)
    add_subdirectory(tests)
endif ()


#=============================================================================
# Adding tools
#=============================================================================
#if (${HDB5IO_build_tools})
#    include(add_cxxopts)
#    add_subdirectory(tools)
#endif ()
