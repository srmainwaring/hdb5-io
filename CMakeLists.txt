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
project(HDB5_IO)

#-----------------------------------------------------------------------------
# Setting the configuration of the build (Debug, Release...)
#-----------------------------------------------------------------------------
#set(CMAKE_BUILD_TYPE Debug)

option(HDB5IO_build_tests "Activate build tests" ON)

#-----------------------------------------------------------------------------
# C++ dialect
#-----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 14)


# Access to local cmake scripts
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#=============================================================================
# Retrieving the current Git revision
#=============================================================================
#include(GetGitRevisionDescription)
#get_git_head_revision(GIT_REFSPEC GIT_SHA1)
#set(FRYDOM_FLAVOR "FRyDoM-CE (Community Edition)" CACHE STRING "The FRyDoM flavor (CE or EE)")
#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/frydom/version.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/src/frydom/version.cpp")

# =============================================================================
# HDB5_IO DEPENDENCIES
# =============================================================================

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/URL.conf.cmake)

include(Add_MathUtils)
include(Add_HDF5)
include(Add_meshoui)

find_package(VTK 6.0 REQUIRED NO_MODULE)
include(${VTK_USE_FILE})


get_target_property(INC eigen INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Found eigen : ${INC}")

set(Eigen3_DIR ${INC})

include(Add_HighFive)

#=============================================================================
# Adding the HDB5_IO project sources
#=============================================================================
add_subdirectory(src/hdb5_io)

#=============================================================================
# Adding tests
#=============================================================================
if (${HDB5IO_build_tests})
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
