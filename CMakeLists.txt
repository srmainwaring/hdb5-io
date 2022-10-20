#=============================================================================
#
# CMake configuration file for the HDB5_IO project
#
#=============================================================================
cmake_minimum_required(VERSION 3.11)

project(hdb5_io)
set(CMAKE_CXX_STANDARD 17)

#-----------------------------------------------------------------------------
# Prevent in-source builds
#-----------------------------------------------------------------------------
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source build is not possible and not recommended. Choose an empty directory for build output.")
endif (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

#-----------------------------------------------------------------------------
# OPTIONS
#-----------------------------------------------------------------------------
cmake_policy(SET CMP0077 NEW)

option(BUILD_SHARED_LIBS "Shared lib version" OFF)
option(HDB5IO_USE_VTK "Enable VTK Support" OFF)
option(HDB5IO_MESH_SUPPORT "Build with mesh support from meshoui" ON)
option(HDB5IO_BUILD_TESTS "Activate build tests" ON)

message(STATUS "SHARED_LIB set to ${BUILD_SHARED_LIBS} in hdb5_io")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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
endif ()

# =============================================================================
# HDB5_IO DEPENDENCIES
# =============================================================================
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/URL.conf.cmake)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
if (HDB5IO_USE_VTK)
    include(Add_VTK)
endif ()

include(Add_Eigen3)
include(Add_mathutils)
include(Add_zlib)
include(Add_hdf5)
include(Add_HighFive)
if (HDB5IO_MESH_SUPPORT)
    include(Add_meshoui)
endif()

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
