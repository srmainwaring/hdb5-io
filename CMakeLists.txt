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
option(USE_VTK "Enable VTK Support" ON)
option(HDB5IO_BUILD_TESTS "Activate build tests" ON)

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
if (USE_VTK)
    include(Add_VTK)
endif ()

include(Add_mathutils)
include(Add_HDF5)
include(Add_HighFive)
include(Add_meshoui)

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
