include(FetchContent)

FetchContent_Declare(highfive
        GIT_REPOSITORY ${highfive_URL}
        GIT_TAG ${highfive_TAG}
        )

FetchContent_GetProperties(highfive)
if (NOT highfive_POPULATED)
    message(STATUS "******* FETCHING highfive dependency from ${PROJECT_NAME} (requested version: ${highfive_TAG}) *******")
    FetchContent_Populate(highfive)

    # Highfive BUILD OPTIONS
    set(HIGHFIVE_USE_BOOST OFF CACHE BOOL "")
    set(HIGHFIVE_USE_EIGEN ON CACHE BOOL "")
    set(HIGHFIVE_USE_XTENSOR OFF CACHE BOOL "")
    set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "")
    set(HIGHFIVE_USE_OPENCV OFF CACHE BOOL "")
    set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "")
    set(HIGHFIVE_PARALLEL_HDF5 OFF CACHE BOOL "")
    set(HIGHFIVE_BUILD_DOCS OFF CACHE BOOL "")

    add_library(HighFive INTERFACE)
    target_compile_definitions(HighFive INTERFACE HIGHFIVE_USE_EIGEN=ON)
    target_include_directories(HighFive INTERFACE ${highfive_SOURCE_DIR}/include)
endif ()

