# Overrides the Add_Eigen from MathUtils as we need to define the variable
# EIGEN3_INCLUDE_DIRS for HighFive

include(FetchContent)

FetchContent_Declare(Eigen3
        GIT_REPOSITORY ${eigen_URL}
        GIT_TAG ${eigen_TAG}
        )

FetchContent_GetProperties(Eigen3)
if(NOT eigen3_POPULATED)
    message(STATUS "******* FETCHING Eigen3 dependency from ${PROJECT_NAME} (requested version: ${eigen_TAG}) *******")
    FetchContent_Populate(Eigen3)
    add_library(eigen INTERFACE)
    target_include_directories(eigen INTERFACE ${eigen3_SOURCE_DIR})
    add_library(Eigen3::Eigen ALIAS eigen)
endif()

set(EIGEN3_INCLUDE_DIRS ${eigen3_SOURCE_DIR})  # Necessary for HighFive
