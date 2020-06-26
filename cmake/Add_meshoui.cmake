
#find_package(meshoui QUIET)

if (NOT meshoui_FOUND)
    include(FetchContent)


    FetchContent_Declare(meshoui
            GIT_REPOSITORY ${meshoui_URL}
            GIT_TAG ${meshoui_TAG}
            )

    FetchContent_GetProperties(meshoui)
    if(NOT meshoui_POPULATED)
        message(STATUS "Downloading, Configuring and Generating 'meshoui' dependency")
        FetchContent_Populate(meshoui)

        # meshoui BUILD OPTIONS
        set(MESHOUI_BUILD_TESTS OFF CACHE BOOL "" FORCE)
        set(MESHOUI_USE_VTK ${HDB5IO_USE_VTK} CACHE BOOL "" FORCE)

        add_subdirectory(${meshoui_SOURCE_DIR} ${meshoui_BINARY_DIR})
    else()
        message(STATUS "meshoui already populated")
    endif()
endif()

if (TARGET meshoui)
    get_target_property(INC meshoui INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS "Found meshoui : ${INC}")
else()
    message(STATUS "meshoui target NOT FOUND")
endif()
