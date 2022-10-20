add_library(hdb5_io)

target_sources(hdb5_io PRIVATE
        containers/HydrodynamicDataBase.cpp
        containers/Body.cpp
        containers/Mask.cpp
        containers/WaveDrift.cpp
        containers/PoleResidue.cpp
        containers/Kochin.cpp
        io/HDBReader.cpp
        io/HDBWriter.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
        )

if (HDB5IO_MESH_SUPPORT)
    target_sources(hdb5_io PRIVATE
            containers/Mesh.cpp
            )
endif ()

target_include_directories(hdb5_io PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)
target_link_libraries(hdb5_io PUBLIC
        MathUtils::MathUtils
        eigen
        HighFive
        )

if (HDB5IO_MESH_SUPPORT)
    target_link_libraries(hdb5_io PUBLIC meshoui)
endif()

if (BUILD_SHARED_LIBS)
    message(STATUS "Using HDF5 SHARED lib into hdb5_io")
    target_link_libraries(hdb5_io PUBLIC hdf5_cpp-shared)
else ()
    message(STATUS "Using HDF5 STATIC lib into hdb5_io")
    target_link_libraries(hdb5_io PUBLIC hdf5_cpp-static)
    set_property(TARGET hdb5_io PROPERTY POSITION_INDEPENDENT_CODE ON)
endif ()

add_dependencies(hdb5_io check_git_${PROJECT_NAME}) # For git_watcher to fetch git informations before effective build

if (HDB5IO_USE_VTK)
    target_link_libraries(hdb5_io PUBLIC ${VTK_LIBRARIES})
    target_compile_definitions(hdb5_io PUBLIC -DOMPI_SKIP_MPICXX)
    target_compile_definitions(hdb5_io PUBLIC -DUSE_VTK)
endif ()

if (HDB5IO_MESH_SUPPORT)
    target_compile_definitions(hdb5_io PUBLIC MESH_SUPPORT)
endif ()
