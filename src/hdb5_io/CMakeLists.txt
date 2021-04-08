if (BUILD_SHARED_LIBS)
    add_library(hdb5_io SHARED)
else ()
    add_library(hdb5_io STATIC)
endif ()

target_sources(hdb5_io PRIVATE
        containers/HydrodynamicDataBase.cpp
        containers/Body.cpp
        containers/Mask.cpp
        containers/Mesh.cpp
        containers/WaveDrift.cpp
        containers/PoleResidue.cpp
        containers/Kochin.cpp
        io/HDBReader.cpp
        io/HDBWriter.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
        )

target_include_directories(hdb5_io PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)
target_link_libraries(hdb5_io PUBLIC
        MathUtils::MathUtils
        eigen
        HighFive
        meshoui
        )

if (BUILD_SHARED_LIBS)
    target_link_libraries(hdb5_io PUBLIC hdf5_cpp-shared)
else ()
    target_link_libraries(hdb5_io PUBLIC hdf5_cpp-static)
endif ()

add_dependencies(hdb5_io check_git_${PROJECT_NAME}) # For git_watcher to fetch git informations before effective build

if (USE_VTK)
    target_link_libraries(hdb5_io PUBLIC ${VTK_LIBRARIES})
    target_compile_definitions(hdb5_io PUBLIC -DOMPI_SKIP_MPICXX)
    target_compile_definitions(hdb5_io PUBLIC -DUSE_VTK)
endif ()
