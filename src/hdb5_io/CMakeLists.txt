
# TODO : ajouter les fichiers source cpp dans la variable suivante :
set(HDB5IO_SRC_FILES
        containers/HydrodynamicDataBase.cpp
        containers/Body.cpp
        containers/Mask.cpp
        containers/Mesh.cpp
        containers/WaveDrift.cpp
        containers/PoleResidue.cpp
        io/HDBReader.cpp
        io/HDBWriter.cpp
        )


add_library(hdb5io SHARED)

target_sources(hdb5io PRIVATE ${HDB5IO_SRC_FILES})
target_link_libraries(hdb5io PUBLIC
        MathUtils::MathUtils
        HDF5
        eigen
        HighFive
        meshoui
        )
target_compile_definitions(hdb5io PUBLIC -DOMPI_SKIP_MPICXX)

if (HDB5IO_USE_VTK)
    target_link_libraries(hdb5io PUBLIC ${VTK_LIBRARIES})
endif ()


target_include_directories(hdb5io PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)
