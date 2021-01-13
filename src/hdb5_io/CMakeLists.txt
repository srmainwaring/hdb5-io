
# TODO : ajouter les fichiers source cpp dans la variable suivante :
set(HDB5IO_SRC_FILES
        containers/HydrodynamicDataBase.cpp
        containers/Body.cpp
        containers/Mask.cpp
        containers/Mesh.cpp
        containers/WaveDrift.cpp
        containers/PoleResidue.cpp
        containers/Kochin.cpp
        io/HDBReader.cpp
        io/HDBWriter.cpp
        )


add_library(hdb5io SHARED)

target_sources(hdb5io PRIVATE ${HDB5IO_SRC_FILES})
target_link_libraries(hdb5io PUBLIC
        MathUtils::MathUtils
#        HDF5  # when using the findHDF5 script into the Add_HDF5.cmake file
        hdf5_cpp-static
        eigen
        HighFive
        meshoui
        )

if (HDB5IO_USE_VTK)
    target_link_libraries(hdb5io PUBLIC ${VTK_LIBRARIES})
    target_compile_definitions(hdb5io PUBLIC -DOMPI_SKIP_MPICXX)
    target_compile_definitions(hdb5io PUBLIC -DH5_USE_VTK)
endif ()


target_include_directories(hdb5io PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)

set_target_properties( hdb5io
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/hdb5io"
        )