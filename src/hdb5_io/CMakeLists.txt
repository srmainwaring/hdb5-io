
# TODO : ajouter les fichiers source cpp dans la variable suivante :
set(HDB5IO_SRC_FILES
        HydrodynamicDataBase.cpp
        Body.cpp
        Discretization1D.cpp
        Mask.cpp
        )


add_library(hdb5io SHARED)

target_sources(hdb5io PRIVATE ${HDB5IO_SRC_FILES})
target_link_libraries(hdb5io PUBLIC
        MathUtils::MathUtils
        HDF5
        eigen
        HighFive
        meshoui
#        ${VTK_LIBRARIES}
        )

if (HDB5IO_USE_VTK)
    target_link_libraries(hdb5io PUBLIC ${VTK_LIBRARIES})
endif()


target_include_directories(hdb5io PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)
