
#add_executable(example example.cpp)
#target_link_libraries(example HighFive)

#add_executable(example_Eigen example_Eigen.cpp)
#target_link_libraries(example_Eigen HDF5 eigen highfive)

#add_executable(test_Eigen test_Eigen.cpp)
#target_link_libraries(test_Eigen HDF5 eigen highfive)

add_executable(test_import_export test_import_export.cpp)
target_link_libraries(test_import_export hdb5_io)

if (HDB5IO_USE_VTK)
    target_compile_definitions(test_import_export PUBLIC -DH5_USE_VTK)
endif()

set_target_properties( test_import_export
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/test_hdb5io"
        )