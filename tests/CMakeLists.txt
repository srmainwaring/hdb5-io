
#add_executable(example example.cpp)
#target_link_libraries(example HighFive)

#add_executable(example_Eigen example_Eigen.cpp)
#target_link_libraries(example_Eigen HDF5 eigen highfive)

#add_executable(test_Eigen test_Eigen.cpp)
#target_link_libraries(test_Eigen HDF5 eigen highfive)

add_executable(test_import test_import_export.cpp)
target_link_libraries(test_import hdb5io)
target_compile_definitions(test_import PUBLIC -DOMPI_SKIP_MPICXX)
