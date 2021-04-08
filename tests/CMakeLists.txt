
add_executable(test_import_export test_import_export.cpp)
target_link_libraries(test_import_export hdb5_io gtest gtest_main)

if (USE_VTK)
    target_compile_definitions(test_import_export PUBLIC -DUSE_VTK)
endif ()

add_test(NAME test_import_export COMMAND test_import_export)

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../data/Boxbarge_Vertices_353_Faces_652.hdb5"
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})