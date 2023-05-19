# Check CMake version.
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)

# ---Build executables--------------------------------------------------------------------
add_executable(make_track make_track.cpp)
add_executable(gas_table gas_table.cpp)

target_link_libraries(make_track Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_link_libraries(gas_table Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)

target_include_directories(make_track PRIVATE ../../include)
target_include_directories(gas_table PRIVATE ../../include)