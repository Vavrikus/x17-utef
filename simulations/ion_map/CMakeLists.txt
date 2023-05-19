# Check CMake version.
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)

# find_package(GSL REQUIRED)
include_directories(/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/include)
link_directories(/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/lib)

# ---Build executables--------------------------------------------------------------------
add_executable(ion_electrons ion_electrons.cpp)
add_executable(make_map make_map.cpp)

target_link_libraries(ion_electrons Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_link_libraries(make_map Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib ROOT::RIO)
target_include_directories(ion_electrons PRIVATE ../../include)
target_include_directories(make_map PRIVATE ../../include)