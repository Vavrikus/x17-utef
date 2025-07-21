# Check CMake version.
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)

# ---Build executables--------------------------------------------------------------------
add_executable(reco_mtracks reco_mtracks.cpp)
add_executable(reco_plots reco_plots.cpp)
add_executable(reco_test reco_test.cpp)
add_executable(map_test map_test.cpp)

target_link_libraries(reco_mtracks Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_link_libraries(reco_plots Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_link_libraries(reco_test Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_link_libraries(map_test Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)

target_include_directories(reco_mtracks PRIVATE ../include)
target_include_directories(reco_plots PRIVATE ../include)
target_include_directories(reco_test PRIVATE ../include)
target_include_directories(map_test PRIVATE ../include)