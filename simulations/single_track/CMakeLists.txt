# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Root---------------------------------------------------------------------------------
message(STATUS "Looking for ROOT")
find_package(ROOT REQUIRED)
message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}")

set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

find_package(GSL REQUIRED)

# ---Build executables--------------------------------------------------------------------
add_executable(make_track make_track.cpp)
add_executable(gas_table gas_table.cpp)

target_link_libraries(make_track Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_link_libraries(gas_table Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)

target_include_directories(make_track PRIVATE ../../include)
target_include_directories(gas_table PRIVATE ../../include)