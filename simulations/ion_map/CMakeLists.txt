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

# find_package(GSL REQUIRED)
include_directories(/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/include)
link_directories(/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/lib)

# ---Build executables--------------------------------------------------------------------
add_executable(ion_electrons ion_electrons.cpp)
add_executable(make_map make_map.cpp)

target_link_libraries(ion_electrons Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_link_libraries(make_map Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_include_directories(ion_electrons PRIVATE ../../include)
target_include_directories(make_map PRIVATE ../../include)