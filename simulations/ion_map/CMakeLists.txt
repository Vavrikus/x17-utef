# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

add_compile_options(-O3) # enable maximum optimization

# ---Root------------------------------------------------------------------------
message(STATUS "Looking for ROOT")
find_package(ROOT REQUIRED)
message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}")

set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

find_package(GSL REQUIRED)

# ---Add subdirectory containing source code----------------------------------------
include_directories(../../include)
file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/../../source/*.cpp")
add_library(X17_source_lib STATIC ${SOURCES})

# ---Build executable------------------------------------------------------------
add_executable(ion_electrons ion_electrons.cpp)

target_link_libraries(ion_electrons Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_include_directories(ion_electrons PRIVATE ../../include)