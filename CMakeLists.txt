# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

add_compile_options(-O3) # enable maximum optimization

# ---Root---------------------------------------------------------------------------------
message(STATUS "Looking for ROOT")
find_package(ROOT REQUIRED)
message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}")

set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

find_package(GSL REQUIRED)

# ---Add subdirectory containing source code----------------------------------------------
include_directories(include)
file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp")
add_library(X17_source_lib STATIC ${SOURCES})

# ---Build executables--------------------------------------------------------------------
# Create the build folder for the root project
set(PROJECT_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build)
file(MAKE_DIRECTORY ${PROJECT_BUILD_DIR})

# Create a list of subdirectories (assuming 'subdir1' and 'subdir2')
set(SUBDIRECTORIES reconstruction simulations/ion_map simulations/single_track simulations/rk_tracks)

# Iterate over each subdirectory
foreach(subdirectory ${SUBDIRECTORIES})
    # Add the subdirectory
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${subdirectory})
endforeach()