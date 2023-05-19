# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)

# Find the Garfield package
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED HINTS "/home/vavrik/work/source/garfield/build/CMake")
endif()
# ########################################################################################

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Enable debug symbols and disable optimization for Debug build type
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
endif()

# Enable optimization for Release build type
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
endif()

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
add_library(X17_source_lib STATIC ${SOURCES} X17_dict.cxx)
target_link_libraries(X17_source_lib PUBLIC ${ROOT_LIBRARIES} Garfield::Garfield GSL::gsl ROOT::RIO)

ROOT_GENERATE_DICTIONARY(X17_dict Field.h Points.h MODULE X17_source_lib LINKDEF LinkDef.h)

# ---Build executables--------------------------------------------------------------------
# Create the build folder for the root project
set(PROJECT_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build)
file(MAKE_DIRECTORY ${PROJECT_BUILD_DIR})

# Create a list of subdirectories
set(SUBDIRECTORIES reconstruction simulations/ion_map simulations/single_track simulations/rk_tracks)

# Iterate over each subdirectory
foreach(subdirectory ${SUBDIRECTORIES})
    # Add the subdirectory
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${subdirectory})
endforeach()

# Specify the desired folders where you want to create the link to the ROOT dictionary
set(LINK_TARGET_FOLDER "${PROJECT_BUILD_DIR}/reconstruction")
set(LINK_TARGET_FOLDER2 "${PROJECT_BUILD_DIR}/simulations/rk_tracks")

# Create a symbolic link to the PCM file in the desired folder
execute_process(
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            "${PROJECT_BUILD_DIR}/libX17_source_lib_rdict.pcm"
            "${LINK_TARGET_FOLDER}/libX17_source_lib_rdict.pcm"
    RESULT_VARIABLE LINK_RESULT
)
execute_process(
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            "${PROJECT_BUILD_DIR}/libX17_source_lib_rdict.pcm"
            "${LINK_TARGET_FOLDER2}/libX17_source_lib_rdict.pcm"
    RESULT_VARIABLE LINK_RESULT2
)

# Check if the link creation was successful
if(LINK_RESULT EQUAL 0)
    message(STATUS "Symbolic link to libX17_source_lib_rdict.pcm created successfully for reconstruction.")
else()
    message(WARNING "Failed to create symbolic to libX17_source_lib_rdict.pcm link for reconstruction.")
endif()
if(LINK_RESULT2 EQUAL 0)
    message(STATUS "Symbolic link to libX17_source_lib_rdict.pcm created successfully for rk_tracks.")
else()
    message(WARNING "Failed to create symbolic to libX17_source_lib_rdict.pcm link for rk_tracks.")
endif()