# ########################################################################################
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)

# We want to export the compile commands for clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# So that repeated builds don't silence warnings
set(CMAKE_CXX_COMPILER_LAUNCHER ccache)

# Set the C++ standard (C++23 not supported in old versions of CMake)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(ROOT_CXX_FLAGS "-std=c++20")

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
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

# clang-tidy for better code quality
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    # Find the clang-tidy executable
    find_program(CLANG_TIDY_EXE NAMES "clang-tidy" "clang-tidy-18" "clang-tidy-19" "clang-tidy-23")
    if(CLANG_TIDY_EXE)
        message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
        # Tell CMake to run clang-tidy on every C++ file during the build
        set(CMAKE_CXX_CLANG_TIDY
            "${CLANG_TIDY_EXE}"
            "--header-filter=${CMAKE_SOURCE_DIR}/.*"
            "--config-file=${CMAKE_SOURCE_DIR}/.clang-tidy"
            "--format-style='file'"
            # "--dump-config"
            "--allow-no-checks")
        set_source_files_properties(
            X17_dict.cxx
            PROPERTIES SKIP_LINTING TRUE
        )
    else()
        message(WARNING "clang-tidy not found!")
    endif()
endif()

# ---Garfield-----------------------------------------------------------------------------
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED HINTS "/home/vavrik/work/source/garfield/build/CMake" "/storage/projects/utefx17/SourceCode/garfield_2603/build/CMake")
endif()

# ---Root---------------------------------------------------------------------------------
# Check if ROOTSYS is defined in the environment
if(NOT DEFINED ENV{ROOTSYS})
    message(FATAL_ERROR "
    ROOT environment not found! 
    For VS Code: it cannot see your ~/.bashrc when launched from the GUI.")
endif()

message(STATUS "Looking for ROOT")
find_package(ROOT REQUIRED)
message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}")

set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")

find_package(GSL REQUIRED)
message(STATUS "CMAKE_CXX_STANDARD: ${CMAKE_CXX_STANDARD}")

# ---Add subdirectory containing source code----------------------------------------------
# include_directories(include)
file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp")
add_library(X17_source SHARED ${SOURCES})

target_include_directories(X17_source PUBLIC include)
target_include_directories(X17_source SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(X17_source 
    PUBLIC 
    ${ROOT_LIBRARIES} 
    Garfield::Garfield 
    GSL::gsl 
    ROOT::RIO
)

ROOT_GENERATE_DICTIONARY(X17_dict Field.h Matrix.h Points.h Track.h Vector.h MODULE X17_source LINKDEF LinkDef.h)

# ---Build executables--------------------------------------------------------------------
# Create the build folder for the root project
set(PROJECT_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build)
file(MAKE_DIRECTORY ${PROJECT_BUILD_DIR})

# Create a list of subdirectories
set(SUBDIRECTORIES
    reconstruction
    simulations/ion_map
    simulations/single_track
    simulations/rk_tracks
    simulations/micro_tracks
    tests/magnetic_field_plot
    tests/thesis_plots
    tests/mardia_test
    tests/driftlines
    tests/ionoel_pad_dist
)

# Iterate over each subdirectory
foreach(subdirectory ${SUBDIRECTORIES})
    # Add the subdirectory
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${subdirectory})
endforeach()

# List of target folders for ROOT dictionary link
set(LINK_TARGET_FOLDERS
    "${PROJECT_BUILD_DIR}/reconstruction"
    "${PROJECT_BUILD_DIR}/simulations/rk_tracks"
    "${PROJECT_BUILD_DIR}/simulations/single_track"
    "${PROJECT_BUILD_DIR}/simulations/ion_map"
    "${PROJECT_BUILD_DIR}/simulations/micro_tracks"
    "${PROJECT_BUILD_DIR}/tests/magnetic_field_plot"
    "${PROJECT_BUILD_DIR}/tests/thesis_plots"
    "${PROJECT_BUILD_DIR}/tests/driftlines"
    "${PROJECT_BUILD_DIR}/tests/ionoel_pad_dist"
)

# Path to the source PCM file
# set(PCM_SOURCE "${PROJECT_BUILD_DIR}/libX17_source_rdict.pcm")

# # Loop through each target folder
# foreach(target_folder ${LINK_TARGET_FOLDERS})
#     execute_process(
#         COMMAND ${CMAKE_COMMAND} -E create_symlink "${PCM_SOURCE}" "${target_folder}/libX17_source_rdict.pcm"
#         RESULT_VARIABLE LINK_RESULT
#     )
    
#     if(LINK_RESULT EQUAL 0)
#         message(STATUS "Symbolic link to libX17_source_rdict.pcm created successfully for ${target_folder}.")
#     else()
#         message(WARNING "Failed to create symbolic link to libX17_source_rdict.pcm for ${target_folder}.")
#     endif()
# endforeach()