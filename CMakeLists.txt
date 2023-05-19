# CMake 3.12.4 required for 20 to be a valid value for CXX_STANDARD
cmake_minimum_required(VERSION 3.12.4)

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.15)
    # Add MSVC runtime library selection flags automatically:
    cmake_policy(SET CMP0091 OLD)

    # Don't override the warning flags in MSVC:
    cmake_policy(SET CMP0092 NEW)
endif()

if(POLICY CMP0076)
    cmake_policy(SET CMP0076 NEW)
endif()

project(akita-inu)
include(CMakeDependentOption)

option(WARNINGS_AS_ERRORS "Treat warnings as errors" ON)
CMAKE_DEPENDENT_OPTION(COMPILE_WITH_DWARF "Add DWARF debugging information" ON "MINGW" OFF)

# Sanity check : Check that all submodules are present
# =======================================================================
function(check_submodules_present)
    file(READ "${PROJECT_SOURCE_DIR}/.gitmodules" gitmodules)
    string(REGEX MATCHALL "path *= *[^ \t\r\n]*" gitmodules ${gitmodules})

    foreach(module ${gitmodules})
        string(REGEX REPLACE "path *= *" "" module ${module})

        if(NOT EXISTS "${PROJECT_SOURCE_DIR}/${module}/.git")
            message(SEND_ERROR "Git submodule ${module} not found."
                "Please run: git submodule update --init --recursive")
        endif()
    endforeach()
endfunction()

check_submodules_present()

# Configure C++ standard
# ===========================
if(MSVC)
    add_compile_options(/std:c++latest)
    add_definitions(-D_HAS_DEPRECATED_RESULT_OF)
else()
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

# set up output paths for executable binaries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# System imported libraries
# ======================

# TODO: Is this necessary?

# Prefer the -pthread flag on Linux.
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Include source code
# ===================
add_subdirectory(externals)
add_subdirectory(src)
