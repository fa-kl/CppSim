# ---------------------------------------------------------------------------------------
#
# @file: CMakeLists.txt
#
# @brief: The simulator's main CMakeLists.txt file
#
# @author: fakl
# @date: December 2025
#
# ---------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------
# Setup & Settings
# ---------------------------------------------------------------------------------------
set(LIB_NAME "MathLib")
set(LIB_VERSION "1.0.0")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


cmake_minimum_required(VERSION 3.20)
project(${LIB_NAME} VERSION ${LIB_VERSION} LANGUAGES CXX)

add_compile_options(
    -Wall                   # Enable common warnings
    -Wextra                 # Enable extra warnings
    -Wpedantic              # Enable strict compliance warnings
    #-Werror                 # Treat warnings as errors
    -Wshadow                # Warn about variable shadowing
    -Wformat=2              # Check for format string issues (use `-Wformat=2` for stricter)
    -Wconversion            # Warn on implicit type conversions
    -Wsign-conversion       # Warn about signed/unsigned comparisons
    -Wuninitialized         # Warn if a variable is used before initialization
    -Wno-unused-parameter   # Optional: disable unused parameter warnings if needed
)

# ---------------------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------------------

# Find SDL2
find_package(PkgConfig REQUIRED)
pkg_check_modules(SDL2 REQUIRED sdl2)

# ---------------------------------------------------------------------------------------
# Sources & Headers
# ---------------------------------------------------------------------------------------
set(MAIN_TARGET Main)
set(INCLUDE_DIR inc/)

# Collect all .cpp sources in src/
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_executable(${MAIN_TARGET} main.cpp ${SOURCES})
target_include_directories(${MAIN_TARGET} PRIVATE ${INCLUDE_DIR})

# Link SDL2
target_link_libraries(${MAIN_TARGET} ${SDL2_LIBRARIES})
target_include_directories(${MAIN_TARGET} PRIVATE ${SDL2_INCLUDE_DIRS})
target_compile_options(${MAIN_TARGET} PRIVATE ${SDL2_CFLAGS_OTHER})
