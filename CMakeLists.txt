cmake_minimum_required(VERSION 3.10)
project(euler)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# define sources
set(SOURCES main.cpp utils/ equations/ )

# Create executable
add_executable(fvm ${SOURCES})

