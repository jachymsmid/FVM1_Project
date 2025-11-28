cmake_minimum_required(VERSION 3.10)
project(euler)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# define sources
set(SOURCES main.cpp utils/ equations/ )

# Create executable
add_executable(fvm ${SOURCES})

