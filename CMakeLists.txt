cmake_minimum_required(VERSION 3.10)
project(euelr_system)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# define sources
set(SOURCES main.cpp utils.cpp utils.h)
# Create executable
add_executable(advection ${SOURCES})

