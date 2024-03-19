#Set the output directory of the library
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

#Find the external libraries installed in the system (MPI, OpenMP)

# Set a definition to include the mpi header
add_definitions(-DUSE_MPI)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(OpenMP REQUIRED)
#include_directories(${OpenMP_INCLUDE_PATH})

#Include the directories to look for header files
set(INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)

#Include the directories
include_directories(${INCLUDE_DIRS})

#Create the library
add_library(mandelbrot_hybrid SHARED mandelbrot.c)

#Link the library with the external libraries
target_link_libraries(mandelbrot_hybrid PRIVATE MPI::MPI_C)

set_target_properties(mandelbrot_hybrid PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(mandelbrot_hybrid PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")

#Set the version of the library
set_target_properties(mandelbrot_hybrid PROPERTIES VERSION 1.0.0 SOVERSION 1)
