#Set the output path of the executables: bin/
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

#Include the header file for the executables
set(INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)

# Set a definition to include the mpi header
add_definitions(-DUSE_MPI)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(OpenMP REQUIRED)

#Include the directories that will be linked to the executables
include_directories(${INCLUDE_DIRS})

#Create the executables
add_executable(main main.c)

#Link the library to the executables
target_link_libraries(main mandelbrot_hybrid)
target_link_libraries(main MPI::MPI_C)

set_target_properties(main PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(main PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
