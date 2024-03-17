#Set the output path of the executables: bin/
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

#Include the header file for the executables
set(INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)

#Include the directories that will be linked to the executables
include_directories(${INCLUDE_DIRS})

#Create the executables
add_executable(main main.c)

#Link the library to the executables
target_link_libraries(main mandelbrot_hybrid)
target_link_libraries(main MPI::MPI_C)
