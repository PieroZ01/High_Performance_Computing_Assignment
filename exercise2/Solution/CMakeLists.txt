cmake_minimum_required(VERSION 3.28.1)
project(Mandelbrot_hybrid_MPI_OpenMP
        VERSION 1.0
        DESCRIPTION "HPC exam - Exercise2C"
        LANGUAGES C)

# C Standard to use
set(CMAKE_C_STANDARD 11)

#C Flags
set(CMAKE_C_FLAGS "-O3 -march=native -Wall -Wextra -Wpedantic -Werror -lm -fopenmp")

# Define the local variables
set(INCLUDE_DIR "include/")
set(SRC_DIR "src/")
set(APPS_DIR "apps/")
set(INCLUDE_FILE "mandelbrot.h")
set(SRC_FILE "mandelbrot.c")
set(APPS_FILE "main.c")
set(EXE_FILE "main")

# Find the required external libraries
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
find_package(OpenMP REQUIRED)
include_directories(${OpenMP_INCLUDE_PATH})

#Set the output path of the executables: bin/
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# Add the executable
add_executable(${EXE_FILE} ${SRC_DIR}${SRC_FILE} ${APPS_DIR}${APPS_FILE})

# Set the compile options for openMP
target_compile_options(${EXE_FILE} PRIVATE ${OpenMP_C_FLAGS})

# Include the directories
target_include_directories(${EXE_FILE} PUBLIC ${INCLUDE_DIR} ${MPI_INCLUDE_PATH} ${OpenMP_INCLUDE_PATH})

# Link the executable with the MPI and OpenMP libraries
target_link_libraries(${EXE_FILE} ${MPI_LIBRARIES} ${OpenMP_LIBRARIES} gomp)
