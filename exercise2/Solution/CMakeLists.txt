cmake_minimum_required(VERSION 3.28.1)
project(Mandelbrot_hybrid_MPI_OpenMP
        VERSION 1.0
        DESCRIPTION "HPC exam - Exercise2C"
        LANGUAGES CXX)

# C Standard to use
set(CMAKE_C_STANDARD 11)


#C Flags
set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native -Wall -Wextra -Wpedantic")
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

#Add subdirectories with CMakeLists.txt files to compile:
add_subdirectory(src)
add_subdirectory(apps)