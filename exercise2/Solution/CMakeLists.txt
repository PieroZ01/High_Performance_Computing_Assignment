project(Mandelbrot_hybrid_MPI_OpenMP
        VERSION 1.0
        DESCRIPTION "HPC exam - Exercise2C"
        LANGUAGES CXX)
cmake_minimum_required(VERSION 3.12)

#C Flags
set(CMAKE_C_FLAGS "-O3 -march=native")

#Add subdirectories with CMakeLists.txt files to compile:
add_subdirectory(src)
add_subdirectory(apps)
