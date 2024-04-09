# High Performance Computing final project - Exercise 2C

## Student's Info

| Name | Surname | Student ID | UniTs email | Personal email | Master course |
|:---:|:---:|:---:|:---:|:---:|:---:|
| Piero | Zappi | SM3600004 | `piero.zappi@studenti.units.it` | `piero.z.2001@gmail.com` | SDIC |

## Folder Organization

This folder is organized with the following structure:

```bash
.
├── CMakeLists.txt
├── Image # Mandelbrot set's image
│   ├── mandelbrot_5000x5000x65535.pgm
│   └── mandelbrot_5000x5000x65535.png
├── Report_and_slides # Final report and slides for the presentation
│   ├── slides_ex2C
│   └── Zappi_ex2C_report.pdf
├── README.MD # This file
├── Results # Data collected on the ORFEO cluster (csv files)
│   ├── mpi_scaling_mandelbrot_1000x1000x65535_core_norobin.csv
│   ├── mpi_strongscaling_1000x1000x65535_core_robin.csv
│   ├── mpi_strongscaling_1000x1000x65535_socket_robin.csv
│   ├── mpi_weakscaling_sizex1000x65535_core_robin.csv
│   ├── omp_scaling_mandelbrot_1000x1000x65535_close_cores_norobin.csv
│   ├── omp_strongscaling_1000x1000x65535_close_cores_dynamic_robin.csv
│   ├── omp_strongscaling_1000x1000x65535_close_cores_static_robin.csv
│   ├── omp_strongscaling_1000x1000x65535_spread_cores_dynamic_robin.csv
│   ├── omp_weakscaling_sizex1000x65535_close_cores_dynamic_robin.csv
│   ├── serial_1000x1000x65535_robin.csv
│   └── serial_100x1000x65535_robin.csv
├── apps # Main program
│   └── main.c
├── build.sh # Bash script for quick and automatic compilation
├── ex2C_scaling_analysis.ipynb # Notebook for the analysis of the results
├── include # Header file
│   └── mandelbrot.h
│ # Some bash scripts to automate the data collection process
├── mpi_strongscaling.sh
├── mpi_weakscaling.sh
├── omp_strongscaling.sh
├── omp_weakscaling.sh
├── serial.sh
└── src # Source file
    └── mandelbrot.c
```

## Contents

This folder contains the solution for exercise 2C, which requested to implement an hybrid MPI+OpenMP code to compute the Mandelbrot set. Further details can be found in the **Zappi_ex2C_report.pdf** file.

## Compilation

From the root folder of this project it is possible to compile it with the `cmake` build system.\
By executing the following commands (and eventually specifying the correct `-DCCMAKE_PREFIX_PATH` for eventual third party libraries if not found) all the `C` codes will be compiled in the `build/` folder.

```bash
cmake -S . -B build/

make -C build/ -j<N>
```
The **build.sh** bash script is provided for quick and automatic compilation. To build on the ORFEO cluster, simply run the following command:

```bash
sbatch build.sh
```
