#ifndef MANDELBROT_H
#define MANDELBROT_H

// Libraries
#include <complex.h>
#include <time.h>

// Constants
#define XWIDTH 12288
#define YWIDTH 12288
#define MAXVAL 65535

// measure the wall-clock time
#define CPU_TIME (clock_gettime( CLOCK_REALTIME, &ts ), (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9)

// Function declarations

// Function to write a PGM file (grey scale image)
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);

// Declare the function f_c(z)=z^2+c
double complex f_c(const double complex z, const double complex c);

// Declare the function to compute the mandelbrot set
int mandelbrot(const double complex c, const int max_iter);

// Declare function to compute the mean value
double mean(const double *, const int n);

// Declare function to compute the standard deviation
double std_dev(const double *, const int n);

#endif // MANDELBROT_H
