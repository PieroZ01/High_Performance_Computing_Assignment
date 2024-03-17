#ifndef MANDELBROT_H
#define MANDELBROT_H

// Libraries
#include <complex.h>

// Constants
#define XWIDTH 12288
#define YWIDTH 12288
#define MAXVAL 65535

// Function declarations

// Function to write a PGM file (grey scale image)
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);

// Declare the function to compute the mandelbrot set
const int mandelbrot(const double complex c, const int max_iter);

#endif // MANDELBROT_H
