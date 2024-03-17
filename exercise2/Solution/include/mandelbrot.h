#ifndef MANDELBROT_H
#define MANDELBROT_H

// Libraries
#include <complex.h>

// Constants
#define XWIDTH 12288
#define YWIDTH 12288
#define MAXVAL 65535

// Function declarations
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);

inline double complex f_c(const double complex z, const double complex c);

const int mandelbrot(const double complex c, const int max_iter);

#endif // MANDELBROT_H
