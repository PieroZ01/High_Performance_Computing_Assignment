// Libraries
#include <complex.h>

// Constants
#define XWIDTH 12288
#define YWIDTH 12288
#define MAXVAL 65535

// Function declarations
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);

double complex f_c(double complex z, double complex c);

int mandelbrot(double complex c, int max_iter);