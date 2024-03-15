// Libraries
#include <complex.h>

// Constants
#define XWIDTH 10000
#define YWIDTH 10000
#define MAXVAL 65535

// Function declarations
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);

double complex f_c(double complex z, double complex c);

int mandelbrot(double complex c, int max_iter);
