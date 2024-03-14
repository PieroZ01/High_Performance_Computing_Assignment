// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <complex.h>

// Constants
#define XWIDTH 256
#define YWIDTH 256
#define MAXVAL 65535

// Function to write a PGM file (grey scale image)
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name)
{
  FILE* image_file; 
  image_file = fopen(image_name, "w"); 
  
  // Writing header
  // The header's format is as follows, all in ASCII.
  // "whitespace" is either a blank or a TAB or a CF or a LF
  // - The Magic Number (see below the magic numbers)
  // - the image's width
  // - the height
  // - a white space
  // - the image's height
  // - a whitespace
  // - the maximum color value, which must be between 0 and 65535
  //
  // if he maximum color value is in the range [0-255], then
  // a pixel will be expressed by a single byte; if the maximum is
  // larger than 255, then 2 bytes will be needed for each pixel
  //

  int color_depth = 1 + ( maxval > 255 );
  fprintf(image_file, "P5\n%d %d\n%d\n", xsize, ysize, maxval);
  
  // Writing file
  fwrite(image, 1, xsize*ysize*color_depth, image_file);  

  fclose(image_file);
  return ;
}

// Define the complex 
