/* Copyright (c) 2023 CROSS Space

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. */

#include <stdio.h>

typedef unsigned int int32;
typedef short int16;
typedef unsigned char byte;

int read_image_bmp(const char *fileName,byte **pixels, int32 *width, int32 *height);
int identify_star_centre_pixels(byte *pixels, int32 width, int32 height,
		int bright_thresh, int **star_center_pixels);
int gaussian_grid(int *star_center_pixels, int N,  byte *pixels, int32 height,
		int32 width, double **star_centroids);
int moment(int *star_center_pixels, int N,  byte *pixels, int32 height,
		int32 width, double **star_centroids);
