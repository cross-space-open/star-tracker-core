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

#include "ImageProcess.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Define Bitmap constants */
#define DATA_OFFSET_OFFSET 0x000A
#define WIDTH_OFFSET 0x0012
#define HEIGHT_OFFSET 0x0016
#define BITS_PER_PIXEL_OFFSET 0x001C
#define HEADER_SIZE 14
#define INFO_HEADER_SIZE 40
#define NO_COMPRESION 0
#define MAX_NUMBER_OF_COLORS 0
#define ALL_COLORS_REQUIRED 0

/* Define image filtering terms */
#define FILTER_RADIUS 2 // median filter window radius in down-sampled image pixels
#define WINDOW_RADIUS 9 // centroiding window radius around a star's center pixel
#define MIN_NUMBER_OF_PIXELS 4

/* Define thresholding grid */
#define SUB_ROWS 4
#define SUB_COLS 6

/* Define star region sructure for sorting regions by brightness */
struct StarRegion {
	double x_c;
	double y_c;
	double sum;
};

/* Function used by qsort to compare the sum of two StarRegion elements */
int compare_region_sum(const void * a, const void * b) {
	return ( (*(struct StarRegion*)b).sum - (*(struct StarRegion*)a).sum );
}

/**
 * \brief Reads image bitmap and loads into live memory as a byte array.
 *
 * \param 	*filename		Name of the bitmap file to be loaded.
 * \param	**pixels		Pixel array to be filled from the raw image file.
 * \param	*width			Size of the bitmap image along the width.
 * \param	*height			Size of the bitmap image along the height.
 * \param	*bytesPerPixel	Size of each pixel (RGB or Grayscale).
 *
 * \return	void
 */
int read_image_bmp(const char *fileName,byte **pixels, int32 *width, int32 *height)
{
	byte L,*new_pixels;
	int i, j;
	int32 bytesPerPixel;

	FILE *imageFile = fopen(fileName, "rb");
	if (!imageFile){
		printf("Unable to open image file: %s.\n",fileName);
		return 0;
	};
	int32 dataOffset;
	fseek(imageFile, DATA_OFFSET_OFFSET, SEEK_SET);
	fread(&dataOffset, 4, 1, imageFile);
	fseek(imageFile, WIDTH_OFFSET, SEEK_SET);
	fread(width, 4, 1, imageFile);
	fseek(imageFile, HEIGHT_OFFSET, SEEK_SET);
	fread(height, 4, 1, imageFile);
	int16 bitsPerPixel;
	fseek(imageFile, BITS_PER_PIXEL_OFFSET, SEEK_SET);
	fread(&bitsPerPixel, 2, 1, imageFile);
	bytesPerPixel = ((int32)bitsPerPixel) / 8;

	int paddedRowSize = (int)(4 * ceil((float)(*width) / 4.0f))*bytesPerPixel;
	int unpaddedRowSize = (*width)*bytesPerPixel;
	int totalSize = unpaddedRowSize*(*height);
	*pixels = (byte*)malloc(totalSize);
	byte *currentRowPointer = *pixels+((*height-1)*unpaddedRowSize);
	for (i = 0; i < *height; i++)
	{
		fseek(imageFile, dataOffset+(i*paddedRowSize), SEEK_SET);
		fread(currentRowPointer, 1, unpaddedRowSize, imageFile);
		currentRowPointer -= unpaddedRowSize;
	}

	fclose(imageFile);

	/* Check if image is grayscale and convert if in RGB */
	if (bytesPerPixel != 1) {
		if (bytesPerPixel == 3) {
			new_pixels = (byte*) malloc(*height*(*width));
			for (i=0;i<*height;i++) {
				for (j=0;j<*width;j++) {
					/* Perform ITU-R 601-2 luma transform */
					L = (*pixels)[i*3*(*width)+j]*299/1000 + (*pixels)[i*3*(*width)+j+1]*587/1000 +
							(*pixels)[i*3*(*width)+j+2]*114/1000;
					new_pixels[i*(*width)+j] = L;
				}
			}
			free(*pixels);
			*pixels = (byte*) malloc(*width*(*height));
			memcpy(*pixels,new_pixels,*width*(*height));
			free(new_pixels);
		} else {
			printf("Incorrect image format with %d bytes per pixel.\n", bytesPerPixel);
			return 0;
		}
	}

	return 1;
}


/**
 * \brief 	Goes through the image and identifies regions of interest by generating
 * 			an array of potential star centre pixels exceeding a certain threshold
 * 			region size.
 *
 * \param 	*pixels				Image array.
 * \param	width				Size of the bitmap image along the width.
 * \param	height				Size of the bitmap image along the height.
 * \param	star_center_pixels	Array to be filled of the location of each star centre
 * 								pixel (x_1,y_1,...,x_N,y_N).
 *
 * \return	Number of identified regions of interest.
 */
int identify_star_centre_pixels(byte *pixels, int32 width, int32 height,
		int bright_thresh, int **star_center_pixels)
{
	double sum = 0.0, mean, std = 0.0;
	int i, j, k, l, rr, cc, *id_groups, length, n = 0, count = 0, N = 0, a, b, len = 0,
			*bright_pixels, *pixel_to_groups, *group_lengths, in_left_group = 0,
			in_up_group = 0, p, q, old_group_id, sub_width, sub_height;

	/* Calculate length of pixel sub grids */
	length = width*height/SUB_ROWS/SUB_COLS;
	sub_width = width/SUB_COLS;
	sub_height = height/SUB_ROWS;

	/* Find all pixels which are brighter than 10 sigma and not too close
	* to the image edge */
	bright_pixels = malloc(sizeof(int)*(width*height));

	/* Loop through each sub grid */
	for (rr = 0; rr < SUB_ROWS; rr++) {
		for (cc = 0; cc < SUB_COLS; cc++) {

				/* Calculate image standard deviation */
				std = 0;
				j = 0;
				for (i = 0; i < sub_height; i++) {
					for (j = j%sub_width; j < sub_width; j+=20) {
						sum += (double) pixels[(rr*sub_height+i)*width + (cc*sub_width+j)];
					}
				}

				mean = sum/(double) (length/20.);
				j = 0;
				for (i = 0; i < sub_height; i++) {
					for (j = j%sub_width; j < sub_width; j+=20) {
						std += pow((double) pixels[(rr*sub_height+i)*width + (cc*sub_width+j)]-mean, 2);
					}
				}
				std = sqrt(std/(double) (length/20.));

				for (i = 0; i < sub_height; i++) {
					for (j = 0; j < sub_width; j++) {
						if (pixels[(rr*sub_height+i)*width + (cc*sub_width+j)] > ((double) bright_thresh)*std) {
							bright_pixels[n*2] = rr*sub_height+i;
							bright_pixels[n*2+1] = cc*sub_width+j;
							n++;
						}
					}
				}			
		}
	}
	
	if (2*n > length) {
		free(bright_pixels);
		return 0;
	}


	/* Creates groups of bright pixels */
	pixel_to_groups = (int*) malloc(sizeof(int)*n*n);
	group_lengths = (int*) malloc(sizeof(int)*n);
	memset(group_lengths,0,sizeof(int)*n);
	id_groups = (int*) malloc(sizeof(int)*n);
	for (i = 0; i < n; i++) {
		/* Assign the left and upper pixel to the bright pixel */
		int left_pixel[2] = {bright_pixels[i*2],bright_pixels[i*2+1]-1};
		int up_pixel[2] = {bright_pixels[i*2]-1,bright_pixels[i*2+1]};

		/* Check if the left and up pixels have already been allocated to a group */
		for (j = 0; j < i; j++) {
			in_left_group = (bright_pixels[j*2] == left_pixel[0] &&
								bright_pixels[j*2+1] == left_pixel[1]);
			if (in_left_group) break;
		}
		for (k = 0; k < i; k++) {
			in_up_group = (bright_pixels[k*2] == up_pixel[0] &&
							bright_pixels[k*2+1] == up_pixel[1]);
			if (in_up_group) break;
		}

		/* If both are part of a group, add the current pixel and combine the groups */
		if (in_up_group && in_left_group && id_groups[j] != id_groups[k]) {
			/* Add the current pixel to the upper pixel's group */
			pixel_to_groups[id_groups[k]*n+group_lengths[id_groups[k]]] = i;
			group_lengths[id_groups[k]]++;
			/* Append the upper pixel group onto the left pixel group */
			for (l = 0; l < group_lengths[id_groups[k]]; l++) {
				pixel_to_groups[id_groups[j]*n+group_lengths[id_groups[j]]] =
						pixel_to_groups[id_groups[k]*n+l];
				group_lengths[id_groups[j]]++;
			}
			/* Replace all of the upper pixel group's entries
			 * with references to the left pixel group */
			old_group_id = id_groups[k];
			for (l = 0; l < group_lengths[id_groups[k]]; l++) {
				id_groups[pixel_to_groups[id_groups[k]*n+l]] =
						id_groups[j];
			}
			group_lengths[old_group_id] = 0;

	    /* If exactly one of the left or upper pixels is part of an existing group,
	     * add the current pixel to that group */
		} else if (in_left_group) {
			/* Add the current pixel to the left pixel's group */
			pixel_to_groups[id_groups[j]*n+group_lengths[id_groups[j]]] = i;
			group_lengths[id_groups[j]]++;
			/* Add the current pixel with the left pixel group entries */
			id_groups[i] = id_groups[j];
		} else if (in_up_group) {
			/* Add the current pixel to the upper pixel's group */
			pixel_to_groups[id_groups[k]*n+group_lengths[id_groups[k]]] = i;
			group_lengths[id_groups[k]]++;
			/* Add the current pixel with the upper pixel group entries */
			id_groups[i] = id_groups[k];
		/* Otherrwise, add the current pixel as a single entry */
		} else {
			pixel_to_groups[count*n] = i;
			id_groups[i] = count;
			group_lengths[count]++;
			count++;
		}
	}

	/* Identify the brightest pixel of each, making sure that each is more than the
	 * minimum number of pixels */
	for (i = 0; i < count; i++) {
		if (group_lengths[i] >= MIN_NUMBER_OF_PIXELS)
			N++;
	}
	if (N > 0) {
		*star_center_pixels = (int*) malloc(sizeof(int)*N*2);
		for (i = 0; i < count; i++) {
			if (group_lengths[i] >= MIN_NUMBER_OF_PIXELS) {
				p = pixel_to_groups[i*n];
				a = pixels[bright_pixels[p*2]*width + bright_pixels[p*2+1]];
				for (j = 0; j < group_lengths[i]; j++) {
					q = pixel_to_groups[i*n+j];
					b = pixels[bright_pixels[q*2]*width + bright_pixels[q*2+1]];
					if (b > a) {
						p = q;
						a = b;
					}
				}
				(*star_center_pixels)[len*2] = bright_pixels[p*2];
				(*star_center_pixels)[len*2+1] = bright_pixels[p*2+1];
				len++;
			}
		}
	}
	
	free(bright_pixels);
	free(pixel_to_groups);
	free(group_lengths);
	free(id_groups);
	return N;
}

/**
 * \brief 	Filters a local pixel bytes using a Gaussian.
 *
 * \param	y				Location of pixel in bitmap by vertical.
 * \param	x				Location of pixel in bitmap by horizontal.
 * \param 	*pixels			Image array.
 * \param	*width			Size of the bitmap image along the width.
 *
 * \return	void
 */
double local_filter(int y, int x, byte *pixels, int32 width) {

	int i, j;
	double gaussian_pixel = 0.0, min;

	/* Gaussian Kernal with for 5x5 filter window and sigma=2.0 */
	double K[25] = {0.0146374578810798,	0.0212973755488066,	0.0241330881575135,	0.0212973755488066,	0.0146374578810798,
				0.0212973755488066,	0.0309874985774132,	0.0351134360774063,	0.0309874985774132,	0.0212973755488066,
				0.0241330881575135,	0.0351134360774063,	0.0397887357729738,	0.0351134360774063,	0.0241330881575135,
				0.0212973755488066,	0.0309874985774132,	0.0351134360774063,	0.0309874985774132,	0.0212973755488066,
				0.0146374578810798,	0.0212973755488066,	0.0241330881575135,	0.0212973755488066,	0.0146374578810798};

	/* Apply Python-based 2D Gaussian filter */
	for (i=-2;i<3;i++) {
		for (j=-2;j<3;j++) {
			gaussian_pixel += K[(i+2)*5+j+2]*(double)pixels[(y+i)*width+x+j];
		}
	}

	/* Normalise image by subtract the minimum of the image pixel and the local Gaussian so > 0 */
	min = fmin(gaussian_pixel,(double)pixels[i*width+j]);
	gaussian_pixel -= min;

	return gaussian_pixel;
}

/**
 * \brief Calculate the centroid of stars according to Moment method.
 *
 * \param	*star_center_pixels	Location of each star centre pixel
 * 								(x_1,y_1,...,x_N,y_N).
 * \param	N					Number of identified regions of interest.
 * \param 	*pixels				Image array.
 * \param	width				Size of the bitmap image along the width.
 * \param	height				Size of the bitmap image along the height.
 * \param 	**star_centroids	Array to be filled with location of each star
 * 								centroid (x_c1,y_c1,...,x_cN,y_cN).
 *
 * \return	Number of calculated centroids.
 */
int moment(int *star_center_pixels, int N,  byte *pixels, int32 height,
		int32 width, double **star_centroids) {

	int i,j,k,x,y,num=0;
	double x_c,y_c,sumX,sumY,sum,V;

	*star_centroids = (double*) malloc(sizeof(double)*N*2);
	struct StarRegion star_regions[N];

	for (i = 0; i < N; i++) {

		/* Throw out star if it's too close to the edge of the image */
		y = star_center_pixels[i*2];
		x = star_center_pixels[i*2+1];
		if (y < WINDOW_RADIUS + FILTER_RADIUS || y >= (height - WINDOW_RADIUS - FILTER_RADIUS) ||
				x < WINDOW_RADIUS + FILTER_RADIUS || x >= (width  - WINDOW_RADIUS - FILTER_RADIUS))
			continue;

		/* Set initial centroids at brightest pixel position */
		x_c = (double) x;
		y_c = (double) y;

    	/* Apply the 7x7 grid approach to find horizontal centroid */
		sumX = 0.0; sumY = 0.0; sum = 0.0;
		for (j = -WINDOW_RADIUS; j <= WINDOW_RADIUS; j++) {
			for (k = -WINDOW_RADIUS; k <= WINDOW_RADIUS; k++) {
				V = local_filter(y+j, x+k, pixels, width);
				sumX += (double)k*V;
				sumY += (double)j*V;
				sum += V;

			}
		}
		x_c += sumX/sum + 0.5;
		y_c += sumY/sum + 0.5;

		star_regions[num].x_c = x_c;
		star_regions[num].y_c = y_c;
		star_regions[num].sum = sum;

		// (*star_centroids)[num*2] = y_c;
		// (*star_centroids)[num*2+1] = x_c;
		
		num++;
	}

	/* Store the star centroids in order of brightness */
	qsort(star_regions, N, sizeof(struct StarRegion), compare_region_sum);
	for (i = 0; i < num; i++) {
		(*star_centroids)[i*2] = star_regions[i].y_c;
		(*star_centroids)[i*2+1] = star_regions[i].x_c;
	}

	return num;
}
