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

void calibrate_nlls_fov(double *fov, int height, int width, double ***image_vectors,
		double **catalog_vectors, double *k_radial, double *offset, int n);
void compute_vectors(double fov, int height, int width, double* star_centroids, int n,
		double* offset, double ***image_vectors) ;
void compute_calibrated_vectors(double fov, int height, int width, double* star_centroids, int n,
		double* k_radial, double* offset, double ***image_vectors);
void calibrate_nlls_focal(double *fov, int height, int width, double ***image_vectors,
		double **catalog_vectors, double *k_radial, double **offset, int n);
void calibrate_lls_radial(double fov, int height, int width, double **image_vectors,
		double* q, double **catalog_vectors, double **k_radial, double *offset, int n);
void calibrate_lm_orion(double *fov, int height, int width, double **image_vectors,
		double *q, double **catalog_vectors, double **k_radial, double *offset, int n);
