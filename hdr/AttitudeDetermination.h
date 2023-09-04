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

void rot2quat(double **A, double *q);
void calculate_attitude_svd(double **A, double **image_vectors, double **catalog_vectors, int N);
void calculate_attitude_q_method(double *q, double **image_vectors, double **catalog_vectors, int N);
void img2cat(double **A, double *image_vector, double *catalog_vector);
void cat2img(double **A, double *catalog_vector, double *image_vector);
void ekf_time_update(double *x, double **P, double **Q, double dt);
void ekf_state_time_update(double *x, double **P, double **Q, double dt);
void ekf_measurement_update(double *x, double **P, double sigma_meas,
              double **image_vectors, double **catalog_vectors, int N);
void ekf_state_measurement_update(double *x, double **P, double sigma_meas,
              double **image_vectors, double **catalog_vectors, int N);
void quat2rot(double *q, double **A);
void estimate_angular_velocity(double **A1, double **A2, double dt, double *omega);
