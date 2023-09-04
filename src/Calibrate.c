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

#include "Calibrate.h"
#include "matrix.h"
#include "AttitudeDetermination.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.1415926
#define TOL_F 1E-6
#define TOL_OFFSET 1E-6
#define TOL_LM_ORION 1E-10
#define MAX_NLLS 50

/**
 * \brief 	Recursive function that calculates the binomial coefficient
 *
 * \param 	n	Number of values to consider.
 * \param	k	Number in a single combination.
 *
 * \return	Binomial coefficient
 */
int nchoosek(int n, int k) {

	double r;

	if (k == 0)
		r = 1.0;
	else
		r = (double) n / (double) k * (double) nchoosek(n-1,k-1);

	return (int) round(r);
}

/**
 * \brief 	Calculates the image unit vectors using the measured centroids.
 *
 * \param 	fov 				Image field of view [deg]
 * \param	height				Size of the bitmap image along the height.
 * \param	width				Size of the bitmap image along the width.
 * \param	n					Total number of determined star centroids.
 * \param	*offset				Image offset from center [px]. (dx,dy)
 * \param	***image_vectors	Set of identified vectors to be calculated in the image frame.
 *
 * \return	void
 */
void compute_vectors(double fov, int height, int width, double* star_centroids, int n,
		double* offset, double ***image_vectors) {

	int i,u_0,v_0;
	double s,u,v;

	/* Calculate center pixel */
	u_0 = width/2.0 + offset[0];
	v_0 = height/2.0 + offset[1];

	/* Calculate scale factor */
	fov = fov*PI/180;
	s = tan(fov/2.0)/(double)(width/2.0);

	/* Calculate star vectors (i,j,k)*/
	for (i = 0; i < n; i++) {
		v = star_centroids[i*2];
		u = star_centroids[i*2+1];
		(*image_vectors)[i][0] = 1/sqrt(1+pow(s*(u_0-u),2)+pow(s*(v_0-v),2));
		(*image_vectors)[i][1] = s*(u_0-u)*(*image_vectors)[i][0];
		(*image_vectors)[i][2] = s*(v_0-v)*(*image_vectors)[i][0];
	}
}

/**
 * \brief 	Calculates the image unit vectors with radial distortion corrections using the measured centroids.
 *
 * \param 	fov 				Image field of view [deg]
 * \param	height				Size of the bitmap image along the height.
 * \param	width				Size of the bitmap image along the width.
 * \param	*star_centroids		Star centroid positions.
 * \param	n					Total number of determined star centroids.
 * \param	*k_radial			Radial distortion parameters. (k_1,k_2)
 * \param	*offset				Image offset from center [px]. (dx,dy)
 * \param	***image_vectors	Set of identified vectors to be calculated in the image frame.
 *
 * \return	void
 */
void compute_calibrated_vectors(double fov, int height, int width, double* star_centroids, int n,
		double* k_radial, double* offset, double ***image_vectors) {

	int iter,i,u_0,v_0;
	double s,x,y,r,x_d,y_d,det;
	double J[4],fx,fy;

	/* Set tolerance */
	double tol = 1e-30;
	int max_iter = 50;

	/* Calculate center pixel */
	u_0 = width/2.0 + offset[0];
	v_0 = height/2.0 + offset[1];

	/* Calculate scale factor */
	fov = fov*PI/180;
	s = tan(fov/2.0)/(double)(width/2.0);

	/* Calculate star vectors (i,j,k)*/
	for (i = 0; i < n; i++) {
		y = v_0-star_centroids[i*2];
		x = u_0-star_centroids[i*2+1];
		y_d = v_0-star_centroids[i*2];
		x_d = u_0-star_centroids[i*2+1];

		/* Undistort star position in radial */
		r = x*x + y*y;
		fx = (1 + k_radial[0]*r + k_radial[1]*r*r)*x - x_d;
		fy = (1 + k_radial[0]*r + k_radial[1]*r*r)*y - y_d;
		iter = 0;
		while ((fx*fx > tol || fy*fy > tol) && iter < max_iter) {
			J[0] = (1 + k_radial[0]*r + k_radial[1]*r*r) + x*x*(2*k_radial[0] + 4*r*k_radial[1]);
			J[1] = x*y*(2*k_radial[0] + 4*r*k_radial[1]);
			J[2] = x*y*(2*k_radial[0] + 4*r*k_radial[1]);
			J[3] = (1 + k_radial[0]*r + k_radial[1]*r*r) + y*y*(2*k_radial[0] + 4*r*k_radial[1]);
			det = J[0]*J[3] - J[1]*J[2];
			x = x - J[3]/det*fx + J[1]/det*fy;
			y = y + J[2]/det*fx - J[0]/det*fy;
			r = x*x + y*y;
			fx = (1 + k_radial[0]*r + k_radial[1]*r*r)*x - x_d;
			fy = (1 + k_radial[0]*r + k_radial[1]*r*r)*y - y_d;
			iter++;
		}

		/* Calculate unit vector */
		(*image_vectors)[i][0] = 1/sqrt(1+pow(s*x,2)+pow(s*y,2));
		(*image_vectors)[i][1] = s*x*(*image_vectors)[i][0];
		(*image_vectors)[i][2] = s*y*(*image_vectors)[i][0];
	}

}

/**
 * \brief 	Estimates the image radial distortion parameters using a basic LLS
 *
 * \param 	fov 				Estimate image field of view [deg]
 * \param	height				Size of the bitmap image along the height.
 * \param	width				Size of the bitmap image along the width.
 * \param	**image_vectors		Set of measured star vectors in the image frame.
 * \param	**catalog_vectors	Set of catalog vectors in the celestial frame.
 * \param	*q					Quaternion vector set (q0,q1,q2,q3).
 * \param	**k_radial			Radial distortion parameters. (k_1,k_2)
 * \param	*offset				Image offset from center [px]. (dx,dy)
 * \param	n					Total number of stars to be used in the estimation.
 *
 * \return	void
 */
void calibrate_lls_radial(double fov, int height, int width, double **image_vectors,
		double* q, double **catalog_vectors, double **k_radial, double *offset, int n) {

	int i;
	double xhat_i, yhat_i, x_i, y_i, x_0, y_0, f, r, *star_centroids;
	mat X, deltaz, beta, A, true_vector, catalog_vector;
	mat XT, XTX, invXTX, invXTXXT;

	/* Assign vectors */
	X = mat_create(n*2,2);
	deltaz = mat_create(n*2,1);
	beta = mat_create(2,1);
	XT = mat_create(2,n);
	XTX = mat_create(2,2);
	invXTX = mat_create(2,2);
	invXTXXT = mat_create(2,n);
	A = mat_create(3,3);
	true_vector = mat_create(3,1);
	catalog_vector = mat_create(3,1);
	star_centroids = malloc(2*n*sizeof(double));


	/* Calculate focal parameters */
	x_0 = (double) width/2.0 + offset[0];
	y_0 = (double) height/2.0 + offset[1];
	f = width/2.0/tan(fov*PI/180.0/2.0);

 	/* Retrieve measured or distortion calibrated coordinates */
	for (i = 0; i < n; i++) {
		star_centroids[i*2] = y_0 - image_vectors[i][2]/image_vectors[i][0]*f;
		star_centroids[i*2+1] = x_0 - image_vectors[i][1]/image_vectors[i][0]*f;
	}

	/* Set attitude matrix */
	quat2rot(q,A);

	for (i = 0; i < n; i++) {

		/* Set measured coordinates */
		yhat_i = star_centroids[i*2];
		xhat_i = star_centroids[i*2+1];

		/* Calculate measurement matrix */
		r = (x_0-xhat_i)*(x_0-xhat_i) + (y_0-yhat_i)*(y_0-yhat_i);
		X[i*2][0] = r*(x_0-xhat_i);
		X[i*2][1] = r*r*(x_0-xhat_i);
		X[i*2+1][0] = r*(y_0-yhat_i);
		X[i*2+1][1] = r*r*(y_0-yhat_i);

		/* Calculate measurement vector */
		catalog_vector[0][0] = catalog_vectors[i][0];
		catalog_vector[1][0] = catalog_vectors[i][1];
		catalog_vector[2][0] = catalog_vectors[i][2];
		mat_mult(A,3,3,catalog_vector,3,1,true_vector);
		y_i = y_0 - true_vector[2][0]/true_vector[0][0]*f;
		x_i = x_0 - true_vector[1][0]/true_vector[0][0]*f;
		deltaz[i*2][0] = x_i - xhat_i;
		deltaz[i*2+1][0] = y_i - yhat_i;
	}

	/* Calculate radial parameters */
	mat_transpose(X,n,2,XT);
	mat_mult(XT,2,n,X,n,2,XTX);
	mat_inverse(XTX,2,invXTX);
	mat_mult(invXTX,2,2,XT,2,n,invXTXXT);
	mat_mult(invXTXXT,2,n,deltaz,n,1,beta);

	/* Assign radial parameters */
	*k_radial[0] = beta[0][0];
	*k_radial[1] = beta[1][0];

	free(star_centroids);
	mat_delete(X,n*2,2);
	mat_delete(deltaz,n*2,1);
	mat_delete(beta,2,1);
	mat_delete(XT,2,n);
	mat_delete(XTX,2,2);
	mat_delete(invXTX,2,2);
	mat_delete(invXTXXT,2,n);
	mat_delete(A,3,3);
	mat_delete(true_vector,3,1);
	mat_delete(catalog_vector,3,1);

}

/**
 * \brief 	Estimates the image FOV using a basic NLLS
 *
 * \param 	*fov 				Estimate image field of view [deg]
 * \param	height				Size of the bitmap image along the height.
 * \param	width				Size of the bitmap image along the width.
 * \param	***image_vectors	Set of recalculated identified vectors according to the new fov
 * 								in the image frame.
 * \param	**catalog_vectors	Set of identified vectors in the catalog frame.
 * \param	*k_radial			Radial distortion parameters. (k_1,k_2)
 * \param	*offset				Image offset from center [px]. (dx,dy)
 * \param	n					Total number of stars to be used in the estimation.
 *
 * \return	void
 */
void calibrate_nlls_fov(double *fov, int height, int width, double ***image_vectors,
		double **catalog_vectors, double *k_radial, double *offset, int n) {

	int i, j, m, pos, count = 0;
	double x_i, y_i, x_j, y_j, N, D_1, D_2, x_0, y_0, f, a, deltaf = 1.0;
	double *z, *zhat, *H, *star_centroids;

	/* Calculate number of star pairs */
	m = nchoosek(n,2);

	/* Assign vectors */
	z = malloc(m*sizeof(double));
	zhat = malloc(m*sizeof(double));
	H = malloc(m*sizeof(double));
	star_centroids = malloc(2*n*sizeof(double));

	/* Calculate focal parameters */
	x_0 = (double) width/2.0 + offset[0];
	y_0 = (double) height/2.0 + offset[1];
	f = width/2.0/tan(*fov*PI/180.0/2.0);

	/* Retrieve measured or distortion calibrated coordinates */
	for (i = 0; i < n; i++) {
		star_centroids[i*2] = y_0 - (*image_vectors)[i][2]/(*image_vectors)[i][0]*f;
		star_centroids[i*2+1] = x_0 - (*image_vectors)[i][1]/(*image_vectors)[i][0]*f;
	}

	/* Iterate through until the iteration limit is reached or the tolerance is met */
	while (fabs(deltaf) > TOL_F && count < MAX_NLLS) {

		/* Go through each star pair */
		pos = 0;
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {

				y_i = star_centroids[i*2];
				x_i = star_centroids[i*2+1];
				y_j = star_centroids[j*2];
				x_j = star_centroids[j*2+1];

				/* Calculate angular separation variables */
				N = (x_0-x_i)*(x_0-x_j) + (y_0-y_i)*(y_0-y_j) + f*f;
				D_1 = sqrt((x_0-x_i)*(x_0-x_i) + (y_0-y_i)*(y_0-y_i) + f*f);
				D_2 = sqrt((x_0-x_j)*(x_0-x_j) + (y_0-y_j)*(y_0-y_j) + f*f);

				/* Calculate focal length partial derivative */
				H[pos] = 2.0*f/D_1/D_2 - N*f/(D_1*D_1*D_1*D_2) - N*f/(D_1*D_2*D_2*D_2);

				/* Calculate expected and estimation vectors */
				zhat[pos] = N/D_1/D_2;
				z[pos] = catalog_vectors[i][0]*catalog_vectors[j][0] + catalog_vectors[i][1]*catalog_vectors[j][1] +
						catalog_vectors[i][2]*catalog_vectors[j][2];

				pos++;
			}
		}

		/* Perform non-linear least squares correction */
		a = 0.0;
		for (i = 0; i < m; i++)
			a += H[i]*H[i];
		deltaf = 0.0;
		for (i = 0; i < m; i++)
			deltaf +=  H[i]*(z[i]-zhat[i])/a;
		f += deltaf;

		count++;
	}

	/* Calculate the new fov */
	*fov = 2.0*atan(x_0/f)*180.0/PI;

	/* Recalculate image vectors using the new fov */
	compute_calibrated_vectors(*fov,height,width,star_centroids,n,k_radial,offset,image_vectors);

	free(z);
	free(zhat);
	free(H);
	free(star_centroids);

}

/**
 * \brief 	Estimates the image FOV using Levenberg-Marquardt
 *
 * \param 	*fov 				Estimate image field of view [deg]
 * \param	height				Size of the bitmap image along the height.
 * \param	width				Size of the bitmap image along the width.
 * \param	***image_vectors	Set of recalculated identified vectors according to the new fov
 * 								in the image frame.
 * \param	**catalog_vectors	Set of identified vectors in the catalog frame.
 * \param	*k_radial			Radial distortion parameters. (k_1,k_2)
 * \param	*offset				Image offset from center [px]. (dx,dy)
 * \param	*n					Total number of stars to be used in the estimation.
 *
 * \return	void
 */
//void calibrate_lm_fov(double *fov, int height, int width, double ***image_vectors,
//		double **catalog_vectors, double *k_radial, double *offset, int n) {
//
//	int i, j, m, pos, count = 0;
//	double x_i, y_i, x_j, y_j, N, D_1, D_2, x_0, y_0, f, a, deltaf = 1.0;
//	double nu = 2.0;
//	double *z, *zhat, *H, *star_centroids;
//
//	/* Calculate number of star pairs */
//	m = nchoosek(n,2);
//
//	/* Assign vectors */
//	z = malloc(m*sizeof(double));
//	zhat = malloc(m*sizeof(double));
//	H = malloc(m*sizeof(double));
//	star_centroids = malloc(2*n*sizeof(double));
//
//	/* Calculate focal parameters */
//	x_0 = (double) width/2.0 + offset[0];
//	y_0 = (double) height/2.0 + offset[1];
//	f = width/2.0/tan(*fov*PI/180.0/2.0);
//
//	/* Iterate through until the iteration limit is reached or the tolerance is met */
//	while (deltaf > TOL_F || count < MAX_NLLS) {
//
//		/* Retrieve measured or distortion calibrated coordinates */
//		for (i = 0; i < n; i++) {
//			star_centroids[i*2] = y_0 - (*image_vectors)[i][2]/(*image_vectors)[i][0]*f;
//			star_centroids[i*2+1] = x_0 - (*image_vectors)[i][1]/(*image_vectors)[i][0]*f;
//		}
//
//		/* Go through each star pair */
//		pos = 0;
//		for (i = 0; i < n-1; i++) {
//			for (j = i+1; j < n; j++) {
//
//				y_i = star_centroids[i*2];
//				x_i = star_centroids[i*2+1];
//				y_j = star_centroids[j*2];
//				x_j = star_centroids[j*2+1];
//
//				/* Calculate angular separation variables */
//				N = (x_0-x_i)*(x_0-x_j) + (y_0-y_i)*(y_0-y_j) + f*f;
//				D_1 = sqrt((x_0-x_i)*(x_0-x_i) + (y_0-y_i)*(y_0-y_i) + f*f);
//				D_2 = sqrt((x_0-x_j)*(x_0-x_j) + (y_0-y_j)*(y_0-y_j) + f*f);
//
//				/* Calculate focal length partial derivative */
//				H[pos] = 2.0*f/D_1/D_2 - N*f/(D_1*D_1*D_1*D_2) - N*f/(D_1*D_2*D_2*D_2);
//
//				/* Calculate expected and estimation vectors */
//				zhat[pos] = N/D_1/D_2;
//				z[pos] = catalog_vectors[i][0]*catalog_vectors[j][0] + catalog_vectors[i][1]*catalog_vectors[j][1] +
//						catalog_vectors[i][2]*catalog_vectors[j][2];
//
//				pos++;
//			}
//		}
//
//		/* Perform non-linear least squares correction */
//		a = 0.0;
//		for (i = 0; i < m; i++)
//			a += H[i]*H[i];
//		deltaf = 0.0;
//		for (i = 0; i < m; i++)
//			deltaf +=  H[i]*(z[i]-zhat[i])/a;
//		f += deltaf;
//
//		count++;
//	}
//
//	/* Calculate the new fov */
//	*fov = 2.0*atan(x_0/f)*180.0/PI;
//
//	/* Recalculate image vectors using the new fov */
//	compute_calibrated_vectors(*fov,height,width,star_centroids,n,k_radial,offset,image_vectors);
//
//	free(z);
//	free(zhat);
//	free(H);
//	free(star_centroids);
//
//}

/**
 * \brief 	Estimates the image FOV and image centre offset using a basic NLLS
 *
 * \param 	*fov 				Estimate image field of view [deg]
 * \param	height				Size of the bitmap image along the height.
 * \param	width				Size of the bitmap image along the width.
 * \param	***image_vectors	Set of recalculated identified vectors according to the new fov
 * 								in the image frame.
 * \param	**catalog_vectors	Set of identified vectors in the catalog frame.
 * \param	*k_radial			Radial distortion parameters. (k_1,k_2)
 * \param	**offset				Image offset from center [px]. (dx,dy)
 * \param	*n					Total number of stars to be used in the estimation.
 *
 * \return	void
 */
void calibrate_nlls_focal(double *fov, int height, int width, double ***image_vectors,
		double **catalog_vectors, double *k_radial, double **offset, int n) {

	int i, j, m, pos=0, count=0;
	double x_i, y_i, x_j, y_j, N, D_1, D_2, x_0, y_0, f;
	double z, zhat, *star_centroids;
	mat H, deltaz, deltax;
	mat HT, HTH, invHTH, invHTHHT;

	/* Calculate number of star pairs */
	m = nchoosek(n,2);

	/* Assign vectors */
	H = mat_create(m,3);
	deltaz = mat_create(m,1);
	deltax = mat_create(3,1);
	HT = mat_create(3,m);
	HTH = mat_create(3,3);
	invHTH = mat_create(3,3);
	invHTHHT = mat_create(3,m);
	star_centroids = malloc(2*n*sizeof(double));
	deltax[0][0] = 1.0; deltax[1][0] = 1.0; deltax[2][0] = 1.0;

	/* Calculate focal parameters */
	x_0 = (double) width/2.0 + offset[0][0];
	y_0 = (double) height/2.0 + offset[0][1];
	f = width/2.0/tan(*fov*PI/180.0/2.0);

	/* Retrieve measured or distortion calibrated coordinates */
	for (i = 0; i < n; i++) {
		star_centroids[i*2] = y_0 - (*image_vectors)[i][2]/(*image_vectors)[i][0]*f;
		star_centroids[i*2+1] = x_0 - (*image_vectors)[i][1]/(*image_vectors)[i][0]*f;
	}

	/* Iterate through until the iteration limit is reached or the tolerance is met */
	while (fabs(deltax[0][0]) > TOL_F &&
			sqrt(deltax[1][0]*deltax[1][0] + deltax[2][0]*deltax[2][0]) > TOL_OFFSET &&
			count < MAX_NLLS) {

		/* Go through each star pair */
		pos = 0;
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {

				y_i = star_centroids[i*2];
				x_i = star_centroids[i*2+1];
				y_j = star_centroids[j*2];
				x_j = star_centroids[j*2+1];

				/* Calculate angular separation variables */
				N = (x_0-x_i)*(x_0-x_j) + (y_0-y_i)*(y_0-y_j) + f*f;
				D_1 = sqrt((x_0-x_i)*(x_0-x_i) + (y_0-y_i)*(y_0-y_i) + f*f);
				D_2 = sqrt((x_0-x_j)*(x_0-x_j) + (y_0-y_j)*(y_0-y_j) + f*f);

				/* Calculate focal length partial derivative */
				H[pos][0] = 2.0*f/D_1/D_2 - N*f/(D_1*D_1*D_1*D_2) - N*f/(D_1*D_2*D_2*D_2);
				H[pos][1] = (2*x_0-x_i-x_j)/D_1/D_2 - N*(x_0-x_i)/D_1/D_1/D_1/D_2 -
						N*(x_0-x_i)/D_1/D_2/D_2/D_2;
				H[pos][2] = (2*y_0-y_i-y_j)/D_1/D_2 - N*(y_0-y_i)/D_1/D_1/D_1/D_2 -
						N*(y_0-y_i)/D_1/D_2/D_2/D_2;

				/* Calculate expected and estimation vectors */
				zhat = N/D_1/D_2;
				z = catalog_vectors[i][0]*catalog_vectors[j][0] + catalog_vectors[i][1]*catalog_vectors[j][1] +
						catalog_vectors[i][2]*catalog_vectors[j][2];
				deltaz[pos][0] = z - zhat;
				pos++;
			}
		}

		/* Perform non-linear least squares correction */
		mat_transpose(H,m,3,HT);
		mat_mult(HT,3,m,H,m,3,HTH);
		mat_inverse(HTH,3,invHTH);
		mat_mult(invHTH,3,3,HT,3,m,invHTHHT);
		mat_mult(invHTHHT,3,m,deltaz,m,1,deltax);
		f = f + deltax[0][0];
		x_0 = x_0 + deltax[1][0];
		y_0 = y_0 + deltax[2][0];
		count++;
	}

	/* Calculate the new fov */
	*fov = 2.0*atan(x_0/f)*180.0/PI;
	*offset[0] = x_0 - width/2;
	*offset[0] = y_0 - height/2;

	/* Recalculate image vectors using the new fov */
	compute_calibrated_vectors(*fov,height,width,star_centroids,n,k_radial,*offset,image_vectors);

	free(star_centroids);
	mat_delete(H,m,3);
	mat_delete(deltaz,m,1);
	mat_delete(deltax,3,1);
	mat_delete(HT,3,m);
	mat_delete(HTH,3,3);
	mat_delete(invHTH,3,3);
	mat_delete(invHTHHT,3,m);

}

 /**
  * \brief 	Estimates the image FOV and image centre offset using a basic NLLS
  *
  * \param 	*fov 				Estimate image field of view [deg]
  * \param	height				Size of the bitmap image along the height.
  * \param	width				Size of the bitmap image along the width.
  * \param	**image_vectors		Set of measured star vectors in the image frame.
  * \param	*q					Estimated attitude quaternion.
  * \param	**catalog_vectors	Set of identified vectors in the catalog frame.
  * \param	**k_radial			Set of recalibrated radial distortion parameters. (k_1,k_2)
  * \param	**offset				Image offset from center [px]. (dx,dy)
  * \param	*n					Total number of stars to be used in the estimation.
  *
  * \return	void
  */
 void calibrate_lm_orion(double *fov, int height, int width, double **image_vectors,
 		double *q, double **catalog_vectors, double **k_radial, double *offset, int n) {

 	int i, count=0;
 	double u_0, v_0, f, y_i, x_i, r_i, new_f, new_k_radial[2];
 	double nbar, lambda, res_norm=1.0, new_res_norm, *star_centroids;
 	mat A, J, res, delta;
 	mat JT, N, invN, invNJT, true_vectors;

 	/* Assign vectors */
 	A = mat_create(3,3);
 	J = mat_create(n*2,3);
 	JT = mat_create(3,n*2);
 	N = mat_create(3,3);
 	invN = mat_create(3,3);
 	invNJT = mat_create(3,n*2);
 	res = mat_create(n*2,1);
 	delta = mat_create(3,1);
 	true_vectors = mat_create(3,n);
	star_centroids = malloc(2*n*sizeof(double));
	delta[0][0] = 1.0; delta[1][0] = 1.0; delta[2][0] = 1.0;

 	/* Calculate focal parameters */
 	u_0 = (double) width/2.0 + offset[0];
 	v_0 = (double) height/2.0 + offset[1];
 	f = width/2.0/tan(*fov*PI/180.0/2.0);

 	/* Retrieve measured or distortion calibrated coordinates */
	for (i = 0; i < n; i++) {
		star_centroids[i*2] = v_0 - image_vectors[i][2]/image_vectors[i][0]*f;
		star_centroids[i*2+1] = u_0 - image_vectors[i][1]/image_vectors[i][0]*f;
	}

 	/* Calculate attitude and transform catalog vectors */
 	quat2rot(q,A);
 	for (i = 0; i < n; i++) {
 		true_vectors[0][i] = A[0][0]*catalog_vectors[i][0] + A[0][1]*catalog_vectors[i][1] + A[0][2]*catalog_vectors[i][2];
 		true_vectors[1][i] = A[1][0]*catalog_vectors[i][0] + A[1][1]*catalog_vectors[i][1] + A[1][2]*catalog_vectors[i][2];
 		true_vectors[2][i] = A[2][0]*catalog_vectors[i][0] + A[2][1]*catalog_vectors[i][1] + A[2][2]*catalog_vectors[i][2];
 	}

 	/* Iterate through until the iteration limit is reached or the tolerance is met */
 	while (count < MAX_NLLS && res_norm > TOL_LM_ORION) {

 		/* Go through each star pair and compute the measurement vector and Jacobian */
 		for (i = 0; i < n; i++) {

			/* Calculate residual vector */
			y_i = -true_vectors[2][i]/true_vectors[0][i]*f;
			x_i = -true_vectors[1][i]/true_vectors[0][i]*f;
			r_i = x_i*x_i + y_i*y_i;
			res[i*2][0] = star_centroids[i*2] - (v_0 + y_i*(1 + k_radial[0][0]*r_i + k_radial[0][1]*r_i*r_i));
			res[i*2+1][0] = star_centroids[i*2+1] - (u_0 + x_i*(1 + k_radial[0][0]*r_i + k_radial[0][1]*r_i*r_i));

			/* Calculate Jacobian */
			J[i*2][0] = -true_vectors[2][0]/true_vectors[0][0];
			J[i*2+1][0] = -true_vectors[1][0]/true_vectors[0][0];
			J[i*2][1] = -f/true_vectors[0][0]*y_i*r_i;
			J[i*2+1][1] = -f/true_vectors[0][0]*x_i*r_i;
			J[i*2][2] = -f/true_vectors[0][0]*y_i*r_i*r_i;
			J[i*2+1][2] = -f/true_vectors[0][0]*x_i*r_i*r_i;
 		}

 		/* Calculate N matrix */
 		mat_transpose(J,n*2,3,JT);
 		mat_mult(JT,3,n*2,J,n*2,3,N);

 		/* Initialise LM if first iteration */
 		if (count == 0) {

 			/* Calculate mean of N */
 			nbar = (N[0][0] + N[1][1] + N[2][2])/3;

 			/* Initialise lambda */
 			lambda = 0.001*nbar;
 		}

 		/* Solve the linear system */
 		for (i = 0; i < 3; i++) {
 			N[i][i] += lambda*N[i][i];
 		}
 		mat_inverse(N,3,invN);
 		mat_mult(invN,3,3,JT,3,n*2,invNJT);
 		mat_mult(invNJT,3,n*2,res,n*2,1,delta);
 		new_f = f + delta[0][0];
 		new_k_radial[0] = k_radial[0][0] + delta[1][0];
 		new_k_radial[1] = k_radial[0][1] + delta[2][0];

 		/* Calculate old and new residual norms */
 		res_norm = 0.0;
 		new_res_norm = 0.0;
 		for (i = 0; i < n; i++) {

 			/* Calculate old residual norm */
 			res_norm += res[i*2][0]*res[i*2][0];
 			res_norm += res[i*2+1][0]*res[i*2+1][0];

			/* Calculate new residual norm */
			y_i = -true_vectors[2][i]/true_vectors[0][i]*new_f;
			x_i = -true_vectors[1][i]/true_vectors[0][i]*new_f;
			r_i = x_i*x_i + y_i*y_i;
			new_res_norm += pow(star_centroids[i*2] - (v_0 + y_i*(1 + new_k_radial[0]*r_i + new_k_radial[1]*r_i*r_i)),2);
			new_res_norm += pow(star_centroids[i*2+1] - (u_0 + x_i*(1 + new_k_radial[0]*r_i + new_k_radial[1]*r_i*r_i)),2);
		}
 		res_norm = sqrt(res_norm);
 		new_res_norm = sqrt(new_res_norm);

		/* Check size and re-initialise appropriately */
 		if (res_norm > new_res_norm) {
 			f = new_f;
 			k_radial[0][0] = new_k_radial[0];
 			k_radial[0][1] = new_k_radial[1];
 			lambda = lambda/10;
 		} else {
 			lambda = lambda*10;
 		}

 		count++;
 	}

 	/* Calculate the new fov */
 	*fov = 2.0*atan(u_0/f)*180.0/PI;

 	/* Unassign memory */
 	free(star_centroids);
 	mat_delete(A,3,3);
 	mat_delete(J,n*2,3);
 	mat_delete(JT,3,n*2);
 	mat_delete(N,3,3);
 	mat_delete(invN,3,3);
 	mat_delete(invNJT,3,n*2);
 	mat_delete(res,n*2,1);
 	mat_delete(delta,3,1);
 	mat_delete(true_vectors,3,n);
 }
