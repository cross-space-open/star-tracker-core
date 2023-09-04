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

#include "AttitudeDetermination.h"
#include "matrix.h"
#include "svdcmp.h"
#include "jacobi.h"
#include <stdlib.h>
#include <math.h>

/* Mathematical constant pi's approximate value. */
#define PI 3.1415926

/**
 * \brief 	Calculate quaternion from the rotation matrix
 *
 * \param 	**A		Attitude matrix (2D).
 * \param	*q		Quaternion vector to be calculated
 *
 * \return	void
 */
void rot2quat(double **A, double *q) {

	double tr, s;

    tr = A[0][0] + A[1][1] + A[2][2];
    if (tr > 0.0) {
    	s = sqrt(tr + 1.0);
        q[0] = 0.25*s;
        q[1] = (A[1][2] - A[2][1]) / s;
        q[2] = (A[2][0] - A[0][2]) / s;
        q[3] = (A[0][1] - A[1][0]) / s;
    } else if ((A[0][0] > A[1][1]) && (A[0][0] > A[2][2])) {
		s = 2.0*sqrt(1.0 + A[0][0] - A[1][1] - A[2][2]);
		q[1] = 0.25*s;
		q[0] = (A[1][2] - A[2][1]) / s;
		q[3] = (A[2][0] + A[0][2]) / s;
		q[2] = (A[0][1] + A[1][0]) / s;
    } else if (A[1][1] > A[2][2]) {
		s = 2.0*sqrt(1.0 - A[0][0] + A[1][1] - A[2][2]);
		q[2] = 0.25*s;
		q[3] = (A[1][2] + A[2][1]) / s;
		q[0] = (A[2][0] - A[0][2]) / s;
		q[1] = (A[0][1] + A[1][0]) / s;
    } else {
		s = 2.0*sqrt(1.0 - A[0][0] - A[1][1] + A[2][2]);
		q[3] = 0.25*s;
		q[2] = (A[1][2] + A[2][1]) / s;
		q[1] = (A[2][0] + A[0][2]) / s;
		q[0] = (A[0][1] - A[1][0]) / s;
    }
}

/**
 * \brief 	Calculate rotation matrix from the quaternion
 *
 * \param	*q		Quaternion vector
 * \param 	**A		Attitude matrix (2D) to be calculated
 *
 * \return	void
 */
void quat2rot(double *q, double **A) {

    A[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
    A[0][1] = 2*(q[1]*q[2] + q[0]*q[3]);
    A[0][2] = 2*(q[1]*q[3] - q[0]*q[2]);
    A[1][0] = 2*(q[1]*q[2] - q[0]*q[3]);
    A[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
    A[1][2] = 2*(q[2]*q[3] + q[0]*q[1]);
    A[2][0] = 2*(q[1]*q[3] + q[0]*q[2]);
    A[2][1] = 2*(q[2]*q[3] - q[0]*q[1]);
    A[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];

    return;
}

/**
 * \brief 	Transform star vectors to the catalog frame.
 *
 * \param 	**A				Attitude matrix (2D).
 * \param	*image_vector	Vector in image frame.
 * \param	*catalog_vector	Vector in catalog frame.
 *
 * \return	void
 */
void img2cat(double **A, double *image_vector, double *catalog_vector) {

	mat AT;
	AT = mat_create(3,3);

	/* Determine rotation matrix transpose and then multiply by image vector */
	mat_transpose(A,3,3,AT);
	catalog_vector[0] = AT[0][0]*image_vector[0] + AT[0][1]*image_vector[1] +
							AT[0][2]*image_vector[2];
	catalog_vector[1] = AT[1][0]*image_vector[0] + AT[1][1]*image_vector[1] +
							AT[1][2]*image_vector[2];
	catalog_vector[2] = AT[2][0]*image_vector[0] + AT[2][1]*image_vector[1] +
							AT[2][2]*image_vector[2];
	mat_delete(AT,3,3);
	return;
}

/**
 * \brief 	Transform star vectors to the image frame.
 *
 * \param 	**A				Attitude matrix (2D).
 * \param	*catalog_vector	Vector in catalog frame.
 * \param	*image_vector	Vector in image frame.
 *
 * \return	void
 */
void cat2img(double **A, double *catalog_vector, double *image_vector) {

	/* Multiply rotation matrix by catalog vector */
	image_vector[0] = A[0][0]*catalog_vector[0] + A[0][1]*catalog_vector[1] +
							A[0][2]*catalog_vector[2];
	image_vector[1] = A[1][0]*catalog_vector[0] + A[1][1]*catalog_vector[1] +
							A[1][2]*catalog_vector[2];
	image_vector[2] = A[2][0]*catalog_vector[0] + A[2][1]*catalog_vector[1] +
							A[2][2]*catalog_vector[2];

	return;
}

/**
 * \brief 	Calculate determinant of 3x3 matrix.
 *
 * \param 	**A		Attitude matrix (2D).
 *
 * \return	Determinant.
 */
double determinant_3matrix(double** A) {
	double det;

	det = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) -
			A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]) +
			A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);

	return det;
}

/**
 * \brief 	Calculate attitude using the SVD approach.
 *
 * Each image (u_x,u_y,_u_z)^T and catalog (v_x,v_y,v_z)^T vectors are used to
 * estimate the star image attitude by Singular Value Decomposition.
 *
 * Markley, F. L., & Crassidis, J. L. (2014). Fundamentals of spacecraft attitude
 * determination and control (pp. 196-197). New York, NY, USA:: Springer New York.
 *
 * \param 	**A					Attitude matrix (2D).
 * \param	**image_vectors		Set of vectors in the image frame.
 * \param	**catalog_vectors	Set of vectors in the catalog frame.
 * \param	N					Total number of measured vectors.
 *
 * \return	void
 */
void calculate_attitude_svd(double **A, double **image_vectors, double **catalog_vectors, int N) {

	int i;
	mat H,V,VT;
	vect W;
	int det;

	/* Construct covariance matrix between image and catalog star vectors */
	H = mat_create(3,3);
	mat_zeroize(H,3,3);
	for (i = 0; i < N; i++) {
		H[0][0] += image_vectors[i][0]*catalog_vectors[i][0];
		H[0][1] += image_vectors[i][0]*catalog_vectors[i][1];
		H[0][2] += image_vectors[i][0]*catalog_vectors[i][2];
		H[1][0] += image_vectors[i][1]*catalog_vectors[i][0];
		H[1][1] += image_vectors[i][1]*catalog_vectors[i][1];
		H[1][2] += image_vectors[i][1]*catalog_vectors[i][2];
		H[2][0] += image_vectors[i][2]*catalog_vectors[i][0];
		H[2][1] += image_vectors[i][2]*catalog_vectors[i][1];
		H[2][2] += image_vectors[i][2]*catalog_vectors[i][2];
	}

	/* Compute SVD where U is returned as H */
	V = mat_create(3,3);
	VT = mat_create(3,3);
	W = vect_create(3);
	svdcmp(H,3,3,W,V);
	mat_transpose(V,3,3,VT);

	/* Calculate attitude matrix by multiplying A = U*V */
	mat_mult(H,3,3,VT,3,3,A);

	/* Calculate determinant, and if less than 0, flip third column */
	det = determinant_3matrix(A);
	if (det < 0) {
		A[0][2] = -A[0][2];
		A[1][2] = -A[1][2];
		A[2][2] = -A[2][2];
	}

	/* Clear vectors and matrices, and return */
	mat_delete(H,3,3);
	mat_delete(V,3,3);
	mat_delete(VT,3,3);
	vect_delete(W);
	return;
}

/**
 * \brief 	Calculate attitude using the q method approach.
 *
 * Each image (u_x,u_y,_u_z)^T and catalog (v_x,v_y,v_z)^T vectors are used to
 * estimate the star image attitude by Davenport's q method.
 *
 * Markley, F. L., & Crassidis, J. L. (2014). Fundamentals of spacecraft attitude
 * determination and control (pp. 196-197). New York, NY, USA:: Springer New York.
 *
 * \param 	*q					Quaternion vector
 * \param	**image_vectors		Set of vectors in the image frame.
 * \param	**catalog_vectors	Set of vectors in the catalog frame.
 * \param	N					Total number of measured vectors.
 *
 * \return	void
 */
void calculate_attitude_q_method(double *q, double **image_vectors, double **catalog_vectors, int N) {

	int i,j;
	mat B,K,U;
	vect lambda;
	double trB, norm_q;

	/* Construct attitude profile matrix */
	B = mat_create(3,3);
	mat_zeroize(B,3,3);
	for (i = 0; i < N; i++) {
		B[0][0] += image_vectors[i][0]*catalog_vectors[i][0];
		B[0][1] += image_vectors[i][0]*catalog_vectors[i][1];
		B[0][2] += image_vectors[i][0]*catalog_vectors[i][2];
		B[1][0] += image_vectors[i][1]*catalog_vectors[i][0];
		B[1][1] += image_vectors[i][1]*catalog_vectors[i][1];
		B[1][2] += image_vectors[i][1]*catalog_vectors[i][2];
		B[2][0] += image_vectors[i][2]*catalog_vectors[i][0];
		B[2][1] += image_vectors[i][2]*catalog_vectors[i][1];
		B[2][2] += image_vectors[i][2]*catalog_vectors[i][2];
	}

	/* Construct K matrix */
	K = mat_create(4,4);
	trB = B[0][0] + B[1][1] + B[2][2]; K[0][0] = trB;
	K[1][0] = B[1][2] - B[2][1]; K[2][0] = B[2][0] - B[0][2]; K[3][0] = B[0][1] - B[1][0];
	K[0][1] = K[1][0]; K[0][2] = K[2][0]; K[0][3] = K[3][0];
	for (i = 1; i < 4; i++) {
		for (j = 1; j < 4; j++) {
			K[i][j] = B[i-1][j-1] + B[j-1][i-1];
			if (i == j) K[i][j] -= trB;
		}
	}

	/* Compute eigenvalues/eigenvectors */
	U = mat_create(4,4);
	lambda = vect_create(4);
	jacobi(K,4,lambda,U);

	/* Identify quaternion */
	if (lambda[0] > lambda[1] && lambda[0] > lambda[2] &&
	        lambda[0] > lambda[3]) {
	    q[0] = U[0][0]; q[1] = U[1][0]; q[2] = U[2][0]; q[3] = U[3][0];
	} else if (lambda[1] > lambda[2] && lambda[1] > lambda[3]) {
	    q[0] = U[0][1]; q[1] = U[1][1]; q[2] = U[2][1]; q[3] = U[3][1];
	} else if (lambda[2] > lambda[3]) {
	    q[0] = U[0][2]; q[1] = U[1][2]; q[2] = U[2][2]; q[3] = U[3][2];
	} else {
	    q[0] = U[0][3]; q[1] = U[1][3]; q[2] = U[2][3]; q[3] = U[3][3];
	}

	/* Normalise */
    norm_q = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    for (i = 0; i < 4; i++)
            q[i] = q[i]/norm_q;

	/* Clear vectors and matrices, and return */
	mat_delete(K,4,4);
	mat_delete(B,3,3);
	mat_delete(U,4,4);
	vect_delete(lambda);
	return;
}

/**
 * \brief 	Estimates the angular velocity a pair of attitude measurements.
 *
 * \param 	**A1	Attitude matrix at time step 1 (2D).
 * \param	**A2	Attitude matrix at time step 2 (2D).
 * \param	dt		Length of time step [s].
 * \param	**omega Angular velocity vector to be filled (omega_x,omega_y,omega_z) [rad/s].
 *
 * \return	void
 */
void estimate_angular_velocity(double **A1, double **A2, double dt, double *omega) {

	omega[0] = (fmod(atan2(A2[0][1],A2[0][0]),2*PI) - fmod(atan2(A1[0][1],A1[0][0]),2*PI)) / dt;
	omega[1] = (atan2(A2[0][2], sqrt(A2[1][2]*A2[1][2] + A2[2][2]*A2[2][2])) -
			atan2(A1[0][2], sqrt(A1[1][2]*A1[1][2] + A1[2][2]*A1[2][2]))) / dt;
	omega[2] = (fmod(atan2(A2[1][2],A2[2][2]),2*PI) - fmod(atan2(A1[1][2],A1[2][2]),2*PI)) / dt;

	return;
}

/**
 * \brief 	Perform a state error time update of the EKF
 *
 * \param 	x			State vector of the EKF (q_0,q_1,q_2,q_3,omega_x,omega_y,omega_z)
 * \param	P			Covariance matrix
 * \param	Q			Process noise matrix
 * \param	dt			Length of time step [s]
 *
 * \return	void
 */
void ekf_time_update(double *x, double **P, double **Q, double dt) {

	int i, j;
    mat q0, q1, Omega, Q_tran, omegaCross, eye4;
    mat P0, P1, Phi, PhiT, omegaCross2, eye3, temp;
    double norm_omega;

    /* Define states matrix */
    q0 = mat_create(4,1);
    q0[0][0] = x[0]; q0[1][0] = x[1]; q0[2][0] = x[2]; q0[3][0] = x[3];
    q1 = mat_create(4,1);
    P0 = mat_create(6,6);
    for (i = 0; i < 6; i++)
        for (j = 0; j < 6; j++)
            P0[i][j] = P[i][j];

    /* Define matrix Omega */
    Omega = mat_create(4,4);
    mat_zeroize(Omega,4,4);
    Omega[0][1] = -x[4]; Omega[0][2] = -x[5]; Omega[0][3] = -x[6];
    Omega[1][0] =  x[4]; Omega[1][2] =  x[6]; Omega[1][3] = -x[5];
    Omega[2][0] =  x[5]; Omega[2][1] = -x[6]; Omega[2][3] =  x[4];
    Omega[3][0] =  x[6]; Omega[3][1] =  x[5]; Omega[3][2] = -x[4];

    /* Define matrix Eyes */
    eye4 = mat_create(4,4);
    mat_zeroize(eye4,4,4);
    eye4[0][0] = 1; eye4[1][1] = 1; eye4[2][2] = 1; eye4[3][3] = 1;
    eye3 = mat_create(3,3);
    mat_zeroize(eye3,3,3);
    eye3[0][0] = 1; eye3[1][1] = 1; eye3[2][2] = 1;

    /* Define matrix omegaCross */
    omegaCross = mat_create(3,3);
    mat_zeroize(omegaCross,3,3);
    omegaCross[0][1] =  x[6]; omegaCross[0][2] = -x[5];
    omegaCross[1][0] = -x[6]; omegaCross[1][2] =  x[4];
    omegaCross[2][0] =  x[5]; omegaCross[2][1] = -x[4];
    omegaCross2 = mat_create(3,3);
    mat_mult(omegaCross,3,3,omegaCross,3,3,omegaCross2);

    /* Calculate updated quaternion */
    Q_tran = mat_create(4,4);
    norm_omega = sqrt(x[4]*x[4] + x[5]*x[5] + x[6]*x[6]);
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            Q_tran[i][j] = eye4[i][j]*cos(norm_omega*dt/2.0);
            if (norm_omega != 0.0)
            	Q_tran[i][j] += Omega[i][j]*2.0/norm_omega*sin(norm_omega*dt/2.0);
        }
	}
	mat_mult(Q_tran,4,4,q0,4,1,q1);
    x[0] = q1[0][0]; x[1] = q1[1][0]; x[2] = q1[2][0]; x[3] = q1[3][0];

    /* Calculate updated covariance */
    Phi = mat_create(6,6);
    PhiT = mat_create(6,6);
    P1 = mat_create(6,6);
    mat_zeroize(Phi,6,6);
    temp = mat_create(6,6);
    for (i = 0; i <3; i++) {
        for (j = 0; j <3; j++) {
            Phi[i][j] = eye3[i][j];
            Phi[i+3][j+3] = eye3[i][j];
            Phi[i][j+3] = eye3[i][j]*dt;
            if (norm_omega != 0.0) {
            	Phi[i][j] += -omegaCross[i][j]*sin(norm_omega*dt)/norm_omega +
                  omegaCross2[i][j]*(1-cos(norm_omega*dt))/norm_omega/norm_omega;
            	Phi[i][j+3] += omegaCross[i][j]*(1 - cos(norm_omega*dt))/norm_omega/norm_omega*(i==j) -
                  omegaCross2[i][j]*(norm_omega*dt + sin(norm_omega*dt))/pow(norm_omega,3);
            }
        }
    }
    mat_mult(Phi,6,6,P0,6,6,temp);
    mat_transpose(Phi,6,6,PhiT);
    mat_mult(temp,6,6,PhiT,6,6,P1);
    for (i = 0; i < 6; i++)
        for (j = 0; j < 6; j++)
            P[i][j] = P1[i][j] + Q[i][j];

    /* Clear matrix */
    mat_delete(Phi,6,6); mat_delete(PhiT,6,6); mat_delete(P1, 6,6);
    mat_delete(omegaCross,3,3); mat_delete(omegaCross2,3,3);
    mat_delete(eye3,3,3); mat_delete(eye4,4,4);
    mat_delete(q0,4,1); mat_delete(q1,4,1); mat_delete(Omega,4,4);
    mat_delete(Q_tran,4,4); mat_delete(P0,6,6);
	mat_delete(temp, 6, 6);

    return;
}

/**
 * \brief 	Perform state time update of the EKF
 *
 * \param 	x			State vector of the EKF (q_0,q_1,q_2,q_3,omega_x,omega_y,omega_z)
 * \param	P			Covariance matrix
 * \param	Q			Process noise matrix
 * \param	dt			Length of time step [s]
 *
 * \return	void
 */
void ekf_state_time_update(double *x, double **P, double **Q, double dt) {

	int i, j, k;
    mat q0, q1, Omega, Q_tran, eye4;
    mat P0, P1, Phi, PhiT, eye3, temp, dOmega;
    double norm_omega;

    /* Define states matrix */
    q0 = mat_create(4,1);
    q0[0][0] = x[0]; q0[1][0] = x[1]; q0[2][0] = x[2]; q0[3][0] = x[3];
    q1 = mat_create(4,1);
    P0 = mat_create(7,7);
    for (i = 0; i < 7; i++)
        for (j = 0; j < 7; j++)
            P0[i][j] = P[i][j];

    /* Define matrix Omega and the derivative matrices*/
    Omega = mat_create(4,4); dOmega = mat_create(12,4);
    mat_zeroize(Omega,4,4); mat_zeroize(dOmega,12,4);
    Omega[0][1] = -x[4]; Omega[0][2] = -x[5]; Omega[0][3] = -x[6];
    Omega[1][0] =  x[4]; Omega[1][2] = -x[6]; Omega[1][3] =  x[5];
    Omega[2][0] =  x[5]; Omega[2][1] =  x[6]; Omega[2][3] = -x[4];
    Omega[3][0] =  x[6]; Omega[3][1] = -x[5]; Omega[3][2] =  x[4];
    dOmega[0][1] = -0.5; dOmega[1][0] = 0.5; dOmega[2][3] = -0.5; dOmega[3][2] = 0.5;
    dOmega[4][2] = -0.5; dOmega[5][3] = 0.5; dOmega[6][0] = 0.5; dOmega[7][1] = -0.5;
    dOmega[8][3] = -0.5; dOmega[9][2] = -0.5; dOmega[10][1] = 0.5; dOmega[11][0] = 0.5;

    /* Define matrix Eyes */
    eye4 = mat_create(4,4);
    mat_zeroize(eye4,4,4);
    eye4[0][0] = 1; eye4[1][1] = 1; eye4[2][2] = 1; eye4[3][3] = 1;
    eye3 = mat_create(3,3);
    mat_zeroize(eye3,3,3);
    eye3[0][0] = 1; eye3[1][1] = 1; eye3[2][2] = 1;

    /* Calculate updated quaternion */
    Q_tran = mat_create(4,4);
    norm_omega = sqrt(x[4]*x[4] + x[5]*x[5] + x[6]*x[6]);
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            Q_tran[i][j] = eye4[i][j]*cos(norm_omega*dt/2.0);
            if (norm_omega != 0.0)
            	Q_tran[i][j] += 0.5*Omega[i][j]*2.0/norm_omega*sin(norm_omega*dt/2.0);
        }
	}
	mat_mult(Q_tran,4,4,q0,4,1,q1);
    x[0] = q1[0][0]; x[1] = q1[1][0]; x[2] = q1[2][0]; x[3] = q1[3][0];

    /* Calculate updated covariance */
    Phi = mat_create(7,7);
    PhiT = mat_create(7,7);
    P1 = mat_create(7,7);
    mat_zeroize(Phi,7,7);
    temp = mat_create(7,7);
    for (i = 0; i < 4; i++) {
    	for (j = 0; j < 4; j++) {
			Phi[i][j] = Q_tran[i][j];
		}
	}
    if (norm_omega != 0.0) {
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 4; j++) {
				for (k = 0; k < 4; k++) {
					temp[j][k] = -x[i+4]*dt/2/norm_omega*sin(norm_omega*dt/2)*eye4[j][k] +
							(x[i+4]*dt/norm_omega/norm_omega*cos(norm_omega*dt/2) -
									2*x[i+4]/pow(norm_omega,3)*sin(norm_omega*dt/2))*0.5*Omega[j][k] +
							2/norm_omega*sin(norm_omega*dt/2)*dOmega[j+4*i][k];
				}
			}
			for (j = 0; j < 4; j++) {
				for (k = 0; k < 4; k++) {
					Phi[j][i+4] += temp[j][k]*q0[k][0];
				}
			}
		}
	}
    Phi[4][4] = 1; Phi[5][5] = 1; Phi[6][6] = 1;
    mat_mult(Phi,7,7,P0,7,7,temp);
    mat_transpose(Phi,7,7,PhiT);
    mat_mult(temp,7,7,PhiT,7,7,P1);
    for (i = 0; i < 7; i++)
        for (j = 0; j < 7; j++)
            P[i][j] = P1[i][j] + Q[i][j];

    /* Clear matrix */
    mat_delete(Phi,7,7); mat_delete(PhiT,7,7); mat_delete(P1,7,7);
    mat_delete(eye4,4,4); mat_delete(temp,7,7);
    mat_delete(q0,4,1); mat_delete(q1,4,1); mat_delete(Omega,4,4);
    mat_delete(dOmega,12,4); mat_delete(Q_tran,4,4); mat_delete(P0,7,7);
	mat_delete(eye3, 3, 3);

    return;
}

/**
 * \brief     Rotate vectors using quaternion defined attitude matrix.
 *
 * \param     q		Quaternion (q_0,q_1,q_2,q_3)
 * \param     *u    Original vector
 * \param     *v    Transformed vector
 *
 * \return    void
 */
void rotq(mat q, double* u, double* v) {

      mat C;

      /* Define rotation matrix */
      C = mat_create(3,3);
      C[0][0] = q[0][0]*q[0][0] + q[1][0]*q[1][0] -
                              q[2][0]*q[2][0] - q[3][0]*q[3][0];
      C[0][1] = 2*(q[1][0]*q[2][0] + q[0][0]*q[3][0]);
      C[0][2] = 2*(q[1][0]*q[3][0] - q[0][0]*q[2][0]);
      C[1][0] = 2*(q[1][0]*q[2][0] - q[0][0]*q[3][0]);
      C[1][1] = q[0][0]*q[0][0] - q[1][0]*q[1][0] +
                              q[2][0]*q[2][0] - q[3][0]*q[3][0];
      C[1][2] = 2*(q[2][0]*q[3][0] + q[0][0]*q[1][0]);
      C[2][0] = 2*(q[1][0]*q[3][0] + q[0][0]*q[2][0]);
      C[2][1] = 2*(q[2][0]*q[3][0] - q[0][0]*q[1][0]);
      C[2][2] = q[0][0]*q[0][0] - q[1][0]*q[1][0] -
                              q[2][0]*q[2][0] + q[3][0]*q[3][0];

      /* Transform u to v */
      v[0] = C[0][0]*u[0] + C[0][1]*u[1] + C[0][2]*u[2];
      v[1] = C[1][0]*u[0] + C[1][1]*u[1] + C[1][2]*u[2];
      v[2] = C[2][0]*u[0] + C[2][1]*u[1] + C[2][2]*u[2];

	  mat_delete(C, 3, 3);

      return;
}

/**
 * \brief     Perform a state error measurement update of the EKF.
 *
 * \param     x                                       State vector of the EKF (q_0,q_1,q_2,q_3,omega_x,omega_y,omega_z)
 * \param     P                                       Co-variance matrix
 * \param     R                                       Measurement noise matrix
 * \param     **image_vectors         Set of vectors in the image frame.
 * \param     **catalog_vectors       Set of vectors in the catalog frame.
 * \param     N                                       Total number of measured vectors.
 *
 * \return    void
 */
void ekf_measurement_update(double *x, double **P, double sigma_meas,
              double **image_vectors, double **catalog_vectors, int N) {

      int i, j;
      mat q0, P0, P1, H, y, yhat;
      mat K, HT, temp1, temp2, temp3, temp4, eye6, inv, Deltax;
      double norm_q,*v;

      v = malloc(3*sizeof(double));

      /* Define states matrix */
      q0 = mat_create(4,1);
      q0[0][0] = x[0]; q0[1][0] = x[1]; q0[2][0] = x[2]; q0[3][0] = x[3];
      P0 = mat_create(6,6);
      for (i = 0; i < 6; i++)
              for (j = 0; j < 6; j++)
                      P0[i][j] = P[i][j];
      P1 = mat_create(6,6);

      /* Define matrix Eyes */
      eye6 = mat_create(6,6);
      mat_zeroize(eye6,6,6);
      eye6[0][0] = 1; eye6[1][1] = 1; eye6[2][2] = 1; eye6[3][3] = 1;
      eye6[4][4] = 1; eye6[5][5] = 1;

      /* Define Jacobian and measurement vectors */
      yhat = mat_create(3*N,1); y = mat_create(3*N,1);
      H = mat_create(3*N,6); mat_zeroize(H,3*N,6);
      for (i = 0; i < N; i++) {

		  /* Calculate measurement vector */
		  yhat[3*i][0] = image_vectors[i][0];
		  yhat[3*i+1][0] = image_vectors[i][1];
		  yhat[3*i+2][0] = image_vectors[i][2];

		  /* Calculate expected measurement */
		  rotq(q0,catalog_vectors[i],v);
		  y[3*i][0] = v[0];
		  y[3*i+1][0] = v[1];
		  y[3*i+2][0] = v[2];

		  /* Specify Jacobian */
		  H[3*i][1] = -v[2];
		  H[3*i][2] =  v[1];
		  H[3*i+1][0] = v[2];
		  H[3*i+1][2] = -v[0];
		  H[3*i+2][0] = -v[1];
		  H[3*i+2][1] = v[0];
      }

      /* Calculate gain matrix */
      K = mat_create(6,3*N);
      HT = mat_create(6,3*N); mat_transpose(H,3*N,6,HT);
      temp1 = mat_create(3*N,6);
      temp2 = mat_create(3*N,3*N); inv = mat_create(3*N,3*N);
      temp3 = mat_create(6,3*N);
      mat_mult(H,3*N,6,P0,6,6,temp1);
      mat_mult(temp1,3*N,6,HT,6,3*N,temp2);
      for (i = 0; i < 3*N; i++)
              temp2[i][i] += sigma_meas*sigma_meas;
      mat_inverse(temp2,3*N,inv);
      mat_mult(P0,6,6,HT,6,3*N,temp3);
      mat_mult(temp3,6,3*N,inv,3*N,3*N,K);

      /* Calculate covariance update */
      temp4 = mat_create(6,6);
      mat_mult(K,6,3*N,H,3*N,6,temp4);
      for (i = 0; i < 6; i++)
          for (j = 0; j < 6; j++)
              temp4[i][j] = eye6[i][j] - temp4[i][j];
      mat_mult(temp4,6,6,P0,6,6,P1);
      for (i = 0; i < 6; i++)
          for (j = 0; j < 6; j++)
              P[i][j] = P1[i][j];

      /* Calculate state update */
      Deltax = mat_create(6,1);
      for (i = 0; i < 3*N; i++)
              yhat[i][0] -= y[i][0];
      mat_mult(K,6,3*N,yhat,3*N,1,Deltax);
      x[0] = q0[0][0] + 0.5*(-q0[1][0]*Deltax[0][0]
                      - q0[2][0]*Deltax[1][0] - q0[3][0]*Deltax[2][0]);
      x[1] = q0[1][0] + 0.5*(q0[0][0]*Deltax[0][0]
                      - q0[3][0]*Deltax[1][0] + q0[2][0]*Deltax[2][0]);
      x[2] = q0[2][0] + 0.5*(q0[3][0]*Deltax[0][0]
                      + q0[0][0]*Deltax[1][0] - q0[1][0]*Deltax[2][0]);
      x[3] = q0[3][0] + 0.5*(-q0[2][0]*Deltax[0][0]
                      + q0[1][0]*Deltax[1][0] - q0[0][0]*Deltax[2][0]);
      norm_q = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
      for (i = 0; i < 4; i++)
              x[i] = x[i]/norm_q;
      x[4] += Deltax[3][0];
      x[5] += Deltax[4][0];
      x[6] += Deltax[5][0];

      /* Clear matrix */
      free(v);
      mat_delete(q0,4,1); mat_delete(P0,6,6); mat_delete(P1,6,6);
      mat_delete(H,3*N,6); mat_delete(y,3*N,1); mat_delete(yhat,3*N,1);
      mat_delete(K,6,3*N); mat_delete(HT,6,3*N); mat_delete(temp1,3*N,6);
      mat_delete(temp2,3*N,3*N); mat_delete(temp3,6,3*N); mat_delete(inv,3*N,3*N);
      mat_delete(eye6,6,6); mat_delete(Deltax,6,1);
	  mat_delete(temp4, 6, 6);

      return;
}

/**
 * \brief     Perform a state measurement update of the EKF.
 *
 * \param     x                                       State vector of the EKF (q_0,q_1,q_2,q_3,omega_x,omega_y,omega_z)
 * \param     P                                       Co-variance matrix
 * \param     R                                       Measurement noise matrix
 * \param     **image_vectors         Set of vectors in the image frame.
 * \param     **catalog_vectors       Set of vectors in the catalog frame.
 * \param     N                                       Total number of measured vectors.
 *
 * \return    void
 */
void ekf_state_measurement_update(double *x, double **P, double sigma_meas,
              double **image_vectors, double **catalog_vectors, int N) {

      int i, j;
      mat q0, P0, P1, H, y, yhat, dAdq;
      mat K, HT, temp1, temp2, temp3, temp4, eye7, inv, Deltax;
      double norm_q,*v;

      v = malloc(3*sizeof(double));

      /* Define states matrix */
      q0 = mat_create(4,1);
      q0[0][0] = x[0]; q0[1][0] = x[1]; q0[2][0] = x[2]; q0[3][0] = x[3];
      P0 = mat_create(7,7);
      for (i = 0; i < 7; i++)
              for (j = 0; j < 7; j++)
                      P0[i][j] = P[i][j];
      P1 = mat_create(7,7);

      /* Define matrix Eyes */
      eye7 = mat_create(7,7);
      mat_zeroize(eye7,7,7);
      eye7[0][0] = 1; eye7[1][1] = 1; eye7[2][2] = 1; eye7[3][3] = 1;
      eye7[4][4] = 1; eye7[5][5] = 1; eye7[6][6] = 1;

      /* Define Jacobian and measurement vectors */
      yhat = mat_create(2*N,1); y = mat_create(2*N,1);
      H = mat_create(2*N,7); mat_zeroize(H,2*N,7);
      dAdq = mat_create(12,3);
      dAdq[0][0] =  x[0]; dAdq[0][1] =  x[3]; dAdq[0][2] = -x[2];
      dAdq[1][0] = -x[3]; dAdq[1][1] =  x[0]; dAdq[1][2] =  x[1];
      dAdq[2][0] =  x[2]; dAdq[2][1] = -x[1]; dAdq[2][2] =  x[0];
      dAdq[3][0] =  x[1]; dAdq[3][1] =  x[2]; dAdq[3][2] =  x[3];
      dAdq[4][0] =  x[2]; dAdq[4][1] = -x[1]; dAdq[4][2] =  x[0];
      dAdq[5][0] =  x[3]; dAdq[5][1] = -x[0]; dAdq[5][2] = -x[1];
      dAdq[6][0] = -x[2]; dAdq[6][1] =  x[1]; dAdq[6][2] = -x[0];
      dAdq[7][0] =  x[1]; dAdq[7][1] =  x[2]; dAdq[7][2] =  x[3];
      dAdq[8][0] =  x[0]; dAdq[8][1] =  x[3]; dAdq[8][2] = -x[2];
      dAdq[9][0] = -x[3]; dAdq[9][1] =  x[0]; dAdq[9][2] =  x[1];
      dAdq[10][0] =-x[0]; dAdq[10][1] =-x[3]; dAdq[10][2] = x[2];
      dAdq[11][0] = x[1]; dAdq[11][1] = x[2]; dAdq[11][2] = x[3];
      for (i = 0; i < N; i++) {

		  /* Calculate measurement vector */
		  yhat[2*i][0] = -image_vectors[i][1]/image_vectors[i][0];
		  yhat[2*i+1][0] = -image_vectors[i][2]/image_vectors[i][0];

		  /* Calculate expected measurement */
		  rotq(q0,catalog_vectors[i],v);
		  y[2*i][0] = -v[1]/v[0];
		  y[2*i+1][0] = -v[2]/v[0];

		  /* Specify Jacobian */
		  for (j = 0; j < 4; j++) {
			  H[2*i][j] = -2*(dAdq[j*3+1][0]*catalog_vectors[i][0] + dAdq[j*3+1][1]*catalog_vectors[i][1] +
					  	  	  dAdq[j*3+1][2]*catalog_vectors[i][2])*v[0] +
					  2*(dAdq[j*3][0]*catalog_vectors[i][0] + dAdq[j*3][1]*catalog_vectors[i][1] +
							  dAdq[j*3][2]*catalog_vectors[i][2])*v[1]*v[0];
			  H[2*i+1][j] = -2*(dAdq[j*3+2][0]*catalog_vectors[i][0] + dAdq[j*3+2][1]*catalog_vectors[i][1] +
					  	  	  dAdq[j*3+2][2]*catalog_vectors[i][2])*v[0] +
					  2*(dAdq[j*3][0]*catalog_vectors[i][0] + dAdq[j*3][1]*catalog_vectors[i][1] +
							  dAdq[j*3][2]*catalog_vectors[i][2])*v[2]*v[0];
		  }
      }

      /* Calculate gain matrix */
      K = mat_create(7,2*N);
      HT = mat_create(7,2*N); mat_transpose(H,2*N,7,HT);
      temp1 = mat_create(2*N,7);
      temp2 = mat_create(2*N,2*N); inv = mat_create(2*N,2*N);
      temp3 = mat_create(7,2*N);
      mat_mult(H,2*N,7,P0,7,7,temp1);
      mat_mult(temp1,2*N,7,HT,7,2*N,temp2);
      for (i = 0; i < 2*N; i++)
              temp2[i][i] += sigma_meas*sigma_meas;
      mat_inverse(temp2,2*N,inv);
      mat_mult(P0,7,7,HT,7,2*N,temp3);
      mat_mult(temp3,7,2*N,inv,2*N,2*N,K);

      /* Calculate covariance update */
      temp4 = mat_create(7,7);
      mat_mult(K,7,2*N,H,2*N,7,temp4);
      for (i = 0; i < 7; i++)
          for (j = 0; j < 7; j++)
              temp4[i][j] = eye7[i][j] - temp4[i][j];
      mat_mult(temp4,7,7,P0,7,7,P1);
      for (i = 0; i < 7; i++)
          for (j = 0; j < 7; j++)
              P[i][j] = P1[i][j];

      /* Calculate state update */
      Deltax = mat_create(7,1);
      for (i = 0; i < 2*N; i++)
              yhat[i][0] -= y[i][0];
      mat_mult(K,7,2*N,yhat,2*N,1,Deltax);
      x[0] += Deltax[0][0];
      x[1] += Deltax[1][0];
      x[2] += Deltax[2][0];
      x[3] += Deltax[3][0];
      norm_q = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
      for (i = 0; i < 4; i++)
              x[i] = x[i]/norm_q;
      x[4] += Deltax[4][0];
      x[5] += Deltax[5][0];
      x[6] += Deltax[6][0];

      /* Clear matrix */
      free(v);
      mat_delete(q0,4,1); mat_delete(P0,7,7); mat_delete(P1,7,7); mat_delete(dAdq,12,3);
      mat_delete(H,2*N,7); mat_delete(y,2*N,1); mat_delete(yhat,2*N,1);
      mat_delete(K,7,2*N); mat_delete(HT,7,2*N); mat_delete(temp1,2*N,7);
      mat_delete(temp2,2*N,2*N); mat_delete(temp3,7,2*N); mat_delete(inv,2*N,2*N);
      mat_delete(eye7,7,7); mat_delete(Deltax,7,1);
	  mat_delete(temp4, 7, 7);

      return;
}
