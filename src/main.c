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

#define _GNU_SOURCE

#include "ImageProcess.h"
#include "Tetra.h"
#include "AttitudeDetermination.h"
#include "Calibrate.h"
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <string.h>

#define PI 3.1415926
#define IMAGE_DIR "data/"
#define RESULTS_DIR "results.csv"
#define USING_SERIES 0
#define USING_EKF	 0
#define CALIBRATE_RADIAL 0
#define Q1 	1E-6
#define Q2	1E-8
#define P1 	(0.1*PI/180)*(0.1*PI/180)
#define P2 	(0.2*PI/180/3600)* (0.2*PI/180/3600)
#define Ps1 6E-7
#define Ps2 1E-8
#define Ps3	8E-4

/* To be set by configuration file */
#define FOV 19.2270974722634
#define MIN_MIS_PROB 1.0e-5
#define OFFSET_X 0//248.79 //0
#define OFFSET_Y 0//-112.99 //0
#define RADIAL_1 -7.9238e-09//-2.2184e-09 //9.44084710655260e-48
#define RADIAL_2 1.0549e-15// -4.6619e-15 //4.46732482015949e-32
#define MEAS_NOISE 0.01
#define BRIGHT_THRESH 10

int main(int argc,char *argv[] ) {

	/* Declare variables */
    byte *pixels;
    int32 width, height;
    int len, i, num_bright, num_centroid, n = 0, *star_center_pixels, num_matches, *star_ids, s_off, num_succ = 0;
#ifndef WIN32
    int num_images;
#endif
    double *star_centroids, **image_vectors, **catalog_vectors, mismatch_prob, **A1, **A2, *x, succ_perc;
    double **P, **Q, fov = FOV, ra, dec, roll, f;
    double **Ps, **Qs, min_mis_prob = MIN_MIS_PROB;
    double *offset, *k_radial, bright_thresh = BRIGHT_THRESH;
    double meas_noise = MEAS_NOISE;
    int hh, mm, ss, dt = 0, t2 = 0, t1 = 0, img_counter = 0;
    char image_filename[300];
    char *fname;
    DIR *image_dir;
#ifndef WIN32
	struct dirent **name_list;
#endif
	struct dirent *in_file;
    FILE *rptr;

    /* Allocate memory and set up attitude arrays */
    k_radial = (double*) calloc(2,sizeof(double));
    offset = (double*) calloc(2,sizeof(double));
	A1 = (double**) calloc(3,sizeof(double*)); A2 = (double**) calloc(3,sizeof(double*));
	for (i = 0; i < 3; i++) A1[i] = (double*) calloc(3,sizeof(double));
	for (i = 0; i < 3; i++) A2[i] = (double*) calloc(3,sizeof(double));
	x = (double*) calloc(7,sizeof(double));
	P = (double**) calloc(6,sizeof(double*));
	Q = (double**) calloc(6,sizeof(double*));
	for (i = 0; i < 6; i++) {
		P[i] = (double*) calloc(6,sizeof(double));
		Q[i] = (double*) calloc(6,sizeof(double));
	}
	P[0][0]=P1;P[1][1]=P1;P[2][2]=P1;P[3][3]=P2;P[4][4]=P2;P[5][5]=P2;
	Q[0][0]=Q1;Q[1][1]=Q1;Q[2][2]=Q1;Q[3][3]=Q2;Q[4][4]=Q2;Q[5][5]=Q2;
	Ps = (double**) calloc(7,sizeof(double*));
	Qs = (double**) calloc(7,sizeof(double*));
	for (i = 0; i < 7; i++) {
		Ps[i] = (double*) calloc(7,sizeof(double));
		Qs[i] = (double*) calloc(7,sizeof(double));
	}
    k_radial[0] = RADIAL_1; k_radial[1] = RADIAL_2;
    offset[0] = OFFSET_X; offset[1] = OFFSET_Y;
	Ps[0][0]=Ps1;Ps[1][1]=Ps2;Ps[2][2]=3*Ps2;Ps[3][3]=3*Ps2;Ps[4][4]=Ps3;Ps[5][5]=Ps3;Ps[6][6]=Ps3;
	Qs[0][0]=Q1;Qs[1][1]=Q1;Qs[2][2]=Q1;Qs[3][3]=Q2;Qs[4][4]=Q2;Qs[5][5]=Q2;Qs[6][6]=Q2;

    /* Open image directory */
#ifdef WIN32
	if ((image_dir = opendir(IMAGE_DIR)) == NULL) {
    	printf("Failed to open image directory");
        closedir(image_dir);
    	return 0;
    }
#else
	num_images = scandir(IMAGE_DIR, &name_list, 0, alphasort);
	if (num_images < 0){
		perror("scandir");
		return 0;
	}
#endif


    /* Create text file to read */
	rptr = fopen(RESULTS_DIR,"w+");
    fprintf(rptr,"Year,Month,Day,Hour,Minute,Second,q_0,q_1,q_2,q_3,"
    		"omega_x [rad/s],omega_y [rad/s],omega_z [rad/s],f,P,Image#\n");


    /* Go through each image in directory */
#ifdef WIN32
    while ((in_file = readdir(image_dir)))  {
    	printf("%s\n",in_file->d_name);
		len = in_file->d_namlen;
		fname = in_file->d_name;
#else
	while (img_counter < num_images) {
		printf("%s\n",name_list[img_counter]->d_name);
    	len = strlen(name_list[img_counter]->d_name);
    	fname = name_list[img_counter]->d_name;
#endif

		img_counter++;
    	/* Skip parent and current directories */
    	if((!strcmp(fname,".")) || (!strcmp(fname,".."))) {
#ifndef WIN32
			free(name_list[img_counter]);
#endif
    		continue;
    	}

    	/* Skip if not bmp file */
    	if (!strcmp(fname + len - 4, ".Bmp") && !strcmp(fname + len - 4, ".bmp")) {
#ifndef WIN32
			free(name_list[img_counter]);
#endif
			continue;
    	}

    	if (USING_SERIES) {
    		s_off = strcspn(fname,"-");
    		fprintf(rptr,"%.4s,%.2s,%.2s,%.2s,%.2s,%.2s,",fname+s_off+5,fname+s_off+3,
    			fname+s_off+1,fname+s_off+9,fname+s_off+11,
				fname+s_off+13);
    		sscanf(fname+s_off+9,"%2d",&hh);
    		sscanf(fname+s_off+11,"%2d",&mm);
    		sscanf(fname+s_off+13,"%2d",&ss);
            t2 = hh*3600 + mm*60 + ss;
            if (t2 > 86400) t2 = 0;
            dt = t2 - t1;
            t1 = t2;
    	} else {
    		fprintf(rptr,"0,0,0,0,0,0,");
    	}


		/*** PARTIAL IMAGE DOWNLOAD ***/
		/* Load test image */
		sprintf(image_filename, "%s%s", IMAGE_DIR, fname);
		if (!read_image_bmp(image_filename, &pixels, &width, &height)) {
			printf("Star image failed to process.\n\n\n");
	    	fprintf(rptr,"0,0,0,0,0,0,0,%.3f,1.0,%.*s\n",fov,3,fname + 4);
#ifndef WIN32
			free(name_list[img_counter]);
#endif
			continue;
		}

		/*** STAR TRACKING ***/
		/* Identify star centre pixels and calculate centroids */
		if ((num_bright = identify_star_centre_pixels(pixels, width, height,
				bright_thresh,&star_center_pixels)) == 0) {
			printf("Star image was unsuccessful in identifying stars.\n\n\n");
			free(pixels);
			f = width/2.0/tan(fov*PI/180.0/2.0);
	    	fprintf(rptr,"0,0,0,0,0,0,0,%.3f,1.0,%.*s\n",f,3,fname + 4);
#ifndef WIN32
			free(name_list[img_counter]);
#endif
			continue;
		};
		num_centroid = moment(star_center_pixels, num_bright, pixels, height,
				width, &star_centroids);

		/* Perform TETRA star identification */
		if ((mismatch_prob = run_tetra(star_centroids, &fov, min_mis_prob, k_radial,
				&offset, num_centroid, height, width, &image_vectors,
				&catalog_vectors, &num_matches, &star_ids)) > min_mis_prob) {
			printf("Likely Mismatch : %.4e\n",mismatch_prob);
			printf("Star image was unsuccessful in identifying stars.\n\n\n");
			free(pixels);
			free(star_center_pixels);
			free(star_centroids);
			f = width/2.0/tan(fov*PI/180.0/2.0);
	    	fprintf(rptr,"0,0,0,0,0,0,0,%.3f,1.0,%.*s\n",f,3,fname + 4);
#ifndef WIN32
			free(name_list[img_counter]);
#endif
			continue;
		}

		// /* Print star tracking specific data */
		// printf("---STAR TRACKING---\nFOV %.4f\n",fov);
		// for (i = 0; i < num_matches; i++) {
		// 	printf("Star %d (Cat # %d): (x,y)=(%.2f,%.2f)\n",i+1,star_ids[i], star_centroids[i*2+1],
		// 			star_centroids[i*2]);
		// }
		// printf("\n\n");

		/*** AUTONOMOUS ATTITUDE DETERMINATION ***/
		/* Determine if this is the initial attitude */
		if (n > 1 && USING_EKF ){
			ekf_state_time_update(x, Ps, Qs, (double) dt);
			ekf_state_measurement_update(x, Ps, sqrt(meas_noise), image_vectors, catalog_vectors, num_matches);
			quat2rot(x,A2);
			n++;
		} else {
			calculate_attitude_q_method(x, image_vectors, catalog_vectors, num_matches);
			quat2rot(x,A2);
			n++;
			if (n > 1) {
				estimate_angular_velocity(A1, A2, (double) dt, x + 4);
			}
		}

		/* Calculate */
		ra = fmod(atan2(A2[0][1],A2[0][0]),2*PI) * 180 / PI;
		dec = atan2(A2[0][2], sqrt(A2[1][2]*A2[1][2] + A2[2][2]*A2[2][2])) * 180 / PI;
		roll = fmod(atan2(A2[1][2],A2[2][2]),2*PI) * 180 / PI;
		A1[0][0] = A2[0][0]; A1[0][1] = A2[0][1]; A1[0][2] = A2[0][2];
		A1[1][0] = A2[1][0]; A1[1][1] = A2[1][1]; A1[1][2] = A2[1][2];
		A1[2][0] = A2[2][0]; A1[2][1] = A2[2][1]; A1[2][2] = A2[2][2];

		/* Determine radial calibration */
		if (CALIBRATE_RADIAL && n==1) {
			calibrate_lls_radial(fov, height, width, image_vectors,
					x, catalog_vectors, &k_radial, offset, num_matches);
//			calibrate_lm_orion(&fov, height, width, image_vectors,
//			 		x, catalog_vectors, &k_radial, offset, num_matches);
		}

		/* Calculate percentage success */
		num_succ += 1;
		succ_perc = num_succ/(double)(img_counter-2)*100;

		/* Print star tracking specific data */
		printf("---AUTONOMOUS ATTITUDE DETERMINATION---\n");
		printf("Attitude: ra = %.4f deg dec = %.4f deg roll = %.4f deg\n",ra,dec,roll);
		printf("Mismatch probability: %.4e\n",mismatch_prob);
		printf("Success Percent: %.2f%%\n",succ_perc);
		printf("\n\n");
		f = width/2.0/tan(fov*PI/180.0/2.0);
    	fprintf(rptr,"%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.3f,%.4e,%.*s\n",x[0],x[1],x[2],x[3],x[4],x[5],x[6],f,
    			mismatch_prob,3,fname + 4);

		/* Free memory */
		free(pixels);
		free(star_ids);
		free(star_center_pixels);
		free(star_centroids);
		for (i = 0; i < num_centroid; i++) {
			free(image_vectors[i]);
			free(catalog_vectors[i]);
	   }
		free(image_vectors);
		free(catalog_vectors);
#ifndef WIN32
		free(name_list[img_counter]);
#endif
	}
	for (i = 0; i < 3; i++) {
		free(A1[i]);
		free(A2[i]);
	}
	for (i = 0; i < 6; i++) {
		free(P[i]);
		free(Q[i]);
	}
	for (i = 0; i < 7; i++) {
		free(Ps[i]);
		free(Qs[i]);
	}
	free(A1); free(A2); free(x); free(P); free(Q); free(Ps); free(Qs);
#ifndef WIN32
	free(name_list);
#endif
    fclose(rptr);	

    return 1;
};
