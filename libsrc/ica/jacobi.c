#include <math.h>
#include "matrix.h"
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);

/** Calculate eigenvalues and eigenvectors by the Jacobi method.
 *
 * Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. On
 * output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
 * v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of a.
 *
 *
 * Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (1988). Numerical recipes in C.*/
int jacobi(mat a, int n, vect d, mat v)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c;
	vect b,z;
	b=vect_create(n);
	z=vect_create(n);

	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}

	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	for (i=1;i<=50;i++) {
		sm=0.0;

		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}

		if (sm == 0.0) {
			vect_delete(z);
			vect_delete(b);

			return 0;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);

				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
						&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));

						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;

					tau=s/(1.0+c);

					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;

					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;

					for (j=0;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];

			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}

	return -1;
}
