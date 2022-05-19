/*
 * TinyEKF: Extended Kalman Filter for embedded processors
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "tiny_ekf.h"

/*	Performance Ninja -- Loop Interchange 1 Summary
	https://www.youtube.com/watch?v=G6BbPB37sYg&list=PLRWO2AL1QAV6bJAU2kgB4xfodGID43Y5d */

/*	Cholesky-decomposition matrix-inversion code, adapated from
	http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt */

static status_t choldc1(number_t * a, number_t * p, const dim_t n) {
	dim_t i,j,k;
	number_t sum;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			sum = a[i*n+j];
			for (k = i - 1; k >= 0; k--) {
				sum -= a[i*n+k] * a[j*n+k];
			}
			if (i == j) {
				if (sum <= (number_t)0) {
					return ERROR; /* error */
				}
				p[i] = sqrtf(sum);
			}
			else {
				a[j*n+i] = sum / p[i];
			}
		}
	}

	return SUCCESS; /* success */
}

static status_t choldcsl(const number_t * A, number_t * a, number_t * p, const dim_t n) 
{
	dim_t i,j,k; number_t sum;
	for (i = 0; i < n; i++) 
		for (j = 0; j < n; j++) 
			a[i*n+j] = A[i*n+j];
	if (choldc1(a, p, n)) return ERROR;
	for (i = 0; i < n; i++) {
		a[i*n+i] = 1 / p[i];
		for (j = i + 1; j < n; j++) {
			sum = (number_t)0;
			for (k = i; k < j; k++) {
				sum -= a[j*n+k] * a[k*n+i];
			}
			a[j*n+i] = sum / p[j];
		}
	}

	return SUCCESS; /* success */
}


static status_t cholsl(const number_t * A, number_t * a, number_t * p, const dim_t n) 
{
	dim_t i,j,k;
	if (choldcsl(A,a,p,n)) return ERROR;
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			a[i*n+j] = (number_t)0;//0.0;
		}
	}
	for (i = 0; i < n; i++) {
		a[i*n+i] *= a[i*n+i];
		for (k = i + 1; k < n; k++) {
			a[i*n+i] += a[k*n+i] * a[k*n+i];
		}
		for (j = i + 1; j < n; j++) {
			for (k = j; k < n; k++) {
				a[i*n+j] += a[k*n+i] * a[k*n+j];
			}
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			a[i*n+j] = a[j*n+i];
		}
	}

	return SUCCESS; /* success */
}

#ifdef __DEBUG__
static void dump(number_t * a, const dim_t m, const dim_t n, const char * fmt)
{
	dim_t i,j;

	char f[100];
	sprintf(f, "%s ", fmt);
	for(i=0; i<m; ++i) {
		for(j=0; j<n; ++j)
			printf(f, a[i*n+j]);
		printf("\n");
	}
}
#endif

#if 0
static void zerosk(number_t * a, const  dim_t k)
{
	dim_t j;
	for (j=0; j<k; ++j)
		a[j] = (number_t)0;	
}
#endif

static void zeros(number_t * a, const  dim_t m, const  dim_t n)
{
	dim_t j;
	for (j=0; j<m*n; ++j)
		a[j] = (number_t)0;
}

//C <- C + A * B 	//l,j interchanged
static void macmat(const number_t * a, const number_t * b, number_t * c, const dim_t arows, const dim_t acols, const dim_t bcols)
{
	dim_t i,j,l;

	for(i=0; i<arows; ++i)
		for(l=0; l<acols; ++l)
			for(j=0; j<bcols; ++j)
				c[i*bcols+j] += a[i*acols+l] * b[l*bcols+j];
}

//C <- C - A*B
static void macmats(const number_t * a, const number_t * b, number_t * c, const dim_t arows, const dim_t acols, const dim_t bcols)
{
	dim_t i,j,l;

	for(i=0; i<arows; ++i)
		for(l=0; l<acols; ++l)
			for(j=0; j<bcols; ++j)
				c[i*bcols+j] -= a[i*acols+l] * b[l*bcols+j];
}

static void macvec(const number_t * a, const  number_t * x, number_t * y, const dim_t m, const dim_t n)
{
	dim_t i, j;

	for(i=0; i<m; ++i) {
		for(j=0; j<n; ++j)
			y[i] += x[j] * a[i*n+j];
	}
}

static void transpose(const number_t * a, number_t * at, const dim_t m, const dim_t n)
{
	dim_t i,j;

	for(i=0; i<m; ++i)
		for(j=0; j<n; ++j) {
			at[j*m+i] = a[i*n+j];
		}
}

/* A <- B */
static void copymat(number_t * a, const number_t * b, const dim_t m, const dim_t n)
{        
	dim_t i,j;

	for(i=0; i<m; ++i)
		for(j=0; j<n; ++j)
			a[i*n+j] = b[i*n+j];
}

static void copyvec(number_t * a, const number_t * b, const dim_t n)
{
	dim_t j;

	for(j=0; j<n; ++j)
		a[j] = b[j];
}

/* C <- A - B */
static void sub(const number_t * a, const number_t * b, number_t * c, const dim_t n)
{
	dim_t j;

	for(j=0; j<n; ++j)
		c[j] = a[j] - b[j];
}

static void eye(number_t * a, const dim_t n)
{
	zeros(a,n,n);
	dim_t i;
	for (i=0; i<n; ++i)
		a[i*n+i] = 1;
}

/* TinyEKF code ------------------------------------------------------------------- */

//__attribute__((packed)) on void v
static void unpack(const void * v, unpacked_ekf_t * ekf, const dim_t n, const dim_t m)
{
	/* skip over n, m in data structure */
	byte_t * cptr = (byte_t *)v;
	cptr += 2*sizeof(dim_t);
	ekf->m = m;
	ekf->n = n;

	/* de-index arrays */
	number_t * dptr = (number_t *)cptr;
	ekf->x = dptr;
	dptr += n;
	ekf->P = dptr;
	dptr += n*n;
	ekf->Q = dptr;
	dptr += n*n;
	ekf->R = dptr;
	dptr += m*m;
	ekf->G = dptr;
	dptr += n*m;
	ekf->F = dptr;
	dptr += n*n;
	ekf->H = dptr;
	dptr += m*n;
	ekf->Ht = dptr;
	dptr += n*m;
	ekf->Ft = dptr;
	dptr += n*n;
	ekf->Pp = dptr;
	dptr += n*n;
	ekf->fx = dptr;
	dptr += n;
	ekf->hx = dptr;
	dptr += m;
	ekf->tmp0 = dptr;
	dptr += n*n;
	ekf->tmp1 = dptr;
	dptr += n*m;
	ekf->tmp2 = dptr;
	dptr += m*n;
	ekf->tmp3 = dptr;
	dptr += m*m;
	ekf->tmp4 = dptr;
	dptr += m*m;
	ekf->tmp5 = dptr;
}

void ekf_init(const void * v, const dim_t n, const dim_t m)
{
	/* retrieve n, m and set them in incoming data structure */
	dim_t * ptr = (dim_t *)v;
	*ptr = n;
	ptr++;
	*ptr = m;

	/* unpack rest of incoming structure for initlization */
	unpacked_ekf_t ekf;
	unpack(v, &ekf, n, m);

	/* zero-out matrices */
	zeros(ekf.P, n, n);
	zeros(ekf.Q, n, n);
	zeros(ekf.R, m, m);
	zeros(ekf.G, n, m);
	zeros(ekf.F, n, n);
	zeros(ekf.H, m, n);
}

void ekf_reset_tmp(unpacked_ekf_t ekf, const dim_t n, const dim_t m){
	zeros(ekf.G, n, m);
	zeros(ekf.Pp, n, n);
	zeros(ekf.tmp0, n, n);
	zeros(ekf.tmp1, n, m);
	zeros(ekf.tmp2, m, n);
	//zeros(ekf.tmp3, m, m);	//Not needed
	//zeros(ekf.tmp4, m, m);	//Not needed
	//zerosk(ekf.tmp5, m);		//Not needed
}

/*z = State vector */
status_t ekf_step_op(const void * v, const number_t * z,bool_t F_changed,bool_t H_changed)
{
	/* unpack incoming structure */
	dim_t * ptr = (dim_t *)v;
	dim_t n = *ptr;
	ptr++;
	dim_t m = *ptr;

	unpacked_ekf_t ekf;
	unpack(v, &ekf, n, m);
	ekf_reset_tmp(ekf, n, m);
 
	/* Predict : Predicted (a priori) estimate covariance */
	/* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
	if(F_changed) transpose(ekf.F, ekf.Ft, n, n);	//F^T_{k-1}
	macmat(ekf.F, ekf.P, ekf.tmp0, n, n, n);		//tmp0 = F_{k-1} * P_{k-1}
	copymat(ekf.Pp, ekf.Q, n, n);					//P_k = Q_{k-1}
	macmat(ekf.tmp0, ekf.Ft, ekf.Pp, n, n, n);		//P_k = P_k + tmp0 * F^T_{k-1}

	/* Update : Optimal Kalman gain */
	/* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
	if(H_changed) transpose(ekf.H, ekf.Ht, m, n);	//H^T_k
	macmat(ekf.Pp, ekf.Ht, ekf.tmp1, n, n, m);		//tmp1 = P_k * H^T_k
	macmat(ekf.H, ekf.Pp, ekf.tmp2, m, n, n);		//tmp2 = H_k * P_k
	copymat(ekf.tmp3, ekf.R, m, m);					//tmp3 = R
	macmat(ekf.tmp2, ekf.Ht, ekf.tmp3, m, n, m);	//tmp3 = tmp3 + tmp2 * H^T_k	
	if (cholsl(ekf.tmp3, ekf.tmp4, ekf.tmp5, m) == ERROR) return ERROR;	//tmp4 = tmp3^{-1} | tmp5 = ???
	macmat(ekf.tmp1, ekf.tmp4, ekf.G, n, m, m);		//G_k = tmp1 * tmp4	//Optimal Kalman gain

	/* Update : Innovation or measurement pre-fit residual */
	/* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
	sub(z, ekf.hx, ekf.tmp5, m);					//tmp5 = z_k - h(\hat{x}_k)
	copyvec(ekf.x,ekf.fx, n);						//\hat{x}_k = \hat{x_k}
	macvec(ekf.G, ekf.tmp5, ekf.x, n, m);			//\hat{x}_k = = \hat{x}_k = + G*tmp5		
	//Subtile floating point difference here. - mulvec + add

	/* Update : Updated (a posteriori) estimate covariance */
	/* P_k = (I - G_k H_k) P_k */
	eye(ekf.tmp0, n);								//tmp0 = I
	macmats(ekf.G, ekf.H, ekf.tmp0, n, m, n);		//tmp0 -= G_k * H_k
	zeros(ekf.P, n, n);								//P_k = 0
	macmat(ekf.tmp0, ekf.Pp, ekf.P, n, n, n);		//P_k = tmp0 * P_{k-1}

	/* success */
	return SUCCESS;
}

status_t ekf_step(const void * v, const number_t * z){
	return(ekf_step_op(v,z,true,true));
}
