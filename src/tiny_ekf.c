/*
 * TinyEKF: Extended Kalman Filter for embedded processors
 *
 * Original Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "tiny_ekf.h"

/* 	Square root algorithm (Kalman Filter)
	https://stats.stackexchange.com/questions/254749/square-root-algorithm-kalman-filter */

/*	Writing Cache Friendly Code - Comparison of Matrix Multiplication
	https://courses.engr.illinois.edu/cs232/sp2009/lectures/X18.pdf */

/*	Cholesky-decomposition matrix-inversion code, adapated from
	http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt */

static status_t choldc1(number_t *__restrict__ a, number_t *__restrict__ p, const dim_t n) {
	dim_t i,j,k;
	number_t acc;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			acc = a[i*n+j];
			for (k = i - 1; k >= 0; k--) {
				acc -= a[i*n+k] * a[j*n+k];
			}
			if (i == j) {
				if (acc <= (number_t)0) {
					return ERROR; /* error */
				}
				p[i] = sqrtf(acc);
			}
			else {
				a[j*n+i] = acc / p[i];
			}
		}
	}

	return SUCCESS; /* success */
}

static status_t choldcsl(const number_t *__restrict__ A, number_t *__restrict__ a, number_t *__restrict__ p, const dim_t n) 
{
	dim_t i,j,k; number_t acc;
	for (i = 0; i < n; i++) 
		for (j = 0; j < n; j++) 
			a[i*n+j] = A[i*n+j];
	if (choldc1(a, p, n)) return ERROR;
	for (i = 0; i < n; i++) {
		a[i*n+i] = 1 / p[i];
		for (j = i + 1; j < n; j++) {
			acc = (number_t)0;
			for (k = i; k < j; k++) {
				acc -= a[j*n+k] * a[k*n+i];
			}
			a[j*n+i] = acc / p[j];
		}
	}

	return SUCCESS; /* success */
}


static status_t cholsl(const number_t *__restrict__ A, number_t *__restrict__ a, number_t *__restrict__ p, const dim_t n) 
{
	dim_t i,j,k;
	number_t acc;
	if (choldcsl(A,a,p,n)) return ERROR;
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			a[i*n+j] = (number_t)0;
		}
	}
	for (i = 0; i < n; i++) {
		acc = a[i*n+i] * a[i*n+i];
		for (k = i + 1; k < n; k++) {
			acc += a[k*n+i] * a[k*n+i];
		}
		a[i*n+i] = acc;
		for (j = i + 1; j < n; j++) {
			acc = 0;
			for (k = j; k < n; k++) {
				acc += a[k*n+i] * a[k*n+j];
			}
			a[i*n+j] += acc;
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
static void dump(number_t *__restrict__ a, const dim_t m, const dim_t n, const char * fmt)
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

static void zeros(number_t *__restrict__ a, const  dim_t m, const  dim_t n)
{
	dim_t j;
	for (j=0; j<m*n; ++j)
		a[j] = (number_t)0;
}

//The standard (i,j,k) tend to scale well for small matrices as it also avoids a redundant write to memory.
//This loop could also be unrolled easier.
//C <- A * B
static void mulmat(const number_t *__restrict__ a, const number_t *__restrict__ b, number_t *__restrict__ c, const dim_t arows, const dim_t acols, const dim_t bcols)
{
	dim_t i,j,l;
	number_t acc;
	for(i=0; i<arows; ++i)
		for(j=0; j<bcols; ++j) {
			acc = 0;
			for(l=0; l<acols; ++l)
				acc += a[i*acols+l] * b[l*bcols+j];
			c[i*bcols+j] = acc;
		}
}

//C <- -A * B
static void mulmats(const number_t *__restrict__ a, const number_t *__restrict__ b, number_t *__restrict__ c, const dim_t arows, const dim_t acols, const dim_t bcols)
{
	dim_t i,j,l;
	number_t acc;
	for(i=0; i<arows; ++i)
		for(j=0; j<bcols; ++j) {
			acc = 0;
			for(l=0; l<acols; ++l)
				acc -= a[i*acols+l] * b[l*bcols+j];
			c[i*bcols+j] = acc;
		}
}

#if 0
static void mulvec(const number_t *__restrict__ a, const  number_t *__restrict__ x, number_t *__restrict__ y, const dim_t m, const dim_t n)
{
	dim_t i, j;
	number_t acc;

	for(i=0; i<m; ++i) {
		acc = 0;
		for(j=0; j<n; ++j)
			acc += x[j] * a[i*n+j];
		y[i] = acc;
	}
}
#endif

static void macvec(const number_t *__restrict__ a, const  number_t *__restrict__ x, number_t *__restrict__ y, const dim_t m, const dim_t n)
{
	dim_t i, j;
	number_t acc;

	for(i=0; i<m; ++i) {
		acc = 0;
		for(j=0; j<n; ++j)
			acc += x[j] * a[i*n+j];
		y[i] += acc;
	}
}

static void transpose(const number_t *__restrict__ a, number_t *__restrict__ at, const dim_t m, const dim_t n)
{
	dim_t i,j;

	for(i=0; i<m; ++i)
		for(j=0; j<n; ++j) {
			at[j*m+i] = a[i*n+j];
		}
}

/* A <- A + B */
static void accum(number_t *__restrict__ a, const number_t *__restrict__ b, const dim_t m, const dim_t n)
{
	int i,j;

	for(i=0; i<m; ++i)
		for(j=0; j<n; ++j)
			a[i*n+j] += b[i*n+j];
}

/* A <- (A+B)/2*/
// #if 0
static void mean(number_t *__restrict__ a, const number_t *__restrict__ b, const dim_t m, const dim_t n)
{
	int i,j;

	for(i=0; i<m; ++i)
		for(j=0; j<n; ++j){
			a[i*n+j] = (a[i*n+j] + b[i*n+j])*0.5f;	//for floats
			//a[i*n+j] = (a[i*n+j] + b[i*n+j])/2;	//-O1
		}
}
// #endif

static void copyvec(number_t *__restrict__ a, const number_t *__restrict__ b, const dim_t n)
{
	dim_t j;

	for(j=0; j<n; ++j)
		a[j] = b[j];
}

#if 0
/* C <- A + B */
static void add(const number_t *__restrict__ a, const number_t *__restrict__ b, number_t *__restrict__ c, const dim_t n)
{
	int j;

	for(j=0; j<n; ++j)
		c[j] = a[j] + b[j];
}
#endif

/* C <- A - B */
static void sub(const number_t *__restrict__ a, const number_t *__restrict__ b, number_t *__restrict__ c, const dim_t n)
{
	dim_t j;

	for(j=0; j<n; ++j)
		c[j] = a[j] - b[j];
}

#if 0
static void negate(number_t *__restrict__ a, const dim_t m, const dim_t n)
{
	dim_t i, j;

	for(i=0; i<m; ++i)
		for(j=0; j<n; ++j)
		a[i*n+j] = -a[i*n+j];
}
#endif

static void mat_addeye(number_t *__restrict__ a, const dim_t n)
{
	dim_t i;
	for (i=0; i<n; ++i)
		a[i*n+i] += 1;
}

/* TinyEKF code ------------------------------------------------------------------- */

//__attribute__((packed)) on void v
static void unpack(const void *__restrict__ v, unpacked_ekf_t * ekf, const dim_t n, const dim_t m)
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

/*z = State vector */
status_t ekf_step_op(const void * v, const number_t * z,bool_t F_changed,bool_t H_changed)
{
	/* unpack incoming structure */

	int * ptr = (int *)v;
	int n = *ptr;
	ptr++;
	int m = *ptr;

	unpacked_ekf_t ekf;
	unpack(v, &ekf, n, m); 

	/* Predict : Update state vector beforehand in case Kalman gain fails to compute */
	copyvec(ekf.x,ekf.fx, n);							//\hat{x}_k = \hat{x_k}
 
 	/* Predict : Predicted (a priori) estimate covariance */
	/* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
	if(F_changed) transpose(ekf.F, ekf.Ft, n, n);		//F^T_{k-1}
	mulmat(ekf.F, ekf.P, ekf.tmp0, n, n, n);			//tmp0 = F_{k-1}*P_{k-1}
	mulmat(ekf.tmp0, ekf.Ft, ekf.Pp, n, n, n);			//P_k = tmp0*F^T_{k-1}
	accum(ekf.Pp, ekf.Q, n, n);							//P_k += Q_{k-1}

	/* Update : Optimal Kalman gain */
	/* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
	if(H_changed) transpose(ekf.H, ekf.Ht, m, n);		//H^T_k
	mulmat(ekf.Pp, ekf.Ht, ekf.tmp1, n, n, m);			//tmp1 = P_k*H^T_k
	mulmat(ekf.H, ekf.Pp, ekf.tmp2, m, n, n);			//tmp2 = H_k*P_k
	mulmat(ekf.tmp2, ekf.Ht, ekf.tmp3, m, n, m);		//tmp3 = tmp2*H^T_k
	accum(ekf.tmp3, ekf.R, m, m);						//tmp3 += R
	if (cholsl(ekf.tmp3, ekf.tmp4, ekf.tmp5, m)) return ERROR;	//tmp4 = tmp3^-1 / tmp5 = ?
	mulmat(ekf.tmp1, ekf.tmp4, ekf.G, n, m, m);			//G_k = tmp1*temp4

	/* Update : Innovation or measurement pre-fit residual */
	/* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
	sub(z, ekf.hx, ekf.tmp5, m);						//tmp5 = z_k - h(\hat{x}_k)
	macvec(ekf.G, ekf.tmp5, ekf.x, n, m);				//\hat{x}_k += G_k*tmp5

	/* Update : Updated (a posteriori) estimate covariance */
	/* P_k = (I - G_k H_k) P_k */
	mulmats(ekf.G, ekf.H, ekf.tmp0, n, m, n);			//tmp0 = -G_k*H_k
	mat_addeye(ekf.tmp0, n);							//tmp0 += I
	mulmat(ekf.tmp0, ekf.Pp, ekf.P, n, n, n);			//P_k+ = tmp0*P_k

	/* Post Update : A classical hack for ensuring at least symmetry is to do cov_plus = (cov_plus + cov_plus')/2 after the covariance update.*/
	/* P_k =() P_k + P_k^T)/2*/
	transpose(ekf.P, ekf.tmp0, n, n);
	//dump(ekf.P,n,n,"|%f|");
	mean(ekf.P, ekf.tmp0, n, n);

	/* success */
	return SUCCESS;
}

status_t ekf_step(const void * v, const number_t * z){
	return(ekf_step_op(v,z,true,true));
}
