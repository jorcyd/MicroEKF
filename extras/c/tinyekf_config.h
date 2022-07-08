/*
 * TinyEKF: Extended Kalman Filter for embedded processors.
 *
 * tinyekf_config.h: static configuration parameters
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */
//O tipo tem que tá configurado certo tbm
#include "tiny_ekf_types.h"

/* states */
#define Nsta 8

/* observables */
#define Mobs 4

//depende da plataforma, mas no caso do WASM dá pra forçar.
typedef struct __attribute__((packed)){

	dim_t n;           /* number of state values */
	dim_t m;           /* number of observables */

	number_t x[Nsta];     /* state vector */

	number_t P[Nsta][Nsta];  /* prediction error covariance */
	number_t Q[Nsta][Nsta];  /* process noise covariance */
	number_t R[Mobs][Mobs];  /* measurement error covariance */

	number_t G[Nsta][Mobs];  /* Kalman gain; a.k.a. K */

	number_t F[Nsta][Nsta];  /* Jacobian of process model */
	number_t H[Mobs][Nsta];  /* Jacobian of measurement model */

	number_t Ht[Nsta][Mobs]; /* transpose of measurement Jacobian */
	number_t Ft[Nsta][Nsta]; /* transpose of process Jacobian */
	number_t Pp[Nsta][Nsta]; /* P, post-prediction, pre-update */

	number_t fx[Nsta];   /* output of user defined f() state-transition function */
	number_t hx[Mobs];   /* output of user defined h() measurement function */

	/* temporary storage */
	number_t tmp0[Nsta][Nsta];
	number_t tmp1[Nsta][Mobs];
	number_t tmp2[Mobs][Nsta];
	number_t tmp3[Mobs][Mobs];
	number_t tmp4[Mobs][Mobs];
	number_t tmp5[Mobs]; 
	number_t tmp6[Nsta][Nsta];

} ekf_t;
