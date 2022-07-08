#include "tiny_ekf_types.h"
//Unpacked struct.
typedef struct __attribute__((packed)){

	dim_t n;		/* number of state values */
	dim_t m;		/* number of observables */

	number_t * x;   /* state vector */

	number_t * P;  /* prediction error covariance */
	number_t * Q;  /* process noise covariance */
	number_t * R;  /* measurement error covariance */

	number_t * G;  /* Kalman gain; a.k.a. K */

	number_t * F;  /* Jacobian of process model */
	number_t * H;  /* Jacobian of measurement model */

	number_t * Ht; /* transpose of measurement Jacobian */
	number_t * Ft; /* transpose of process Jacobian */
	number_t * Pp; /* P, post-prediction, pre-update */

	number_t * fx;  /* output of user defined f() state-transition function */
	number_t * hx;  /* output of user defined h() measurement function */

	/* temporary storage */
	number_t * tmp0;
	number_t * tmp1;
	number_t * tmp2;
	number_t * tmp3;
	number_t * tmp4;
	number_t * tmp5;
	number_t * tmp6; 

} unpacked_ekf_t;