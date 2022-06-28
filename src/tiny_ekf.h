/*
 * TinyEKF: Extended Kalman Filter for embedded processors.
 *
 * Original Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

#include "tiny_ekf_unpacked_struct.h"
#include "tiny_ekf_types.h" /* type definition */
/**
	* Initializes an EKF structure.
	* @param ekf pointer to EKF structure to initialize
	* @param n number of state variables
	* @param m number of observables
	*
	* <tt>ekf</tt> should be a pointer to a structure defined as follows, where <tt>N</tt> and </tt>M</tt> are 
	* constants:
	* <pre>
				dim_t n;           // number of state values 
				dim_t m;           // number of observables 

				number_t x[N];     // state vector

				number_t P[N][N];  // prediction error covariance
				number_t Q[N][N];  // process noise covariance 
				number_t R[M][M];  // measurement error covariance

				number_t G[N][M];  // Kalman gain; a.k.a. K

				number_t F[N][N];  // Jacobian of process model
				number_t H[M][N];  // Jacobian of measurement model

				number_t Ht[N][M]; // transpose of measurement Jacobian
				number_t Ft[N][N]; // transpose of process Jacobian
				number_t Pp[N][N]; // P, post-prediction, pre-update

				number_t fx[N];   // output of user defined f() state-transition function
				number_t hx[M];   // output of user defined h() measurement function

			&nbsp; // temporary storage
				number_t tmp0[N][N];
				number_t tmp1[N][Msta];
				number_t tmp2[M][N];
				number_t tmp3[M][M];
				number_t tmp4[M][M];
				number_t tmp5[M]; 
		* </pre>
	*/

void ekf_init(const void * ekf, const dim_t n, const dim_t m);

/**
	* Runs one step of EKF prediction and update. Your code should first build a model, setting
	* the contents of <tt>ekf.fx</tt>, <tt>ekf.F</tt>, <tt>ekf.hx</tt>, and <tt>ekf.H</tt> to appropriate values.
	* @param ekf pointer to structure EKF 
	* @param z array of measurement (observation) values
	* @return 0 on success, 1 on failure caused by non-positive-definite matrix.
	*/
status_t ekf_step(const void * ekf, const number_t * z);
