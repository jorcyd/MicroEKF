/* gps_ekf: TinyEKF test case using You Chong's GPS example:
 * 
 * http://www.mathworks.com/matlabcentral/fileexchange/31487-extended-kalman-filter-ekf--for-gps
 * https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/31487/versions/6/previews/EKF/GPS_EKF.m/index.html
 *
 * Example:
 * 
 * Reads file gps.csv of satellite data and writes file ekf.csv of mean-subtracted estimated positions.
 *
 * Kalman filter for GPS positioning
 * This file provide an example of using the Extended_KF function with the 
 * the application of GPS navigation. The pseudorange and satellite position
 * of a GPS receiver at fixed location for a period of 25 seconds is
 * provided. Least squares and Extended KF are used for this task.
 *
 * The following is a brief illustration of the principles of GPS. For more
 * information see reference [2].
 * The Global Positioning System(GPS) is a satellite-based navigation system
 * that provides a user with proper equipment access to positioning
 * information. The most commonly used approaches for GPS positioning are
 * the Iterative Least Square(ILS) and the Kalman filtering(KF) methods. 
 * Both of them is based on the pseudorange equation:
 *                rho = || Xs - X || + b + v
 * in which Xs and X represent the position of the satellite and
 * receiver, respectively, and || Xs - X || represents the distance between 
 * them. b represents the clock bias of receiver, and it need to be solved 
 * along with the position of receiver. rho is a measurement given by 
 * receiver for each satellites, and v is the pseudorange measurement noise 
 * modeled as white noise.
 * There are 4 unknowns: the coordinate of receiver position X and the clock
 * bias b. The ILS can be used to calculate these unknowns and is
 * implemented in this example as a comparison. In the KF solution we use
 * the Extended Kalman filter (EKF) to deal with the nonlinearity of the
 * pseudorange equation, and a CV model (constant velocity)[1] as the process
 * model.
 *
 * References:
 *
 * 1. R G Brown, P Y C Hwang, "Introduction to random signals and applied 
 * Kalman filtering : with MATLAB exercises and solutions",1996
 *
 * 2. Pratap Misra, Per Enge, "Global Positioning System Signals, 
 * Measurements, and Performance(Second Edition)",2006
 * 
 * 3. Greg Welch, Gary Bishop , "An Introduction to the Kalman Filter"
 * Department of Computer Science, University of North Carolina at Chapel Hill, 2006
 *
 * Original Copyright (C) : 2015 Simon D. Levy
 * MIT License
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strings.h>
#include "ekf_math.h"
#include "tinyekf_config.h"		//Local ekf structure.
#include "tiny_ekf.h"

// positioning interval
static const number_t T = 1;

static void blkfill(ekf_t * ekf, const number_t * a, int off)
{
	off *= 2;

	ekf->Q[off]   [off]   = a[0]; 
	ekf->Q[off]   [off+1] = a[1];
	ekf->Q[off+1] [off]   = a[2];
	ekf->Q[off+1] [off+1] = a[3];
}


static void model_init(ekf_t * ekf)
{
	// Set Q, see [1]
	const number_t Sf    = 36;
	const number_t Sg    = 0.01;
	const number_t sigma = 5;         		// state transition variance
	const number_t Qb[4] = {Sf*T+Sg*T*T*T/3, Sg*T*T/2, Sg*T*T/2, Sg*T};
	const number_t Qxyz[4] = {sigma*sigma*T*T*T/3, sigma*sigma*T*T/2, sigma*sigma*T*T/2, sigma*sigma*T};

	blkfill(ekf, Qxyz, 0);
	blkfill(ekf, Qxyz, 1);
	blkfill(ekf, Qxyz, 2);
	blkfill(ekf, Qb,   3);

	// initial covariances of state noise, measurement noise (how to estimat ?)
	number_t P0 = 10;		//state noise
	number_t R0 = 36;		//measurament noise

	dim_t i;

	for (i=0; i<8; ++i){
		ekf->P[i][i] = P0;
	}

	for (i=0; i<4; ++i){
		ekf->R[i][i] = R0;
	}

	// initial position
	ekf->x[0] =  (number_t)-2.168816181271560e+006;
	ekf->x[2] =  (number_t)4.386648549091666e+006;
	ekf->x[4] =  (number_t)4.077161596428751e+006;

	// initial zeroed velocity (zeroed ?)
	ekf->x[1] = (number_t)0;
	ekf->x[3] = (number_t)0;
	ekf->x[5] = (number_t)0;
	// Assuming some initial velocity in this case (?)

	// clock bias
	ekf->x[6] = (number_t)3.575261153706439e+006;
	// clock drift
	ekf->x[7] = (number_t)4.549246345845814e+001;
}

static void update_fx(ekf_t * ekf){
	dim_t j;

	for (j=0; j<8; j+=2) {
		ekf->fx[j] = ekf->x[j] + T * ekf->x[j+1];
		ekf->fx[j+1] = ekf->x[j+1];
	}
}

static void update_F(ekf_t * ekf){
	dim_t j;
	for (j=0; j<8; ++j) {
		ekf->F[j][j] = 1;
	}

	for (j=0; j<4; ++j) {
		ekf->F[2*j][2*j+1] = T;
	}
}

static void update_H(ekf_t * ekf, number_t SV[4][3]){
	dim_t i, j;
	number_t dx[4][3];
	number_t d;
	number_t hx;

	for (i=0; i<4; ++i) {
		hx = 0;
		for (j=0; j<3; ++j) {
			d = ekf->fx[j*2] - SV[i][j];
			dx[i][j] = d;
			hx += d*d;
		}
		ekf->hx[i] = fast_sqrtf(hx) + ekf->fx[6];	//this requires a more precise sqrt (would diverge otherwise)
	}

	for (i=0; i<4; ++i) {
		for (j=0; j<3; ++j) {
			ekf->H[i][j*2]  = dx[i][j] * fast_inv(ekf->hx[i]);
		}
		ekf->H[i][6] = 1;
	} 
}

static void model_step(ekf_t * ekf, number_t SV[4][3]) { 
	update_fx(ekf);
	update_F(ekf);
	update_H(ekf,SV);
}

static void readline(char * line, FILE * fp)
{
	fgets(line, 1000, fp);
}

static void readdata(FILE * fp, number_t SV_Pos[4][3], number_t SV_Rho[4])
{
	char line[1000];

	readline(line, fp);
	//printf("%s",line);

	char * p = strtok(line, ",");

	dim_t i, j;

	for (i=0; i<4; ++i) {
		for (j=0; j<3; ++j) {
			SV_Pos[i][j] = strtof(p,NULL);//strtod(p,NULL);//atof(p);
			p = strtok(NULL, ",");
		}
	}

	for (j=0; j<4; ++j) {
		SV_Rho[j] = strtof(p,NULL);//strtod(p,NULL);//atof(p);
		p = strtok(NULL, ",");
	}
}


static void skipline(FILE * fp)
{
	char line[1000];
	readline(line, fp);
}

void error(const char * msg)
{
	fprintf(stderr, "%s\n", msg);
}

int main(int argc, char ** argv)
{    
	// Do generic EKF initialization
	ekf_t ekf;
	unpacked_ekf_t un_ekf;
	// ekf_init(&ekf, Nsta, Mobs);
	ekf_init_ext(&ekf, Nsta, Mobs, &un_ekf);

	// Do local initialization
	model_init(&ekf);

	// Open input data file
	FILE * ifp = fopen("gps.csv", "r");
	if(!ifp){
		printf("erro ao abrir gps.csv\n");
		//return 0; ?
	}

	// Skip CSV header
	skipline(ifp);

	// Make a place to store the data from the file and the output of the EKF
	//size_t samples = 250;
	//size_t reps = 1000;
	size_t samples = 600;
	//size_t samples = 2;
	size_t samples_in_file = 60;
	number_t SV_Pos[4][3];
	number_t SV_Rho[4];
	number_t Pos_KF[samples][3];
	number_t Vel_KF[samples][3];

	// Open output CSV file and write header
	const char * OUTFILE = "ekf.csv";
	FILE * ofp = fopen(OUTFILE, "w");
	if(!ofp){
		printf("erro ao abrir ekf.csv\n");
		//return 0; ?
	}

	fprintf(ofp, "X,Y,Z\n");

	dim_t j, k;
	status_t chol_status = 0;
	number_t x_pos,y_pos,z_pos;

	// Loop till no more data
	for (j=0; j<samples; ++j) {		
		if(j>=samples_in_file && j%samples_in_file==0){
			rewind(ifp);
			skipline(ifp);
		}

		readdata(ifp, SV_Pos, SV_Rho);	//SV Pos/Rho = Sensor vector ?
		model_step(&ekf, SV_Pos);

		//first_iter = j==0;
		//chol_status = ekf_step(&ekf, SV_Rho);	//Observation matrix always change
		chol_status = ekf_step_ext(&un_ekf, SV_Rho);
		if(chol_status == ERROR){
			printf("Cholesky inversion failed on step %d \n",j);
		}

		// grab positions & velocities - 
		for (k=0; k<3; ++k){
			Pos_KF[j][k] = ekf.x[2*k];
			Vel_KF[j][k] = ekf.x[2*k+1];
		}
	}

	// Compute means of filtered positions
	number_t mean_Pos_KF[3] = {0, 0, 0};
	for (j=0; j<samples; ++j){ 
		for (k=0; k<3; ++k){
			mean_Pos_KF[k] += Pos_KF[j][k];
		}
	}
	for (k=0; k<3; ++k){
		mean_Pos_KF[k] /= samples;
	}

	// Dump filtered positions minus their means
	// Also print position minus means and velocities.
	for (j=0; j<samples; ++j) {
		x_pos = Pos_KF[j][0]-mean_Pos_KF[0];
		y_pos = Pos_KF[j][1]-mean_Pos_KF[1];
		z_pos = Pos_KF[j][2]-mean_Pos_KF[2];
		fprintf(ofp, "%f,%f,%f\n",x_pos,y_pos,z_pos);
		printf("%f %f %f / ", x_pos,y_pos,z_pos);
		//printf("%f %f %f / ", Pos_KF[j][0],Pos_KF[j][1],Pos_KF[j][2]);
		printf("%f %f %f\n", Vel_KF[j][0], Vel_KF[j][1], Vel_KF[j][2]);
	}
	

	printf("ifp %p\n", ifp);
	printf("ofp %p\n", ofp);

	// Done!
	if(ifp){
	   fclose(ifp);
	}
	if(ofp) {
	   fclose(ofp);
	}
	
	printf("Wrote file %s\n", OUTFILE);

	return 0;
}
