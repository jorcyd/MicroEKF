/* gps_ekf: TinyEKF test case using You Chong's GPS example:
 * 
 *   http://www.mathworks.com/matlabcentral/fileexchange/31487-extended-kalman-filter-ekf--for-gps
 * 
 * Reads file gps.csv of satellite data and writes file ekf.csv of mean-subtracted estimated positions.
 *
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
 *
 * MIT License
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

#include "tinyekf_config.h"		//Especificação da estrutura do EKF local
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


static void init(ekf_t * ekf)
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

	// initial covariances of state noise, measurement noise
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
	ekf->x[0] = (number_t)-2.168816181271560e+006;
	// ekf->x[0] = (number_t)-1.168816181271560e+006;
	ekf->x[2] =  (number_t)4.386648549091666e+006;
	ekf->x[4] =  (number_t)4.077161596428751e+006;

	// initial velocity (zeroed ?)
	ekf->x[1] = (number_t)0;
	ekf->x[3] = (number_t)0;
	ekf->x[5] = (number_t)0;
	//Assumindo que a porra do sistema não tem velocidade, tem é que ficar parado(?)
	//colocando alguma velocidade inicial pra dar uma graça nisso
	// ekf->x[1] = (number_t)1e+005;
	// ekf->x[3] = (number_t)-2e+005;
	// ekf->x[5] = (number_t)3e+005;

	// clock bias
	ekf->x[6] = (number_t)3.575261153706439e+006;

	// clock drift
	ekf->x[7] = (number_t)4.549246345845814e+001;
}

//As part of the EKF both F and H matrices can change over time
static void model(ekf_t * ekf, number_t SV[4][3])
{ 
	dim_t i, j;

	for (j=0; j<8; j+=2) {
		ekf->fx[j] = ekf->x[j] + T * ekf->x[j+1];
		ekf->fx[j+1] = ekf->x[j+1];
	}

	for (j=0; j<8; ++j) {
		ekf->F[j][j] = 1;
	}

	for (j=0; j<4; ++j) {
		ekf->F[2*j][2*j+1] = T;
	}

	number_t dx[4][3];
	number_t d;

	for (i=0; i<4; ++i) {
		ekf->hx[i] = 0;
		for (j=0; j<3; ++j) {
			d = ekf->fx[j*2] - SV[i][j];
			dx[i][j] = d;
			ekf->hx[i] += d*d;
		}
		ekf->hx[i] = sqrtf(ekf->hx[i]) + ekf->fx[6];
	}

	for (i=0; i<4; ++i) {
		for (j=0; j<3; ++j) {
			ekf->H[i][j*2]  = dx[i][j] / ekf->hx[i];
		}
		ekf->H[i][6] = 1;
	}   
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
	ekf_init(&ekf, Nsta, Mobs);

	// Do local initialization
	init(&ekf);

	// Open input data file
	FILE * ifp = fopen("gps.csv", "r");
	if(!ifp){
		printf("erro ao abrir gps.csv\n");
		//return 0; ?
	}

	// Skip CSV header
	skipline(ifp);

	// Make a place to store the data from the file and the output of the EKF
	size_t samples = 25;
	//size_t reps = 1000;
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
	bool_t first_iter = false;
	//External loop for profiling only
	// while(reps--){
	// Loop till no more data
	for (j=0; j<samples; ++j) {		
		// if(j>=25 && j%25==0){
		// 	rewind(ifp);
		// 	skipline(ifp);
		// }

		readdata(ifp, SV_Pos, SV_Rho);

		model(&ekf, SV_Pos);

		first_iter = j==0;
		chol_status = ekf_step_op(&ekf, SV_Rho,first_iter,true);	//Observation matrix always change
		if(chol_status == ERROR){
			printf("Cholesky inversion failed on step %d \n",j);
		}

		// grab positions & velocities - 
		for (k=0; k<3; ++k){
			Pos_KF[j][k] = ekf.x[2*k];
			Vel_KF[j][k] = ekf.x[2*k+1];
		}
	}
	// 	rewind(ifp);
	// 	skipline(ifp);
	// 	ekf_init(&ekf, Nsta, Mobs);
	// 	// Do local initialization
	// 	init(&ekf);
	// }

	// Compute means of filtered positions
	number_t mean_Pos_KF[3] = {0, 0, 0};
	for (j=0; j<samples; ++j){ 
		for (k=0; k<3; ++k){
			mean_Pos_KF[k] += Pos_KF[j][k];
		}
	}
	for (k=0; k<3; ++k){
		mean_Pos_KF[k] /= 25;
	}


	// Dump filtered positions minus their means
	for (j=0; j<samples; ++j) {
		fprintf(ofp, "%f,%f,%f\n",Pos_KF[j][0]-mean_Pos_KF[0], Pos_KF[j][1]-mean_Pos_KF[1], Pos_KF[j][2]-mean_Pos_KF[2]);
		printf("%f %f %f / ", Pos_KF[j][0], Pos_KF[j][1], Pos_KF[j][2]);
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

	//EMSCRIPTEN LIBS
	// emscripten::val::global("window").call<void>(
	// 	"offerFileAsDownload",
	// 	string("ekf.csv"),
	// 	string("text/csv")
	// );

	return 0;
}
