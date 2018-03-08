#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>
#include <string.h>

#include "gl_integration.c" // For Gauss-Legendre integration


static int lmin = 2;
static int lmax = 20;
static int npts_l; 			//  = lmax - lmin + 1;
static int npts_mu; 			//  = 3*lmax/2 + 1;

static double *kvec;
static double *step_k;
static double kmin;
static double kmax = 8.466251e-1;	// this comes from kmax of the transfer functions data
static int npts_k; 			//  = kmax/npts_k;
//static int npts_k = 6000;

static double *xvec;
static double *step_x;
static double xmax = 1.68e4;		// = 2 * tau0
static int npts_x;

static double omega = 1e2;
static double pi = M_PI;
static double deltaphi = 1.5714e-8;

int i, j, k, l, mu, n, o;	//loop indicies. i for x, j for y, o for omega, n for vetorised indices etc.

//static char *data_dir = "/fast/space/projects/dp002/jf334/data/";
static char *data_dir = "./data/";
static char *plot_dir = ".";

static int bessel_npts_l, bessel_npts_x;
static int bessel_lmin, bessel_lmax;
static double *bessel_xvec;
static double **bessel;
static double tau0; 	// distance to the surface of last scattering 

static int transfer_npts_l, transfer_npts_k;
static int transfer_lmin, transfer_lmax;
static double *transfer_kvec;
static double **transfer;

static double **cos_tilde;
static double **sin_tilde;

static int cl_lmin, cl_lmax, cl_npts_l;
static double *C;

static int BN_lmin, BN_lmax, BN_npts_l;
static double *beam;
static double *noise;

static double *C_inv;

static double *gl_nodes, *gl_weights;
static double **legendre;

static int **i_l;
static int npts_x_l;
static int start_i_l;
static int end_i_l;

static int **i_j;
static int npts_x_y;
static int start_i_j;
static int end_i_j;

static double N_cos, N_sin;
static double sigma_cos, sigma_sin;

int min(int x, int y){
	 return (x < y) ? x : y;
}
int max(int x, int y){
	return (x > y) ? x : y;
}

double *create_array(int d1){
	double *arr = (double *) malloc(d1 * sizeof(double *));
	if (arr == NULL){
		printf("Insuffcient memory for creating a %d array. Bye bye\n", d1);
		exit(1);
	}
	return arr;
}

double **create_2D_array(int d1, int d2){
	double **arr = (double **) malloc(d1 * sizeof(double *));
	arr[0] = (double *) malloc(d1 * d2 * sizeof(double));
	if (arr[0] == NULL){
		printf("Insuffcient memory for creating a %d X %d array. Bye bye\n", d1, d2);
		exit(1);
	}
	int i;
	for(i=1; i<d1; i++) arr[i] = arr[i-1] + d2;
	return arr;
}

/*
double **create_aligned_2D_array(int d1, int d2){
	int d1_pad = (d1 + 7) & ~7;
	int d2_pad = (d2 + 7) & ~7;
	double **arr = (double **) malloc(d1_pad * sizeof(double *));
	arr[0] = (double *) malloc(d1_pad * d2_pad * sizeof(double));
	if (arr[0] == NULL){
		printf("Insuffcient memory for creating a %d X %d array. Bye bye\n", d1, d2);
		exit(1);
	}
	int i;
	for(i=1; i<d1; i++) arr[i] = arr[i-1] + d2_pad;
	return arr;
}
*/	

void free_array(double *arr){
	free(arr);
}

void free_2D_array(double **arr){
	free(arr[0]);
	free(arr);
}

void write_2D_array(double **data, double size1, double size2, char *output_filename){

	char filename[100];
	strcat(strcpy(filename, plot_dir), output_filename);
	FILE *output_file = fopen(filename, "w");
	if (output_file == NULL){
		printf("Error opening file %s to write data. Bye bye\n", output_filename);
		exit(1);
	}
	
	int ii,jj;
	for(ii=0; ii<size1; ii++){
		for(jj=0; jj<size2; jj++){
			fprintf(output_file, "%e ", data[ii][jj]);
		}
		fprintf(output_file, "\n");
	}
	fclose(output_file);

}

void make_xvec(){
	// Makes a vector of x using steps and cuts defined below.
//	double step3 = 200.0;
//	double step2 = 50.0;
//	double step1 = 10.0;
	double step3 = 60.0;
	double step2 = 30.0;
	double step1 = 15.0;
//	double step3 = 40.0;
//	double step2 = 15.0;
//	double step1 = 10.0;
	double cut1 = 13000.0;
	double cut2 = 13600.0;
	double cut3 = 14200.0;
	double cut4 = 15000.0;

	double x = 0.0;
	n = 0;
	for(n=0; x<xmax; n++){
		if (x < cut1){
			x += step3;
		}else if (x < cut2){
			x += step2;
		}else if (x < cut3){
			x += step1;
		}else if (x < cut4){
			x += step2;
		}else{
			x += step3;
		}
	}

	npts_x = n;
	xvec = (double *) malloc(npts_x * sizeof(double));

	x = 0.0;
	i = 0;
	for(i=0; i<npts_x; i++){
		xvec[i] = x;
		if (x < cut1){
			x += step3;
		}else if (x < cut2){
			x += step2;
		}else if (x < cut3){
			x += step1;
		}else if (x < cut4){
			x += step2;
		}else{
			x += step3;
		}
	}
	xmax = xvec[npts_x-1];

	step_x = create_array(npts_x);
//	for(i=1; i<npts_x-1; i++) step_x[i] = 0.5 * (xvec[i+1] - xvec[i-1]);
//	step_x[0] = 0.5 * (xvec[1] - xvec[0]);
//	step_x[npts_x-1] = 0.5 * (xvec[npts_x-1] - xvec[npts_x-2]);
	for(i=0; i<npts_x-1; i++) step_x[i] = xvec[i+1] - xvec[i];
	step_x[npts_x-1] = step3;

}

int **vectorise_two_indices(int size1, int size2){
	int **vec_ind = (int **) malloc(size1 * size2 * sizeof(int *));
	vec_ind[0] = (int *) malloc(2 * size1 * size2 * sizeof(int));
	n = 0;
	for(i=0; i<size1; i++){
		for(j=0; j<size2; j++){
			vec_ind[n] = vec_ind[0] + 2*n;
			vec_ind[n][0] = i;
			vec_ind[n][1] = j;
			n++;
		}
	}
	return vec_ind;
}

int **vectorise_three_indices(int size1, int size2, int size3){
	int **vec_ind = (int **) malloc(size1 * size2 * size3 * sizeof(int *));
	vec_ind[0] = (int *) malloc(3 * size1 * size2 * size3 * sizeof(int));
	n = 0;
	for(i=0; i<size1; i++){
		for(j=0; j<size2; j++){
			for(k=0; k<size3; k++){
				vec_ind[n] = vec_ind[0] + 3*n;
				vec_ind[n][0] = i;
				vec_ind[n][1] = j;
				vec_ind[n][2] = k;
				n++;
			}
		}
	}
	return vec_ind;
}

void prompt(){
	/* Prints parameters and computed error bars*/
	printf("Computing error bars for sinusodial shape functions\n");
	printf("*** Parameters ***\n");
	printf("omega = %e\n", omega);
	printf("lmin = %d, lmax = %d\n", lmin, lmax);
	printf("xmax = %e, npts_x = %d\n", xmax, npts_x);
	printf("kmax = %e, npts_k = %d\n", kmax, npts_k);
	printf("npts_mu = %d\n", npts_mu);
	printf("Amplitude factor deltaphi = %e\n", deltaphi);
	printf("Used bessel_4000 data set and transfer functions,  power spectrum from Planck\n");
	printf("\n");
	printf("Results: sigma_cos = %e, sigma_sin = %e\n\n", sigma_cos, sigma_sin);
}

void load_bessel(){

	FILE *bessel_size_file, *bessel_data_file;
	char bessel_size_filename[100], bessel_data_filename[100];
	int bessel_sizes[2], bessel_size;

	strcat(strcpy(bessel_size_filename, data_dir), "bessel_4000_size");
	strcat(strcpy(bessel_data_filename, data_dir), "bessel_4000");
	bessel_size_file = fopen(bessel_size_filename, "r");
	bessel_data_file = fopen(bessel_data_filename, "r");
	if (bessel_size_file == NULL || bessel_data_file == NULL){
		printf("Error loading bessel size and data files. Bye bye\n");
		exit(1);
	}

	fread(bessel_sizes, sizeof(int), 2, bessel_size_file);
	fclose(bessel_size_file);

	bessel_size = bessel_sizes[0] * bessel_sizes[1];	// total size
	// printf("Bessel data file has sizes %d X %d = %d\n", bessel_sizes[0], bessel_sizes[1], bessel_size);

	bessel = create_2D_array(bessel_sizes[0], bessel_sizes[1]);
	fread(bessel[0], sizeof(double), bessel_size, bessel_data_file);
	fclose(bessel_data_file);

	// The first row of bessel data contains the tau0 and x values (linear)
	// while the first column contains the l values
	// It is implied that the bessel function vanishes at x = 0

	bessel_npts_l = bessel_sizes[0] - 1;
	bessel_npts_x = bessel_sizes[1];

	tau0 = bessel[0][0];
	bessel[0][0] = 0;
	bessel_xvec = &bessel[0][0];	// A vector containing x values of the bessel data
	bessel_lmin = (int) bessel[1][0];
	bessel_lmax = (int) bessel[bessel_npts_l][0];
	for(l=1; l<=bessel_npts_l; l++) bessel[l][0] = 0;
//	 printf("tau0 = %e, xmin = %e, xmax = %e, lmin = %d, lmax = %d\n", tau0, bessel_xvec[0], bessel_xvec[bessel_npts_x-1], bessel_lmin, bessel_lmax);
}

void load_transfer(){

	FILE *transfer_size_file, *transfer_data_file;
	char transfer_size_filename[100], transfer_data_filename[100];
	int transfer_sizes[2], transfer_size;

	strcat(strcpy(transfer_size_filename, data_dir), "transfer_planck_DDX9_5000_size");
	strcat(strcpy(transfer_data_filename, data_dir), "transfer_planck_DDX9_5000_T");
	transfer_size_file = fopen(transfer_size_filename, "r");
	transfer_data_file = fopen(transfer_data_filename, "r");
	if (transfer_size_file == NULL || transfer_data_file == NULL){
		printf("Error loading transfer size and data files. Bye bye\n");
		exit(1);
	}

	fread(transfer_sizes, sizeof(int), 2, transfer_size_file);
	fclose(transfer_size_file);

	transfer_size = transfer_sizes[0] * transfer_sizes[1]; // total size
	// printf("Transfer function data file has size %d X %d = %d\n", transfer_sizes[0], transfer_sizes[1], transfer_size);

	transfer = create_2D_array(transfer_sizes[0], transfer_sizes[1]);
	fread(transfer[0], sizeof(double), transfer_size, transfer_data_file);
	fclose(transfer_data_file);

	// The first row of transfer functions data contains values of k (exponential),
	// while the first column of the data contains l values
	// It is implied that the transfer function vanishes at k = 0

	transfer_npts_l = transfer_sizes[0] - 1;
	transfer_npts_k = transfer_sizes[1];

	transfer_kvec = &transfer[0][0];	// A vector containing k values of the transfer functions data
	transfer_lmin = (int) transfer[1][0];
	transfer_lmax = (int) transfer[transfer_npts_l][0];
	for(l=1; l<=transfer_npts_l; l++) transfer[l][0] = 1;
	// printf("kmin = %e, kmax = %e, lmin = %d, lmax = %d\n", transfer_kvec[0], transfer_kvec[transfer_npts_k-1], transfer_lmin, transfer_lmax);

}

void load_cls(){

	FILE *angular_power_spectrum_file;
	char angular_power_spectrum_filename[100];
	char line[200], **cptr = (char **) malloc(1 * sizeof(char *));
	C = create_array(3499);

	strcat(strcpy(angular_power_spectrum_filename, data_dir), "cls_planck_DDX9_scalar.txt");
	angular_power_spectrum_file = fopen(angular_power_spectrum_filename, "r");
	if (angular_power_spectrum_file == NULL){
		printf("Error loading angular power spectrum file. Bye bye\n");
		exit(1);
	}

	l = 0;
	while (fgets(line, 200, angular_power_spectrum_file)){
//		printf("%s ", line);
		cl_lmax = strtod(line, cptr);
//		printf("%f %s\n", cl_lmax, cptr);
		if (l == 0) cl_lmin = cl_lmax;
		C[l] = strtod(*cptr, cptr);
		C[l] = 2e0 * pi * C[l] / ((l + cl_lmin) * (l + cl_lmin + 1e0));
		l++;
	}
	cl_npts_l = l;
	fclose(angular_power_spectrum_file);

	//printf("Angular power spectrum data file has %d points in l\n", cl_npts_l);
	//printf("lmin = %d, lmax = %d\n", cl_lmin, cl_lmax);
	//for (i=0; i<100; i++) printf("%e ", C[i]);
	//printf("\n");

}

void load_BN(){

	FILE *BN_file;	// Beam and Noises
	char BN_filename[100];
	char line[200], **cptr = (char **) malloc(1 * sizeof(char *));
	beam = create_array(2501);
	noise = create_array(2501);

	strcat(strcpy(BN_filename, data_dir), "BN_DX11d_smica_case1_HP_T.txt");
	BN_file = fopen(BN_filename, "r");

	if (BN_file == NULL){
		printf("Error loading beam and noise file. Bye bye\n");
		exit(1);
	}

	l = 0;
	while (fgets(line, 200, BN_file)){
		BN_lmax = strtod(line, cptr);
		if (l == 0) BN_lmin = BN_lmax;
		beam[l] = strtod(*cptr, cptr);
		noise[l] = strtod(*cptr, cptr);
		l++;
	}
	BN_npts_l = l;
	fclose(BN_file);

}

void parameter_init(){

	// Load data
	load_bessel();
	load_transfer();
	load_cls();
	load_BN();

	printf("loaded\n");
	// Initialise some parameters and the x, k vectors
	npts_l = lmax - lmin + 1;
	npts_mu = 3*lmax/ 2 + 1;

	kvec = transfer_kvec;
	npts_k = transfer_npts_k;
	kmin = kvec[0];
	kmax = kvec[npts_k-1];
	step_k = create_array(npts_k);
//	for(k=1; k<npts_k-1; k++) step_k[k] = 0.5 * (kvec[k+1] - kvec[k-1]);
//	step_k[0] = kvec[1] - kvec[0];
//	step_k[npts_k-1] = kvec[npts_k-1] - kvec[npts_k-2];
	for(k=0; k<npts_k-1; k++) step_k[k] = kvec[k+1] - kvec[k];
	step_k[npts_k-1] = step_k[npts_k-2];

	printf("kvec made\n");
	make_xvec();

	// Vectorise indicies for more efficient looping
	i_l = vectorise_two_indices(npts_x, npts_l);
	npts_x_l = npts_x * npts_l;
	start_i_l = 0;
	end_i_l = npts_x_l;

	i_j = vectorise_two_indices(npts_x, npts_x);
	npts_x_y = npts_x * npts_x;
	start_i_j = 0;
	end_i_j = npts_x_y;

	// Allocate memory for tilde arrays
	cos_tilde = create_2D_array(npts_x, npts_l);
	sin_tilde = create_2D_array(npts_x, npts_l);

}


void precompute_C_inv(){
	/* Incorporates the beam and noise effects to the power spectrum, and then computes C_inv[l] = (2l+1)/C[l]	*/
	// Note this index l is the main l, varying from lmin to lmax
	l = 0;
	while ((l < cl_npts_l) && (l + cl_lmin - BN_lmin < BN_npts_l)){
		// C[l] = pow(beam[l + cl_lmin - BN_lmin], 2) * C[l] + noise[l + cl_lmin - BN_lmin];
		C[l] = C[l] + noise[l + cl_lmin - BN_lmin] / pow(beam[l + cl_lmin - BN_lmin], 2);
		l++;
	}

	C_inv = create_array(npts_l);
	for(l=0; l<npts_l; l++){
		C_inv[l] = (2*(l+lmin)+1) / C[l + lmin - cl_lmin];
	}

}


void precompute_tilde(){
	/* Evaluate sin_tilde(x,l) and cos_tilde(x,l) */

	// Initialise to zero
	for(i=0; i<npts_x; i++){
		for(l=0; l<npts_l; l++){
			cos_tilde[i][l] = sin_tilde[i][l] = 0;
		}
	}

	// Perform integral over k to evaluate cos_tilde(x,l) and  sin_tilde(x,l)
	// cos_tilde(x,l) = (2/pi) * integral{dk * cos(omega*k) * j_l(k*x) * Delta_l(k)}
	// sin_tilde(x,l) = (2/pi) * integral{dk * sin(omega*k) * j_l(k*x) * Delta_l(k)}

	#pragma omp parallel private(n,i,l,k)
	{
		// Initialise GSL tools for interpolation and integration
		gsl_spline *bessel_spline = gsl_spline_alloc(gsl_interp_linear, bessel_npts_x);
		gsl_interp_accel *bessel_acc = gsl_interp_accel_alloc();
		gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, npts_k);
		gsl_interp_accel *acc = gsl_interp_accel_alloc();
		double bes, trans;	// temporary variables
		double *y1 = create_array(npts_k);
		double *y2 = create_array(npts_k);

		#pragma omp for
		for(n=start_i_l; n<end_i_l; n++){
			i = i_l[n][0];
			l = i_l[n][1];

			// Note that bessel[1][0] contains the first value for l = bessel_lmin
			// and similarly transfer[1][0]	contains the first value for l = transfer_lmin	
			gsl_spline_init(bessel_spline, bessel_xvec, &bessel[l + lmin - bessel_lmin + 1][0], bessel_npts_x);

			for(k=0; k<npts_k; k++){
				bes = gsl_spline_eval(bessel_spline, xvec[i] * kvec[k], bessel_acc);
				trans = transfer[l + lmin - transfer_lmin + 1][k];

				y1[k] = cos(omega * kvec[k]) * trans * bes;
				y2[k] = sin(omega * kvec[k]) * trans * bes;
			}

			gsl_interp_accel_reset(bessel_acc);

			gsl_spline_init(spline, kvec, y1, npts_k);
			cos_tilde[i][l] = (2e0/pi) * gsl_spline_eval_integ(spline, kmin, kmax, acc);
			gsl_interp_accel_reset(acc);

			gsl_spline_init(spline, kvec, y2, npts_k);
			sin_tilde[i][l] = (2e0/pi) * gsl_spline_eval_integ(spline, kmin, kmax, acc);
			gsl_interp_accel_reset(acc);
		}

		// Free the interpolaters
		gsl_spline_free(bessel_spline);
		gsl_interp_accel_free(bessel_acc);
		gsl_spline_free(spline);
		gsl_interp_accel_free(acc);

	}
}

void compute_error_bars(){
	/* Compute P_ss, P_sc, P_cs, P_cc and use them to compute N and hence sigma */

	// Compute N
	// For S = sin,
	// N = (1/8pi) * integral{x^2*dx * y^2*dy * dmu * [3*P_ss*P_cc*P_cc + 6*P_sc*P_cs*P_cc + 1*P_ss*P_ss*P_ss - 3*P_sc*P_sc*P_ss - 3*P_cs*P_cs*P_ss]}
	// For S = cos,
	// N = (1/8pi) * integral{x^2*dx * y^2*dy * dmu * [3*P_cc*P_ss*P_ss + 6*P_cs*P_sc*P_ss + 1*P_cc*P_cc*P_cc - 3*P_cs*P_cs*P_cc - 3*P_sc*P_sc*P_cc]}
	// where
	// P_ss(x,y,mu) = sum_l{((2l+1)/C_l) * sin_tilde(x,l) * sin_tilde(y,l) * P_l(mu)}
	// P_cs(x,y,mu) = sum_l{((2l+1)/C_l) * cos_tilde(x,l) * sin_tilde(y,l) * P_l(mu)}
	// P_sc(x,y,mu) = sum_l{((2l+1)/C_l) * sin_tilde(x,l) * cos_tilde(y,l) * P_l(mu)}
	// P_cc(x,y,mu) = sum_l{((2l+1)/C_l) * cos_tilde(x,l) * cos_tilde(y,l) * P_l(mu)}


	N_cos = N_sin = 0;
	double pref0 = 1e0/(8e0*pi) * pow(6e0 * deltaphi * deltaphi, 2);

	#pragma omp parallel private(n,i,j,mu,l)
	{
		double cc, cs, sc, ss;
		double pref, leg;
		double *CC = create_array(npts_l);
		double *CS = create_array(npts_l);
		double *SC = create_array(npts_l);
		double *SS = create_array(npts_l);

		#pragma omp for reduction(+:N_cos,N_sin)
		for(n=start_i_j; n<end_i_j; n++){
			i = i_j[n][0];
			j = i_j[n][1];

			for(l=0; l<npts_l; l++){
				// Compute these combinations first to save computation time
				CC[l] = C_inv[l] * cos_tilde[i][l] * cos_tilde[i][l];
				CS[l] = C_inv[l] * cos_tilde[i][l] * sin_tilde[i][l];
				SC[l] = C_inv[l] * sin_tilde[i][l] * cos_tilde[i][l];
				SS[l] = C_inv[l] * sin_tilde[i][l] * sin_tilde[i][l];
			}

			pref = step_x[i] * step_x[j] * (xvec[i]*xvec[i]) * (xvec[j]*xvec[j]);	// integration prefactor
			
			for(mu=0; mu<npts_mu; mu++){
				// Sum over l first
				cc = cs = sc = ss = 0;
				for(l=0; l<npts_l; l++){
					cc += CC[l] * legendre[mu][l];
					cs += CS[l] * legendre[mu][l];
					sc += SC[l] * legendre[mu][l];
					ss += SS[l] * legendre[mu][l];
				}	

				N_cos += gl_weights[mu] * pref * ((cc * cc * cc) + 3*(cc * ss * ss) - 3*(sc * sc * cc) - 3*(cs * cs * cc) + 6*(sc * cs * ss));
				N_sin += gl_weights[mu] * pref * ((ss * ss * ss) + 3*(ss * cc * cc) - 3*(cs * cs * ss) - 3*(sc * sc * ss) + 6*(cs * sc * cc));
			}

		}

		free(CC);
		free(CS);
		free(SC);
		free(SS);
	}

	N_cos = pref0 * N_cos; 
	N_sin = pref0 * N_sin;

	sigma_cos = sqrt(6e0/N_cos);
	sigma_sin = sqrt(6e0/N_sin);

}

int multiple_omega(){
	/* Calculates error bars for multiple omegas */

	// Initialise
	parameter_init();

	precompute_C_inv();
	free_array(C);
	free_array(beam);
	free_array(noise);

	// Use the Gauss-Legendre method for integrating mu
	gl_nodes = create_array(npts_mu);
	gl_weights = create_array(npts_mu);
	asy(gl_nodes, gl_weights, npts_mu); // Implemented in gl_integration.c

	// Values of Legendre polynomials at the nodes are computed and stored 
	legendre = create_2D_array(npts_mu, npts_l);
	#pragma omp parallel for private(mu,l)
	for(mu=0; mu<npts_mu; mu++){
		for(l=0; l<npts_l; l++){
			legendre[mu][l] = gsl_sf_legendre_Pl(lmin + l, gl_nodes[mu]);
		}
	}


	int npts_omega = 10;
	double results[npts_omega][2];
	for(o=0; o<npts_omega; o++){
		omega = 10.0 * (o + 1);
		precompute_tilde();
		compute_error_bars();
		results[o][0] = sigma_cos;
		results[o][1] = sigma_sin;
		prompt();

		free_2D_array(cos_tilde);
		free_2D_array(sin_tilde);
	}

	printf("\n\n sigma_cos: ");
	for(o=0; o<npts_omega; o++) printf("%e ", results[o][0]);
	printf("\n sigma_sin: ");
	for(o=0; o<npts_omega; o++) printf("%e ", results[o][1]);
	printf("\n");
	
}

void bispectrum_plot(){
	// Calculates bispectrum for l1=l2=l3=l and saves the results to a file

	// Initialise
	parameter_init();

	precompute_tilde();

	free_2D_array(bessel);
	free_2D_array(transfer);

	double *b_cos = (double *) calloc(npts_l, sizeof(double));
	double *b_sin = (double *) calloc(npts_l, sizeof(double));
	
	double c, s;
	#pragma omp parallel for private(l,i,c,s)
	for(l=0; l<npts_l; l++){
		for(i=0; i<npts_x; i++){
			c = cos_tilde[i][l];
			s = sin_tilde[i][l];
			b_cos[l] += xvec[i] * xvec[i] * step_x[i] * (c*c*c - 3*c*s*s);
			b_sin[l] += xvec[i] * xvec[i] * step_x[i] * (-s*s*s + 3*s*c*c);
		}
	}



	FILE *bispectrum_output_file = fopen("./bispectrum_output.txt", "w");
	if (bispectrum_output_file == NULL){
		printf("Error opening file to write bispectrum data. Bye bye\n");
		exit(1);
	}

	for(l=0; l<npts_l; l++){
		fprintf(bispectrum_output_file, "%d %e %e\n", l+lmin, b_cos[l], b_sin[l]);
	}
	fclose(bispectrum_output_file);

	free(b_cos);
	free(b_sin);
}


int error_bars(){

	clock_t start = clock();
	
	// Initialise and load data
	parameter_init();
	
	// Precompute sin_tilde(x) and cos_tilde(x)
	precompute_tilde();
	
	free_2D_array(bessel);
	free_2D_array(transfer);

	// Incorporate BN effects to the power spectrum and compute (2l+1)/Cl
	precompute_C_inv();
	printf("Finished computing sin_tilde, cos_tilde. Elapsed time: %fs\n", (double)(clock() - start) / CLOCKS_PER_SEC);
	
	free_array(C);
	free_array(beam);
	free_array(noise);
	
	// Use the Gauss-Legendre method for integrating mu
	gl_nodes = create_array(npts_mu);
	gl_weights = create_array(npts_mu);
	asy(gl_nodes, gl_weights, npts_mu); // Implemented in gl_integration.c
	
	// Values of Legendre polynomials at the nodes are computed and stored 
	legendre = create_2D_array(npts_mu, npts_l);
	#pragma omp parallel for private(mu,l)
	for(mu=0; mu<npts_mu; mu++){
		for(l=0; l<npts_l; l++){
			legendre[mu][l] = gsl_sf_legendre_Pl(lmin + l, gl_nodes[mu]);
		}
	}
	
	
	// Main computation
	compute_error_bars();
	
	free_2D_array(cos_tilde);
	free_2D_array(sin_tilde);
	free_array(C_inv);
	free_2D_array(legendre);
	free_array(gl_nodes);
	free_array(gl_weights);
	free(i_j[0]);
	free(i_j);	

	// Print out the answer
	prompt();
	
	clock_t end = clock();
	double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Elapsed time: %fs\n", time_spent);
	
	return 0;
	
}

int mpi_error_bars(int argc, char **argv){

	// MPI parameters
	int rank, nprocs;

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Initialise parameters and load data
	parameter_init();

//	printf("Proc %d finished initialising\n", rank);

	// Distribute the workload by setting start and end points for vectorised indices.
	// Note that for non-mpi codes they are set to 0 and npts, respectively
	start_i_l = rank * (npts_x_l / nprocs) + fmin(rank, npts_x_l % nprocs);
	end_i_l = (rank + 1) * (npts_x_l / nprocs) + fmin(rank + 1, npts_x_l % nprocs);
	start_i_j = rank * (npts_x_y / nprocs) + fmin(rank, npts_x_y % nprocs);
	end_i_j = (rank + 1) * (npts_x_y / nprocs) + fmin(rank + 1, npts_x_y % nprocs);
	
//	printf("Proc %d finished loading bessel & transfer\n", rank);

	MPI_Barrier(MPI_COMM_WORLD);
	
	// Precompute sin_tilde(x) and cos_tilde(x)
	precompute_tilde();
	printf("Proc %d Finished computing sin_tilde, cos_tilde.\n", rank);

//	printf("Proc %d finished precomputing\n", rank);

	MPI_Barrier(MPI_COMM_WORLD);

	// In order to put the data together, we define temporary arrays
	double **temp_cos_tilde = create_2D_array(npts_x, npts_l);
	double **temp_sin_tilde = create_2D_array(npts_x, npts_l);

	// MPI_Allreduce to incorporate calculated results
	MPI_Allreduce(&cos_tilde[0][0], &temp_cos_tilde[0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&sin_tilde[0][0], &temp_sin_tilde[0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	free_2D_array(cos_tilde);
	free_2D_array(sin_tilde);
	cos_tilde = temp_cos_tilde;
	sin_tilde = temp_sin_tilde;

	MPI_Barrier(MPI_COMM_WORLD);

	// We no longer need bessel and transfer functions data
	free_2D_array(bessel);
	free_2D_array(transfer);

	// Incorporate BN effects to the power spectrum and compute (2l+1)/Cl
	precompute_C_inv();
	
	free_array(C);
	free_array(beam);
	free_array(noise);
	
	// Use the Gauss-Legendre method for integrating mu
	gl_nodes = create_array(npts_mu);
	gl_weights = create_array(npts_mu);
	asy(gl_nodes, gl_weights, npts_mu); // Implemented in gl_integration.c
	
	// Values of Legendre polynomials at the nodes are computed and stored 
	legendre = create_2D_array(npts_mu, npts_l);
	#pragma omp parallel for private(mu,l)
	for(mu=0; mu<npts_mu; mu++){
		for(l=0; l<npts_l; l++){
			legendre[mu][l] = gsl_sf_legendre_Pl(lmin + l, gl_nodes[mu]);
		}
	}
	
//	printf("Proc %d finished preparation for the main computation\n", rank);
	
	// Main computation
	compute_error_bars();

//	printf("Proc %d finished the main computation\n", rank);

	MPI_Barrier(MPI_COMM_WORLD);

	// Again need temporary variables
	double final_N_cos, final_N_sin;

	// MPI Reduce N_cos & sins in root. Compute sigma in the root process
	MPI_Reduce(&N_cos, &final_N_cos, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&N_sin, &final_N_sin, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	printf("Proc #%d gets N_cos = %e\n, N_sin = %e\n", rank, N_cos, N_sin);

	if (rank == 0){
		N_cos = final_N_cos;
		N_sin = final_N_sin;
		sigma_cos = sqrt(6.0 / N_cos);
		sigma_sin = sqrt(6.0 / N_sin);
//		printf("Proc #%d has final value of N_cos = %e\n", rank, N_cos);
	}
	
	free_2D_array(cos_tilde);
	free_2D_array(sin_tilde);
	free_array(C_inv);
	free_2D_array(legendre);
	free_array(gl_nodes);
	free_array(gl_weights);
	free(i_j[0]);
	free(i_j);	

	// Print out the answer if root
	if (rank == 0){
		prompt();
	}

	MPI_Finalize();
	
	return 0;

}


int mpi_multiple_omega(int argc, char **argv){
	// Calculates error bars for multiple omegas using mpi

	// MPI parameters
	int rank, nprocs;

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Initialise
	parameter_init();

	// Distribute the workload by setting start and end points for vectorised indices.
	// Note that for non-mpi codes they are set to 0 and npts, respectively
	start_i_l = rank * (npts_x_l / nprocs) + fmin(rank, npts_x_l % nprocs);
	end_i_l = (rank + 1) * (npts_x_l / nprocs) + fmin(rank + 1, npts_x_l % nprocs);
	start_i_j = rank * (npts_x_y / nprocs) + fmin(rank, npts_x_y % nprocs);
	end_i_j = (rank + 1) * (npts_x_y / nprocs) + fmin(rank + 1, npts_x_y % nprocs);

	precompute_C_inv();
	free_array(C);
	free_array(beam);
	free_array(noise);

	// Use the Gauss-Legendre method for integrating mu
	gl_nodes = create_array(npts_mu);
	gl_weights = create_array(npts_mu);
	asy(gl_nodes, gl_weights, npts_mu); // Implemented in gl_integration.c

	// Values of Legendre polynomials at the nodes are computed and stored 
	legendre = create_2D_array(npts_mu, npts_l);
	#pragma omp parallel for private(mu,l)
	for(mu=0; mu<npts_mu; mu++){
		for(l=0; l<npts_l; l++){
			legendre[mu][l] = gsl_sf_legendre_Pl(lmin + l, gl_nodes[mu]);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Temporary variables for synchronising qtildes
	double **temp_cos_tilde = create_2D_array(npts_x, npts_l);
	double **temp_sin_tilde = create_2D_array(npts_x, npts_l);
	
	int npts_omega = 31;
	double **results = create_2D_array(npts_omega, 3);
	for(o=0; o<npts_omega; o++){
		omega = 10.0 * o;
		results[0][0] = omega;

		// Precompute cos_tilde, sin_tilde and synchronise
		precompute_tilde();

		MPI_Allreduce(&cos_tilde[0][0], &temp_cos_tilde[0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&sin_tilde[0][0], &temp_sin_tilde[0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		for(i=0; i<npts_x; i++){
			for(l=0; l<npts_l; l++){
				cos_tilde[i][l] = temp_cos_tilde[i][l];
				sin_tilde[i][l] = temp_sin_tilde[i][l];
			}
		}

		// Compute N_cos, N_sin and incorporate results
		compute_error_bars();

		double final_N_cos, final_N_sin;
		MPI_Reduce(&N_cos, &final_N_cos, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&N_sin, &final_N_sin, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if (rank == 0){
			N_cos = final_N_cos;
			N_sin = final_N_sin;
			sigma_cos = sqrt(6e0/N_cos);
			sigma_sin = sqrt(6e0/N_sin);
			results[o][1] = sigma_cos;
			results[o][2] = sigma_sin;
			prompt();
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (rank == 0){
		write_2D_array(results, npts_omega, 3, "multiple_omega_mpi_output"); 
	}

	MPI_Finalize();
	
	return 0;
}

int main(int argc, char **argv){
//	mpi_error_bars(argc, argv);
	mpi_multiple_omega(argc, argv);
//	error_bars();
	//bispectrum_plot();
	//multiple_omega();

	return 0;
}

