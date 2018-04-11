#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>
#include <string.h>

#include "gl_integration.h" // For Gauss-Legendre integration
#include "arrays.h"
#include "fixed_data.h"
#include "io.h"

static int lmin = 2;
static int lmax = 2000;
static int npts_l; 			//  = lmax - lmin + 1;
static int npts_mu; 			//  = 3*lmax/2 + 1;

static double *kvec;
static double *step_k;
static double kmin;
static double kmax = 8.466251e-1;	// this comes from kmax of the transfer functions data
static int npts_k;

static double *xvec;
static double *step_x;
static double xmax = 1.68e4;		// = 2 * tau0
static int npts_x;

static double omega = 1e2;

static double ***cos_tilde;		// cos_tilde[polarisation][x][l]
static double ***sin_tilde;
static int npts_p;			// Number of polarisation data included
static int npts_pp;			// e.g. 3 for TT, TE, EE

static double **C;

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


void make_xvec(){
	// Makes a vector of x using steps and cuts defined below.
//	double step3 = 200.0;
//	double step2 = 50.0;
//	double step1 = 10.0;
//	double step3 = 60.0;
//	double step2 = 30.0;
//	double step1 = 15.0;
//	double step3 = 30.0;
//	double step2 = 15.0;
//	double step1 = 10.0;
	double step3 = 6.0;
	double step2 = 3.0;
	double step1 = 1.0; 
	double cut1 = 13000.0;
	double cut2 = 13600.0;
	double cut3 = 14200.0;
	double cut4 = 15000.0;

	double x = 0.0;
	int i;
	for(i=0; x<xmax; i++){
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

	npts_x = i;
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

void prepare_C(){
	/* Initialises the covariance matrix and incorporates the beam and noise */

	// Load data. fixed_data.c
	load_C();
	load_BN();

	// Note that Cls starts with l=0
	C = create_2D_array(npts_pp, lmax);

	int pp, l;
	double *beam, *noise;
	for(pp=0; pp<npts_pp; pp++){

		C[pp] = get_C(pp);
		beam = get_beam(pp);
		noise = get_noise(pp);

		// Incorporate the beam and noise
		for(l=0; l<lmax; l++){
			C[pp][l] += noise[l] / (beam[l] * beam[l]); 
		}
	}

	// We no longer need the loaded data
	free_C();
	free_BN();

}

void print_result(){
	/* Prints parameters and computed error bars*/
	printf("Computing error bars for sinusodial shape functions\n");
	printf("*** Parameters ***\n");
	printf("omega = %e\n", omega);
	printf("lmin = %d, lmax = %d\n", lmin, lmax);
	printf("xmax = %e, npts_x = %d\n", xmax, npts_x);
	printf("kmax = %e, npts_k = %d\n", kmax, npts_k);
	printf("npts_mu = %d\n", npts_mu);
	printf("Amplitude factor deltaphi = %e\n", deltaphi);
	printf("Used %s, %s, %s\n", bessel_data_filename, transfer_T_data_filename, transfer_E_data_filename);
	printf("\n");
	printf("Results: N_cos = %e, N_sin = %e\n", N_cos, N_sin);
	printf("Results: sigma_cos = %e, sigma_sin = %e\n\n", sigma_cos, sigma_sin);
}


void initialise(){

	// Load data. fixed_data.c
	load_bessel();
	load_transfer();

	printf("loaded\n");
	// Initialise some parameters and the x, k vectors
	npts_l = lmax - lmin + 1;
	npts_mu = 3*lmax/2 + 1;

	kvec = transfer_kvec;
	npts_k = transfer_npts_k;
	kmin = kvec[0];
	kmax = kvec[npts_k-1];
	step_k = create_array(npts_k);

	int k;
	for(k=0; k<npts_k-1; k++) step_k[k] = kvec[k+1] - kvec[k];
	step_k[npts_k-1] = step_k[npts_k-2];

	printf("kvec made\n");
	make_xvec();

	// Vectorise indicies for more efficient looping. arrays.c
	i_l = vectorise_two_indices(npts_x, npts_l);
	npts_x_l = npts_x * npts_l;
	start_i_l = 0;
	end_i_l = npts_x_l;

	i_j = vectorise_two_indices(npts_x, npts_x);
	npts_x_y = npts_x * npts_x;
	start_i_j = 0;
	end_i_j = npts_x_y;

	// Allocate memory for tilde arrays

	if(do_polarisation == 0){
		cos_tilde = (double ***) malloc(sizeof(double **));
		cos_tilde[T] = create_2D_array(npts_x, npts_l);
		sin_tilde = (double ***) malloc(sizeof(double **));
		sin_tilde[T] = create_2D_array(npts_x, npts_l);

		npts_p = 1;
		npts_pp = 1; 	// Only C_TT used

	}else if(do_polarisation == 1){
		cos_tilde = (double ***) malloc(2 * sizeof(double **));
		cos_tilde[T] = create_2D_array(npts_x, npts_l);
		cos_tilde[E] = create_2D_array(npts_x, npts_l);
		sin_tilde = (double ***) malloc(2 * sizeof(double **));
		sin_tilde[T] = create_2D_array(npts_x, npts_l);
		sin_tilde[E] = create_2D_array(npts_x, npts_l);

		npts_p = 2;
		npts_pp = 3;	// C_TT, C_TE, C_EE used
	}		

	// Initialise the covariance matrix
	prepare_C();

}


void precompute_tilde(){
	/* Evaluate sin_tilde(x,l) and cos_tilde(x,l) */

	// Initialise to zero
	int i, l, p;
	for(p=0; p<npts_p; p++){
		for(i=0; i<npts_x; i++){
			for(l=0; l<npts_l; l++){
					cos_tilde[p][i][l] = sin_tilde[p][i][l] = 0;
			}
		}
	}

	// Perform integral over k to evaluate cos_tilde(x,l) and  sin_tilde(x,l)
	// cos_tilde(x,l) = (2/pi) * integral{dk * cos(omega*k) * j_l(k*x) * Delta_l(k)}
	// sin_tilde(x,l) = (2/pi) * integral{dk * sin(omega*k) * j_l(k*x) * Delta_l(k)}

	double pref = 2e0/pi;

	#pragma omp parallel private(i,l)
	{
		// Initialise GSL tools for interpolation and integration
		gsl_spline *bessel_spline = gsl_spline_alloc(gsl_interp_linear, bessel_npts_x);
		gsl_interp_accel *bessel_acc = gsl_interp_accel_alloc();
		double bes;	// temporary variables
		double **trans = (double **) malloc(npts_p * sizeof(double *));
		int n, k, p;

		#pragma omp for
		for(n=start_i_l; n<end_i_l; n++){
			i = i_l[n][0];
			l = i_l[n][1];

			// Prepare bessel and transfer data vectors
			gsl_spline_init(bessel_spline, bessel_xvec, get_bessel(l+lmin), bessel_npts_x);
			for(p=0; p<npts_p; p++){
				trans[p] = get_transfer(p, l+lmin);
			}

			// Perform the k integral
			for(k=0; k<npts_k; k++){
				bes = gsl_spline_eval(bessel_spline, xvec[i] * kvec[k], bessel_acc);

				for(p=0; p<npts_p; p++){
					cos_tilde[p][i][l] += step_k[k] * cos(omega * kvec[k]) * trans[p][k] * bes; 
					sin_tilde[p][i][l] += step_k[k] * sin(omega * kvec[k]) * trans[p][k] * bes; 
				}
			}

			// Multiply by 2/pi
			for(p=0; p<npts_p; p++){
				cos_tilde[p][i][l] *= pref;
				sin_tilde[p][i][l] *= pref;
			}
			
			gsl_interp_accel_reset(bessel_acc);
		}

		// Free the interpolaters
		gsl_spline_free(bessel_spline);
		gsl_interp_accel_free(bessel_acc);
		free(trans);
	}
}

void orthogonalise_tilde(){
	/* Orthogonalises the computed cos_tilde, sin_tilde so that C = I */

	int x, l;

	// Precompute 1/sqrt(C_TT) for computational efficiency
	double *vec = create_array(npts_l);
	for(l=0; l<npts_l; l++){
		vec[l] = 1e0 / sqrt(C[TT][l+lmin]);
	}

	#pragma omp parallel for private(x,l)
	for(x=0; x<npts_x; x++){
		for(l=0; l<npts_l; l++){
			cos_tilde[T][x][l] *= vec[l];
			sin_tilde[T][x][l] *= vec[l];
		}
	}

	if(do_polarisation == 1){
		// Again precompute the denominator
		for(l=0; l<npts_l; l++){
			vec[l] *= 1e0 / sqrt(C[TT][l+lmin] * C[EE][l+lmin] - C[TE][l+lmin] * C[TE][l+lmin]);
		}

		#pragma omp parallel for private(x,l)
		for(x=0; x<npts_x; x++){
			for(l=0; l<npts_l; l++){
				cos_tilde[E][x][l] = vec[l] * (C[TT][l+lmin] * cos_tilde[E][x][l] - C[TE][l+lmin] * cos_tilde[T][x][l]);
				sin_tilde[E][x][l] = vec[l] * (C[TT][l+lmin] * sin_tilde[E][x][l] - C[TE][l+lmin] * sin_tilde[T][x][l]);
			}
		}
	}

	free_array(vec);

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

	double pref = 1e0/(8e0*pi) * pow(6e0 * deltaphi * deltaphi, 2);

	#pragma omp parallel
	{
		double cc, cs, sc, ss;
		double dxdy, x2y2;
		int i, j, n, l, mu;

		#pragma omp for reduction(+:N_cos,N_sin)
		for(n=start_i_j; n<end_i_j; n++){
			i = i_j[n][0];
			j = i_j[n][1];

			dxdy = step_x[i] * step_x[j];
			x2y2 = (xvec[i]*xvec[i]) * (xvec[j]*xvec[j]);
			
			for(mu=0; mu<npts_mu; mu++){
				// Sum over l first
				cc = cs = sc = ss = 0;
				for(l=0; l<npts_l; l++){
					cc += C_inv[l] * cos_tilde[T][i][l] * cos_tilde[T][j][l] * legendre[mu][l];
					cs += C_inv[l] * cos_tilde[T][i][l] * sin_tilde[T][j][l] * legendre[mu][l];
					sc += C_inv[l] * sin_tilde[T][i][l] * cos_tilde[T][j][l] * legendre[mu][l];
					ss += C_inv[l] * sin_tilde[T][i][l] * sin_tilde[T][j][l] * legendre[mu][l];
				}	

				N_cos += gl_weights[mu] * dxdy * x2y2 * ((cc * cc * cc) + 3*(cc * ss * ss) - 3*(sc * sc * cc) - 3*(cs * cs * cc) + 6*(sc * cs * ss));
				N_sin += gl_weights[mu] * dxdy * x2y2 * ((ss * ss * ss) + 3*(ss * cc * cc) - 3*(cs * cs * ss) - 3*(sc * sc * ss) + 6*(cs * sc * cc));
			}
		}
	}

	N_cos *= pref;
	N_sin *= pref;
	sigma_cos = sqrt(6e0/N_cos) * sqrt(1e0/fskyT);
	sigma_sin = sqrt(6e0/N_sin) * sqrt(1e0/fskyT);

}

void multiple_omega(){
	/* Calculates error bars for multiple omegas */

	// Initialise
	initialise();

	precompute_C_inv();
	free_C();
	free_BN();

	// Use the Gauss-Legendre method for integrating mu
	gl_nodes = create_array(npts_mu);
	gl_weights = create_array(npts_mu);
	asy(gl_nodes, gl_weights, npts_mu); // Implemented in gl_integration.c

	// Values of Legendre polynomials at the nodes are computed and stored 
	legendre = create_2D_array(npts_mu, npts_l);

	int mu, l;
	#pragma omp parallel for private(mu,l)
	for(mu=0; mu<npts_mu; mu++){
		for(l=0; l<npts_l; l++){
			legendre[mu][l] = gsl_sf_legendre_Pl(lmin + l, gl_nodes[mu]);
		}
	}


	int o, npts_omega = 10;
	double results[npts_omega][3];
	for(o=0; o<npts_omega; o++){
		omega = 10.0 * (o + 1);
		results[o][0] = omega;

		precompute_tilde();
		compute_error_bars();

		results[o][1] = sigma_cos;
		results[o][2] = sigma_sin;
		print_result();
	}

	printf("\n\n sigma_cos: ");
	for(o=0; o<npts_omega; o++) printf("%e ", results[o][1]);
	printf("\n sigma_sin: ");
	for(o=0; o<npts_omega; o++) printf("%e ", results[o][2]);
	printf("\n");
	
}

void bispectrum_plot(){
	// Calculates bispectrum for l1=l2=l3=l and saves the results to a file

	// Initialise
	initialise();

	precompute_tilde();

	free_bessel();
	free_transfer();

	double *b_cos = (double *) calloc(npts_l, sizeof(double));
	double *b_sin = (double *) calloc(npts_l, sizeof(double));
	
	double c, s;
	int l, i;
	#pragma omp parallel for private(l,i,c,s)
	for(l=0; l<npts_l; l++){
		for(i=0; i<npts_x; i++){
			c = cos_tilde[T][i][l];
			s = sin_tilde[T][i][l];
			b_cos[l] += xvec[i] * xvec[i] * step_x[i] * (c*c*c - 3*c*s*s);
			b_sin[l] += xvec[i] * xvec[i] * step_x[i] * (-s*s*s + 3*s*c*c);
		}
	}

	write_1D_array(b_cos, npts_l, "cos_bispectrum_equal_l");
	write_1D_array(b_sin, npts_l, "sin_bispectrum_equal_l");
}

void const_model(){
	// Computes the bispectrum for constant shape function S to compare with analytic solution

	lmax = 200;
	initialise();
	double **transfer = create_2D_array(npts_l, npts_k);

	// Large angle approximation for l<<200. Transfer(l,k) ~= (1/3)*bessel(l,tau0*k)
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, bessel_npts_x);
	int l, k;
	for(l=0; l<npts_l; l++){
		gsl_spline_init(spline, bessel_xvec, get_bessel(l+lmin), bessel_npts_x);
		for(k=0; k<transfer_npts_k; k++){
			transfer[l][k] = (1e0/3e0) * gsl_spline_eval(spline, transfer_kvec[k] * tau0, acc);
		}
		gsl_interp_accel_reset(acc);	
	}

	// S=cos(wk) for omega = 0 corresponds to the constant model
	omega = 0;

	// **** THIS NEEDS TO BE CHANGED ****
	precompute_tilde();	// Here large angle approx. to the transfer functions needs to be used
	write_2D_array(cos_tilde[T], npts_x, npts_l, "const_tilde");


	// Bispectrum
	double ***b_const = create_3D_array(npts_l, npts_l, npts_l);
	double ***b_analytic = create_3D_array(npts_l, npts_l, npts_l);
	double threshold = 1e-3;

	double *xint = create_array(npts_x);
	int i;
	for(i=0; i<npts_x; i++){
		xint[i] = step_x[i] * xvec[i] * xvec[i];	// to save computing time later
	}

	#pragma omp parallel private(i)
	{
		int l1, l2, l3;
		int ii, jj, kk;
		double error;

		#pragma omp for
		for(ii=0; ii<npts_l; ii++){
			l1 = ii + lmin;

			for(jj=0; jj<=ii; jj++){
				l2 = jj + lmin;

				// triangle condition and l1>=l2>=l3 enforced here
				for(kk=abs(l1-l2)-lmin; kk<=jj; kk++){
					l3 = kk + lmin;

					b_const[ii][jj][kk] = 0;
					for(i=0; i<npts_x; i++){
						b_const[ii][jj][kk] += xint[i] * cos_tilde[T][i][ii] * cos_tilde[T][i][jj] * cos_tilde[T][i][kk];
					}

					b_analytic[ii][jj][kk] = 1e0 / (27e0 * (2*l1+1) * (2*l2+1) * (2*l3+1)) * (1e0 / (l1+l2+l3+3e0) + 1e0 / (l1+l2+l3));

					error = (b_const[ii][jj][kk] - b_analytic[ii][jj][kk]) / b_analytic[ii][jj][kk];

					if(fabs(error) > threshold){
						printf("b(%d,%d,%d) = %e, analytic %e, error %e\n", l1, l2, l3, b_const[ii][jj][kk], b_analytic[ii][jj][kk], error);
					}
				}
			}
		}
		
	}

}


void error_bars(){

	clock_t start = clock();
	
	// Initialise and load data
	initialise();
	
	// Precompute sin_tilde(x) and cos_tilde(x)
	precompute_tilde();
	
	free_bessel();
	free_transfer();

	// Incorporate BN effects to the power spectrum and compute (2l+1)/Cl
	precompute_C_inv();
	printf("Finished computing sin_tilde, cos_tilde. Elapsed time: %fs\n", (double)(clock() - start) / CLOCKS_PER_SEC);
	
	free_C();
	free_BN();
	
	// Use the Gauss-Legendre method for integrating mu
	gl_nodes = create_array(npts_mu);
	gl_weights = create_array(npts_mu);
	asy(gl_nodes, gl_weights, npts_mu); // Implemented in gl_integration.c
	
	// Values of Legendre polynomials at the nodes are computed and stored 
	legendre = create_2D_array(npts_mu, npts_l);

	int mu, l;
	#pragma omp parallel for private(mu,l)
	for(mu=0; mu<npts_mu; mu++){
		for(l=0; l<npts_l; l++){
			legendre[mu][l] = gsl_sf_legendre_Pl(lmin + l, gl_nodes[mu]);
		}
	}
	
	
	// Main computation
	compute_error_bars();
	
	free_2D_array(cos_tilde[T]);
	free_2D_array(sin_tilde[T]);
	free_array(C_inv);
	free_2D_array(legendre);
	free_array(gl_nodes);
	free_array(gl_weights);
	free(i_j[0]);
	free(i_j);	

	// Print out the answer
	print_result();
	
	clock_t end = clock();
	double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Elapsed time: %fs\n", time_spent);
	
}

void write_qtilde_integrand(){

	initialise();

	/* Save integrands for k the integration */

	int npts_l_save = 10;

	double ***cos_integrand = create_3D_array(npts_l_save, npts_x, npts_k);
	double ***sin_integrand = create_3D_array(npts_l_save, npts_x, npts_k);

	free_2D_array(cos_tilde[T]);
	free_2D_array(sin_tilde[T]);

	cos_tilde[T] = create_2D_array(npts_x, npts_l_save);
	sin_tilde[T] = create_2D_array(npts_x, npts_l_save);

	// Initialise to zero
	int i, l;
	for(i=0; i<npts_x; i++){
		for(l=0; l<npts_l_save; l++){
			cos_tilde[T][i][l] = sin_tilde[T][i][l] = 0;
		}
	}

	// Perform integral over k to evaluate cos_tilde(x,l) and  sin_tilde(x,l)
	// cos_tilde(x,l) = (2/pi) * integral{dk * cos(omega*k) * j_l(k*x) * Delta_l(k)}
	// sin_tilde(x,l) = (2/pi) * integral{dk * sin(omega*k) * j_l(k*x) * Delta_l(k)}

	#pragma omp parallel private(i,l)
	{
		// Initialise GSL tools for interpolation and integration
		gsl_spline *bessel_spline = gsl_spline_alloc(gsl_interp_linear, bessel_npts_x);
		gsl_interp_accel *bessel_acc = gsl_interp_accel_alloc();
		gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, npts_k);
		gsl_interp_accel *acc = gsl_interp_accel_alloc();
		double bes;	// temporary variables
		double *transfer;
		double *y1 = create_array(npts_k);
		double *y2 = create_array(npts_k);
		int k;

		#pragma omp for
		for(l=0; l<npts_l_save; l++){
			for(i=0; i<npts_x; i++){

				gsl_spline_init(bessel_spline, bessel_xvec, get_bessel(l+lmin), bessel_npts_x);
				transfer = get_transfer(l+lmin);

				for(k=0; k<npts_k; k++){
					bes = gsl_spline_eval(bessel_spline, xvec[i] * kvec[k], bessel_acc);

					y1[k] = cos(omega * kvec[k]) * transfer[k] * bes;
					y2[k] = sin(omega * kvec[k]) * transfer[k] * bes;

					cos_integrand[l][i][k] = y1[k];
					sin_integrand[l][i][k] = y2[k];
				}

				gsl_interp_accel_reset(bessel_acc);

				gsl_spline_init(spline, kvec, y1, npts_k);
				cos_tilde[T][i][l] = (2e0/pi) * gsl_spline_eval_integ(spline, kmin, kmax, acc);
				gsl_interp_accel_reset(acc);

				gsl_spline_init(spline, kvec, y2, npts_k);
				sin_tilde[T][i][l] = (2e0/pi) * gsl_spline_eval_integ(spline, kmin, kmax, acc);
				gsl_interp_accel_reset(acc);
			}
		}

		// Free the interpolaters
		gsl_spline_free(bessel_spline);
		gsl_interp_accel_free(bessel_acc);
		gsl_spline_free(spline);
		gsl_interp_accel_free(acc);
	}

	// Write data
	for(l=0; l<npts_l_save; l++){
		char name1[200], name2[200];
		sprintf(name1, "cos_tilde_integrand_l_%d", l+lmin);
		sprintf(name2, "sin_tilde_integrand_l_%d", l+lmin);
		write_2D_array(cos_integrand[l], npts_x, npts_k, name1);
		write_2D_array(sin_integrand[l], npts_x, npts_k, name2);
	}	
}


void mpi_error_bars(int argc, char **argv){

	// MPI parameters
	int rank, nprocs;

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int i, mu, l; //loop indices

	// Initialise parameters and load data
	initialise();

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
	double **temp_cos_tilde[T] = create_2D_array(npts_x, npts_l);
	double **temp_sin_tilde[T] = create_2D_array(npts_x, npts_l);

	// MPI_Allreduce to incorporate calculated results
	MPI_Allreduce(&cos_tilde[T][0][0], &temp_cos_tilde[T][0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&sin_tilde[T][0][0], &temp_sin_tilde[T][0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	for(i=0; i<npts_x; i++){
		for(l=0; l<npts_l; l++){
			cos_tilde[T][i][l] = temp_cos_tilde[T][i][l];
			sin_tilde[T][i][l] = temp_sin_tilde[T][i][l];
		}
	}
	free_2D_array(temp_cos_tilde[T]);
	free_2D_array(temp_sin_tilde[T]);

	MPI_Barrier(MPI_COMM_WORLD);

	// We no longer need bessel and transfer functions data
	free_bessel();
	free_transfer();

	// Incorporate BN effects to the power spectrum and compute (2l+1)/Cl
	precompute_C_inv();
	
	free_C();
	free_BN();
	
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
	printf("Proc #%d gets N_cos = %e\n, N_sin = %e\n", rank, N_cos, N_sin);

	MPI_Barrier(MPI_COMM_WORLD);

	// Again need temporary variables
	double final_N_cos, final_N_sin;

	// MPI Reduce N_cos & sins in root. Compute sigma in the root process
	MPI_Reduce(&N_cos, &final_N_cos, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&N_sin, &final_N_sin, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


	if (rank == 0){
		N_cos = final_N_cos;
		N_sin = final_N_sin;
		sigma_cos = sqrt(6e0 / N_cos) * sqrt(1e0/fskyT);
		sigma_sin = sqrt(6e0 / N_sin) * sqrt(1e0/fskyT);

		// Print out the answer
		print_result();

		// Some intermediate data saving 
		write_2D_array(cos_tilde[T], npts_x, npts_l, "cos_tilde_last");
		write_2D_array(sin_tilde[T], npts_x, npts_l, "sin_tilde_last"); 	
		write_1D_array(xvec, npts_x, "xvec_last");


	}
	
	free_2D_array(cos_tilde[T]);
	free_2D_array(sin_tilde[T]);
	free_array(C_inv);
	free_2D_array(legendre);
	free_array(gl_nodes);
	free_array(gl_weights);
	free(i_j[0]);
	free(i_j);	

	MPI_Finalize();

}


void mpi_multiple_omega(int argc, char **argv){
	// Calculates error bars for multiple omegas using mpi

	// MPI parameters
	int rank, nprocs;

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int i, mu, l; //loop indices

	// Initialise
	initialise();

	// Distribute the workload by setting start and end points for vectorised indices.
	// Note that for non-mpi codes they are set to 0 and npts, respectively
	start_i_l = rank * (npts_x_l / nprocs) + fmin(rank, npts_x_l % nprocs);
	end_i_l = (rank + 1) * (npts_x_l / nprocs) + fmin(rank + 1, npts_x_l % nprocs);
	start_i_j = rank * (npts_x_y / nprocs) + fmin(rank, npts_x_y % nprocs);
	end_i_j = (rank + 1) * (npts_x_y / nprocs) + fmin(rank + 1, npts_x_y % nprocs);

	precompute_C_inv();
	free_C();
	free_BN();

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
	double **temp_cos_tilde[T] = create_2D_array(npts_x, npts_l);
	double **temp_sin_tilde[T] = create_2D_array(npts_x, npts_l);
	
	int o, npts_omega = 41;
	double **results = create_2D_array(npts_omega, 3);
	for(o=0; o<npts_omega; o++){
		omega = 10.0 * o;
		results[o][0] = omega;

		// Precompute cos_tilde, sin_tilde and synchronise
		precompute_tilde();

		MPI_Allreduce(&cos_tilde[T][0][0], &temp_cos_tilde[T][0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&sin_tilde[T][0][0], &temp_sin_tilde[T][0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		for(i=0; i<npts_x; i++){
			for(l=0; l<npts_l; l++){
				cos_tilde[T][i][l] = temp_cos_tilde[T][i][l];
				sin_tilde[T][i][l] = temp_sin_tilde[T][i][l];
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
			sigma_cos = sqrt(6e0/N_cos) * sqrt(1e0/fskyT);
			sigma_sin = sqrt(6e0/N_sin) * sqrt(1e0/fskyT);
			results[o][1] = sigma_cos;
			results[o][2] = sigma_sin;
			print_result();
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (rank == 0){
		write_2D_array(results, npts_omega, 3, "multiple_omega_mpi_lensed_output"); 
	}

	MPI_Finalize();
}

void mpi_bispectrum(int argc, char **argv){

	// MPI parameters
	int rank, nprocs;

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Initialise parameters and load data
	initialise();


//	printf("Proc %d finished initialising\n", rank);

	// Distribute the workload by setting start and end points for vectorised indices.
	// Note that for non-mpi codes they are set to 0 and npts, respectively
	start_i_l = rank * (npts_x_l / nprocs) + fmin(rank, npts_x_l % nprocs);
	end_i_l = (rank + 1) * (npts_x_l / nprocs) + fmin(rank + 1, npts_x_l % nprocs);
	
//	printf("Proc %d finished loading bessel & transfer\n", rank);

	MPI_Barrier(MPI_COMM_WORLD);
	
	// Precompute sin_tilde(x) and cos_tilde(x)
	precompute_tilde();
	printf("Proc %d Finished computing sin_tilde, cos_tilde.\n", rank);

//	printf("Proc %d finished precomputing\n", rank);

	MPI_Barrier(MPI_COMM_WORLD);

	// In order to put the data together, we define temporary arrays
	double **temp_cos_tilde[T] = create_2D_array(npts_x, npts_l);
	double **temp_sin_tilde[T] = create_2D_array(npts_x, npts_l);

	// MPI_Allreduce to incorporate calculated results
	MPI_Allreduce(&cos_tilde[T][0][0], &temp_cos_tilde[T][0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&sin_tilde[T][0][0], &temp_sin_tilde[T][0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	int i, l;	
	for(i=0; i<npts_x; i++){
		for(l=0; l<npts_l; l++){
			cos_tilde[T][i][l] = temp_cos_tilde[T][i][l];
			sin_tilde[T][i][l] = temp_sin_tilde[T][i][l];
		}
	}
	free_2D_array(temp_cos_tilde[T]);
	free_2D_array(temp_sin_tilde[T]);

	MPI_Barrier(MPI_COMM_WORLD);

	// We no longer need bessel and transfer functions data
	free_bessel();
	free_transfer();

	// Bispectrum
	int lstep = 10;
	npts_l = npts_l / lstep;
	lmax = lmin + (npts_l - 1) * lstep;
	double ***bispectrum = create_3D_array(npts_l, npts_l, npts_l);

	// Initialise to zero
	int l1, l2, l3;
	for(l1=0; l1<npts_l; l1++){
		for(l2=0; l2<npts_l; l2++){
			for(l3=0; l3<npts_l; l3++){
				bispectrum[l1][l2][l3]= 0;
			}
		}
	}
			
	// Distribute workload
	int **l_l = vectorise_two_indices(npts_l, npts_l);	
	int npts_l_l = npts_l * npts_l;
	int start_l_l, end_l_l;
	start_l_l = rank * (npts_l_l / nprocs) + fmin(rank, npts_l_l % nprocs);
	end_l_l = (rank + 1) * (npts_l_l / nprocs) + fmin(rank + 1, npts_l_l % nprocs);


	#pragma omp parallel private(i,l1,l2,l3)
	{
		gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, npts_x);
		gsl_interp_accel *acc = gsl_interp_accel_alloc();
		double *y = create_array(npts_x);
		int n;

		#pragma omp for
		for(n=start_l_l; n<end_l_l; n++){
			l1 = l_l[n][0];
			l2 = l_l[n][1];

			for (l3=0; l3<npts_l; l3++){
				if (lstep*(l1+l2)+lmin >= lstep*l3 && lstep*(l2+l3)+lmin >= lstep*l1 && lstep*(l3+l1)+lmin >= lstep*l2){ 
					for(i=0; i<npts_x; i++){
						y[i] = xvec[i] * xvec[i] * (cos_tilde[T][i][lstep*l1]*cos_tilde[T][i][lstep*l2]*cos_tilde[T][i][lstep*l3] - cos_tilde[T][i][lstep*l1]*sin_tilde[T][i][lstep*l2]*sin_tilde[T][i][lstep*l3] - sin_tilde[T][i][lstep*l1]*cos_tilde[T][i][lstep*l2]*sin_tilde[T][i][lstep*l3] - sin_tilde[T][i][lstep*l1]*sin_tilde[T][i][lstep*l2]*cos_tilde[T][i][lstep*l3]);

					}

					gsl_spline_init(spline, xvec, y, npts_x);
					bispectrum[l1][l2][l3] = gsl_spline_eval_integ(spline, xvec[0], xvec[npts_x-1], acc);	 
					gsl_interp_accel_reset(acc);
				}
			}
		}

		gsl_spline_free(spline);
		gsl_interp_accel_free(acc);
		free_array(y);

	}


	MPI_Barrier(MPI_COMM_WORLD);

	// Again need temporary variables
	double ***final_bispectrum = create_3D_array(npts_l, npts_l, npts_l);

	// MPI reduce the bisepctrum in root. 
	MPI_Reduce(&bispectrum[0][0][0], &final_bispectrum[0][0][0], npts_l*npts_l*npts_l, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0){
		char *filename = "cos_bispectrum_3D_2000";
		write_bispectrum(final_bispectrum, npts_l, lstep, lstep, filename);

		printf("Computing bispectrum for feature models\n");
		printf("*** Parameters ***\n");
		printf("omega = %e\n", omega);
		printf("lmin = %d, lmax = %d\n", lmin, lmax);
		printf("xmax = %e, npts_x = %d\n", xmax, npts_x);
		printf("kmax = %e, npts_k = %d\n", kmax, npts_k);
		printf("Used bessel_4000 data set and transfer functions,  power spectrum from Planck\n");
		printf("\n");
		printf("Results saved to file %s\n", filename);
	}
	
	free_2D_array(cos_tilde[T]);
	free_2D_array(sin_tilde[T]);


	MPI_Finalize();
	

}

int main(int argc, char **argv){
//	mpi_error_bars(argc, argv);
	mpi_multiple_omega(argc, argv);
//	error_bars();
	//bispectrum_plot();
	//multiple_omega();
//	mpi_bispectrum(argc, argv);
//	const_model();	
//	write_qtilde_integrand();

	return 0;
}

