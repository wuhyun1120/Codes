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
static double step_o = 10.0;
static int npts_o = 31;

static double ***cos_tilde;		// cos_tilde[polarisation][x][l]
static double ***sin_tilde;
static int npts_p;			// Number of polarisation data included
static int npts_pp;			// e.g. 3 for TT, TE, EE
static double fsky;

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

static double N_cross;
static double sigma_cross;

int min(int x, int y){
	 return (x < y) ? x : y;
}
int max(int x, int y){
	return (x > y) ? x : y;
}


void make_xvec(){
	///* Makes a vector of x which is denser around the last scattering surface ~13900 */  
	
//	double step3 = 200.0;
//	double step2 = 50.0;
//	double step1 = 10.0;
	double step3 = 60.0;
	double step2 = 30.0;
	double step1 = 15.0;
//	double step3 = 30.0;
//	double step2 = 15.0;
//	double step1 = 10.0;
//	double step3 = 6.0;
//	double step2 = 3.0;
//	double step1 = 1.0; 
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
	///* Initialises the covariance matrix and incorporates the beam and noise */

	// Load data. fixed_data.c
	load_C();
	load_BN();

	C = create_2D_array(npts_pp, npts_l);

	int pp, l;
	double *raw_C, *beam, *noise;
	for(pp=0; pp<npts_pp; pp++){

		// Get deep copies of data. fixed_data.c
		raw_C = get_C(pp);
		beam = get_beam(pp);
		noise = get_noise(pp);

		// Incorporate the beam and noise
		for(l=0; l<npts_l; l++){
			// Note the loaded raw_C, beam and noise all start with l=0
			C[pp][l] = raw_C[l+lmin] + noise[l+lmin] / (beam[l+lmin] * beam[l+lmin]); 
		}

		free(raw_C);
		free(beam);
		free(noise);
	}

	// We no longer need the loaded data
	free_C();
	free_BN();

}

void initialise(){

	// Load data. fixed_data.c
	load_bessel();
	load_transfer();

//	printf("loaded\n");
	// Initialise some parameters and the x, k vectors
	npts_l = lmax - lmin + 1;

	kvec = transfer_kvec;
	npts_k = transfer_npts_k;
	kmin = kvec[0];
	kmax = kvec[npts_k-1];
	step_k = create_array(npts_k);

	int k;
	for(k=0; k<npts_k-1; k++) step_k[k] = kvec[k+1] - kvec[k];
	step_k[npts_k-1] = step_k[npts_k-2];

//	printf("kvec made\n");
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
		fsky = fskyT;

	}else if(do_polarisation == 1){
		cos_tilde = (double ***) malloc(2 * sizeof(double **));
		cos_tilde[T] = create_2D_array(npts_x, npts_l);
		cos_tilde[E] = create_2D_array(npts_x, npts_l);
		sin_tilde = (double ***) malloc(2 * sizeof(double **));
		sin_tilde[T] = create_2D_array(npts_x, npts_l);
		sin_tilde[E] = create_2D_array(npts_x, npts_l);

		npts_p = 2;
		npts_pp = 3;	// C_TT, C_TE, C_EE used
		fsky = fskyE;
	}		

	// Initialise the covariance matrix
	prepare_C();

}

void precompute_tilde(){
	///* Evaluate sin_tilde(x,l) and cos_tilde(x,l) */

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

	#pragma omp parallel private(i,l,p)
	{
		// Initialise GSL tools for interpolation and integration
		gsl_spline *bessel_spline = gsl_spline_alloc(gsl_interp_linear, bessel_npts_x);
		gsl_interp_accel *bessel_acc = gsl_interp_accel_alloc();
		double bes;	// temporary variables
		double **trans = (double **) malloc(npts_p * sizeof(double *));
		int n, k;

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
	// Linear transformation 'L' on polarisation indices (T,E):
	// cos_tilde[T] -> (1 / sqrt(C[TT})) * cos_tilde[T],
	// cos_tilde[E] -> (C[TT] * cos_tilde[E] - C[TE] * cos_tilde[T]) / sqrt(C[TT] * (C[TT] * C[EE] - C[TE] * C[TE])

	if(do_polarisation == 0){

		int i, l, n;
		double *L_TT = create_array(npts_l);
		for(l=0; l<npts_l; l++){
			L_TT[l] = sqrt(1e0 / C[TT][l]);
		}

		#pragma omp parallel for private(i,l,n)
		for(n=start_i_l; n<end_i_l; n++){
			i = i_l[n][0];
			l = i_l[n][1];
			cos_tilde[T][i][l] *= L_TT[l];
			sin_tilde[T][i][l] *= L_TT[l];
		}
	
		free_array(L_TT);

	}else if(do_polarisation == 1){

		int i, l, n;
		double *L_TT = create_array(npts_l);
		double *L_ET = create_array(npts_l);
		double *L_EE = create_array(npts_l);
	
		for(l=0; l<npts_l; l++){
			L_TT[l] = sqrt(1e0 / C[TT][l]);
			L_ET[l] = - C[TE][l] * L_TT[l] / sqrt(C[TT][l] * C[EE][l] - C[TE][l] * C[TE][l]);
			L_EE[l] = C[TT][l] * L_TT[l] / sqrt(C[TT][l] * C[EE][l] - C[TE][l] * C[TE][l]);
		}

		#pragma omp parallel for private(i,l,n)
		for(n=start_i_l; n<end_i_l; n++){
			i = i_l[n][0];
			l = i_l[n][1];
			cos_tilde[E][i][l] = L_ET[l] * cos_tilde[T][i][l] + L_EE[l] * cos_tilde[E][i][l];
			sin_tilde[E][i][l] = L_ET[l] * sin_tilde[T][i][l] + L_EE[l] * sin_tilde[E][i][l];
			cos_tilde[T][i][l] *= L_TT[l];
			sin_tilde[T][i][l] *= L_TT[l];
		}

		free_array(L_TT);
		free_array(L_ET);
		free_array(L_EE);
	}
}


void prepare_mu_integration(){
	///* Prepares Gauss-Legendre integration and precomputes the Legendre function values needed*/

	// Gauss-Legendre quadrature is precise for polynomials of degree up to 2n-1.
	// Since we have three Pl's on our integrand, we require (3*lmax+1)/2 points
	npts_mu = 3*lmax/2 + 1;
	gl_nodes = create_array(npts_mu);
	gl_weights = create_array(npts_mu);
	asy(gl_nodes, gl_weights, npts_mu); // Implemented in gl_integration.c

	// Values of Legendre polynomials at the nodes are computed and stored 
	legendre = create_2D_array(npts_mu, npts_l);
	int mu, l;
	#pragma omp parallel for private(mu,l)
	for(mu=0; mu<npts_mu; mu++){
		for(l=0; l<npts_l; l++){
			// A factor of 2l+1 was included here for later convenience
			legendre[mu][l] = (2e0*(l+lmin)+1e0) * gsl_sf_legendre_Pl(l+lmin, gl_nodes[mu]);
		}
	}
}


void compute_error_bars(){
	///* Compute P_ss, P_sc, P_cs, P_cc and use them to compute N and hence sigma */

	// Compute N
	// For S = sin,
	// N = (1/8pi) * integral{x^2*dx * y^2*dy * dmu * [3*P_ss*P_cc*P_cc + 6*P_sc*P_cs*P_cc + 1*P_ss*P_ss*P_ss - 3*P_sc*P_sc*P_ss - 3*P_cs*P_cs*P_ss]}
	// For S = cos,
	// N = (1/8pi) * integral{x^2*dx * y^2*dy * dmu * [3*P_cc*P_ss*P_ss + 6*P_cs*P_sc*P_ss + 1*P_cc*P_cc*P_cc - 3*P_cs*P_cs*P_cc - 3*P_sc*P_sc*P_cc]}
	// and finally for the cross term between sin and cos,
	// N = (1/8pi) * integral{x^2*dx * y^2*dy * dmu * [3*P_cs*P_cc*P_cc - P_cs*P_cs*P_cs - 3*P_cs*P_sc*P_sc - 6*P_cc*P_ss*P_sc + 3*P_cs*P_ss*P_ss]}
	// where
	// P_ss(x,y,mu) = sum_l{((2l+1)/C_l) * sin_tilde(x,l) * sin_tilde(y,l) * P_l(mu)}
	// P_cs(x,y,mu) = sum_l{((2l+1)/C_l) * cos_tilde(x,l) * sin_tilde(y,l) * P_l(mu)}
	// P_sc(x,y,mu) = sum_l{((2l+1)/C_l) * sin_tilde(x,l) * cos_tilde(y,l) * P_l(mu)}
	// P_cc(x,y,mu) = sum_l{((2l+1)/C_l) * cos_tilde(x,l) * cos_tilde(y,l) * P_l(mu)}
	// With polarisation data, use P_ss = P_ss[T] + P_ss[E] etc.

	N_cos = N_sin = N_cross = 0;

	// Factor of 1/(8*pi) from the integral
	// Factor of 6*(delta_phi)^2 from the definition of shape function for the feature model
	// ... and a factor of fsky from the sky coverage. 
	double pref = (1e0/(8e0*pi)) * pow(6e0 * deltaphi * deltaphi, 2) * fsky;

	#pragma omp parallel
	{
		double cc, cs, sc, ss;
		double xyint;
		int i, j, n, l, mu, p;

		#pragma omp for reduction(+:N_cos,N_sin,N_cross)
		for(n=start_i_j; n<end_i_j; n++){
			i = i_j[n][0];
			j = i_j[n][1];

			xyint = step_x[i] * step_x[j] * (xvec[i]*xvec[i]) * (xvec[j]*xvec[j]);
			
			for(mu=0; mu<npts_mu; mu++){
				// Sum over l and polarisation first
				cc = cs = sc = ss = 0;
				for(p=0; p<npts_p; p++){
					for(l=0; l<npts_l; l++){
						cc += cos_tilde[p][i][l] * cos_tilde[p][j][l] * legendre[mu][l];
						cs += cos_tilde[p][i][l] * sin_tilde[p][j][l] * legendre[mu][l];
						sc += sin_tilde[p][i][l] * cos_tilde[p][j][l] * legendre[mu][l];
						ss += sin_tilde[p][i][l] * sin_tilde[p][j][l] * legendre[mu][l];
					}	
				}

				N_cos += gl_weights[mu] * xyint * ((cc * cc * cc) + 3*(cc * ss * ss) - 3*(sc * sc * cc) - 3*(cs * cs * cc) + 6*(sc * cs * ss));
				N_sin += gl_weights[mu] * xyint * ((ss * ss * ss) + 3*(ss * cc * cc) - 3*(cs * cs * ss) - 3*(sc * sc * ss) + 6*(cs * sc * cc));
				N_cross += gl_weights[mu] * xyint * (3*(cs * cc * cc) - (cs * cs * cs) - 3*(cs * sc * sc) - 6*(cc * ss * sc) + 3*(cs * ss * ss));
			}
		}
	}

	N_cos *= pref;
	N_sin *= pref;
	N_cross *= pref;

	sigma_cos = sqrt(6e0/N_cos);
	sigma_sin = sqrt(6e0/N_sin);
}

double error_bar_with_phase(double phi){
	/* Returns the error bar for S = sin(omega*k + phi) */

	double c = cos(phi), s = sin(phi);
	double N = (s*s) * N_cos + (c*c) * N_sin + (2*c*s) * N_cross;
	return sqrt(6e0/N);
}


void print_results(){
	///* Prints parameters and computed error bars*/
	printf("Computing error bars for sinusodial shape functions\n");
	printf("*** Parameters ***\n");
	printf("omega = %e\n", omega);
	printf("lmin = %d, lmax = %d\n", lmin, lmax);
	printf("xmax = %e, npts_x = %d\n", xmax, npts_x);
	printf("kmax = %e, npts_k = %d\n", kmax, npts_k);
	printf("npts_mu = %d\n", npts_mu);
	printf("Amplitude factor deltaphi = %e\n", deltaphi);
	printf("Used %s, %s, %s, %s\n", bessel_data_filename, transfer_T_data_filename, C_data_filename, BN_TT_data_filename );
	printf("Polarisation %s\n", (do_polarisation ? "ON" : "OFF"));
	printf("\n* Results *\n");
	printf("N_cos = %e, N_sin = %e, N_cross = %e\n", N_cos, N_sin, N_cross);
	printf("sigma_cos = %e, sigma_sin = %e\n\n", sigma_cos, sigma_sin);
}


void free_tilde(){
	///* Free up the memory of cos, sin tildes */	
	int p;
	for(p=0; p<npts_p; p++){
		free_2D_array(cos_tilde[p]);
		free_2D_array(sin_tilde[p]);
	}
}


///***** All the ingredients for computing error bars are now ready. Different computational routines to follow. *****/


void error_bars(){
	///* Pure openMP routine to compute error bars for a single value of omega */

	clock_t start = clock();
	
	// Initialise and load data
	initialise();
	
	// Perform the k integral to compute cos_tilde(x) and sin_tilde(x)
	precompute_tilde();
	printf("Finished computing cos_tilde, sin_tilde. Elapsed time: %fs\n", (double)(clock() - start) / CLOCKS_PER_SEC);

	// We no longer need bessel and transfer data
	free_bessel();
	free_transfer();

	// Now modify sin, cos tilde basis so that the covariance matrix is orthonormal
	orthogonalise_tilde(); 

	// Compute the Gauss-Legendre coefficients and store Legendre polynomial values 
	prepare_mu_integration();	

	// Now the main computation	
	compute_error_bars();
	
	// Print out the answer together with parameters used
	print_results();

	// Clean up
	free_tilde();
	free_2D_array(C);
	free_2D_array(legendre);
	free_array(gl_nodes);
	free_array(gl_weights);
	free_2D_int_array(i_j);
	free_2D_int_array(i_l);
	
	clock_t end = clock();
	double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Elapsed time: %fs\n", time_spent);
	
}



void multiple_omega(){
	/* Pure openMP routine to compute error bars for multiple values of omega */

	// Initialise and load data
	initialise();
	
	// Compute the Gauss-Legendre coefficients and store Legendre polynomial values 
	prepare_mu_integration();	

	int o;
	double **results = create_2D_array(npts_o, 3);
	for(o=0; o<npts_o; o++){
		omega = step_o * o;
		results[o][0] = omega;
	
		// Perform the k integral to compute cos_tilde(x) and sin_tilde(x) for this value of omega
		precompute_tilde();

		// Now modify sin, cos tilde basis so that the covariance matrix is orthonormal
		orthogonalise_tilde(); 

		// Now the main computation	
		compute_error_bars();
		
		// Print out the answer together with parameters used
		print_results();

		// Record the results
		results[o][1] = sigma_cos;
		results[o][2] = sigma_sin;
	}

	write_2D_array(results, npts_o, 3, "multiple_omega_results");
	
	// Clean up
	free_tilde();
	free_2D_array(C);
	free_2D_array(legendre);
	free_2D_array(results);
	free_array(gl_nodes);
	free_array(gl_weights);
	free_2D_int_array(i_j);
	free_2D_int_array(i_l);
	free_bessel();
	free_transfer();
}


void mpi_error_bars(int argc, char **argv){
	/* MPI + openMP routine to compute error bars for a single value of omega */

	// MPI parameters
	int rank, nprocs;

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Initialise and load data
	initialise();

	printf("Rank #%d finished initialising\n", rank);
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Distribute the workload by setting start and end points for vectorised indices.
	start_i_l = rank * (npts_x_l / nprocs) + fmin(rank, npts_x_l % nprocs);
	end_i_l = (rank + 1) * (npts_x_l / nprocs) + fmin(rank + 1, npts_x_l % nprocs);
	start_i_j = rank * (npts_x_y / nprocs) + fmin(rank, npts_x_y % nprocs);
	end_i_j = (rank + 1) * (npts_x_y / nprocs) + fmin(rank + 1, npts_x_y % nprocs);

	// Perform the k integral to compute cos_tilde(x) and sin_tilde(x)
	precompute_tilde();

	// We no longer need bessel and transfer data
	free_bessel();
	free_transfer();

	// Now modify sin, cos tilde basis so that the covariance matrix is orthonormal
	orthogonalise_tilde(); 

	printf("Rank #%d finished precomputing and orthogonalising\n", rank);
	MPI_Barrier(MPI_COMM_WORLD);

	// Now incorporate the results from different ranks

	int p, i, l;
	double **temp_cos_tilde = create_2D_array(npts_x, npts_l);
	double **temp_sin_tilde = create_2D_array(npts_x, npts_l);

	for(p=0; p<npts_p; p++){

		MPI_Allreduce(&cos_tilde[p][0][0], &temp_cos_tilde[0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&sin_tilde[p][0][0], &temp_sin_tilde[0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		for(i=0; i<npts_x; i++){
			for(l=0; l<npts_l; l++){
				cos_tilde[p][i][l] = temp_cos_tilde[i][l];
				sin_tilde[p][i][l] = temp_sin_tilde[i][l];
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	free_2D_array(temp_cos_tilde);
	free_2D_array(temp_sin_tilde);

	// Compute the Gauss-Legendre coefficients and store Legendre polynomial values 
	prepare_mu_integration();	

	// Now the main computation	
	compute_error_bars();

	printf("Rank #%d finished computing error bars\n", rank);
	MPI_Barrier(MPI_COMM_WORLD);

	// Incorporate the results in the root process
	double final_N_cos, final_N_sin, final_N_cross;

	MPI_Reduce(&N_cos, &final_N_cos, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&N_sin, &final_N_sin, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&N_cross, &final_N_cross, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0){

		N_cos = final_N_cos;
		N_sin = final_N_sin;
		N_cross = final_N_cross;
		sigma_cos = sqrt(6e0/N_cos);
		sigma_sin = sqrt(6e0/N_sin);

		// Print out the answer together with parameters used
		print_results();

		// Save xvec and kvec data
		write_1D_array(kvec, npts_k, "kvec");
		write_1D_array(xvec, npts_x, "xvec");

	}

	// Clean up
	free_tilde();
	free_2D_array(C);
	free_2D_array(legendre);
	free_array(gl_nodes);
	free_array(gl_weights);
	free_2D_int_array(i_j);
	free_2D_int_array(i_l);

	MPI_Finalize();
}


void mpi_multiple_omega(int argc, char **argv){
	/* MPI + openMP routine to compute error bars for multiple values of omega */

	// MPI parameters
	int rank, nprocs;

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Initialise and load data
	initialise();

	// Compute the Gauss-Legendre coefficients and store Legendre polynomial values 
	prepare_mu_integration();	

	printf("Rank #%d finished initialising\n", rank);
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Distribute the workload by setting start and end points for vectorised indices.
	start_i_l = rank * (npts_x_l / nprocs) + fmin(rank, npts_x_l % nprocs);
	end_i_l = (rank + 1) * (npts_x_l / nprocs) + fmin(rank + 1, npts_x_l % nprocs);
	start_i_j = rank * (npts_x_y / nprocs) + fmin(rank, npts_x_y % nprocs);
	end_i_j = (rank + 1) * (npts_x_y / nprocs) + fmin(rank + 1, npts_x_y % nprocs);

	// Temporary storage spaces for MPI communications
	double **temp_cos_tilde = create_2D_array(npts_x, npts_l);
	double **temp_sin_tilde = create_2D_array(npts_x, npts_l);
	double final_N_cos, final_N_sin, final_N_cross;

	// Phases specification
	int npts_phi = 10;
	double *phase = create_array(npts_phi);
	int phi;
	for(phi=0; phi<npts_phi; phi++){
		phase[phi] = (pi / npts_phi) * phi;
	}

	double **results = create_2D_array(npts_o, npts_phi+1);
	int o, p, i, l;

	for(o=0; o<npts_o; o++){

		omega = step_o * o;
		results[o][0] = omega;

		// Perform the k integral to compute cos_tilde(x) and sin_tilde(x)
		precompute_tilde();

		// Modify sin, cos tilde basis so that the covariance matrix is orthonormal
		orthogonalise_tilde(); 

		MPI_Barrier(MPI_COMM_WORLD);

		// Now incorporate the results from different ranks
		for(p=0; p<npts_p; p++){

			MPI_Allreduce(&cos_tilde[p][0][0], &temp_cos_tilde[0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&sin_tilde[p][0][0], &temp_sin_tilde[0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			for(i=0; i<npts_x; i++){
				for(l=0; l<npts_l; l++){
					cos_tilde[p][i][l] = temp_cos_tilde[i][l];
					sin_tilde[p][i][l] = temp_sin_tilde[i][l];
				}
			}

			MPI_Barrier(MPI_COMM_WORLD);
		}

		// Now the main computation	
		compute_error_bars();

		MPI_Barrier(MPI_COMM_WORLD);

		// Incorporate the results in the root process

		MPI_Reduce(&N_cos, &final_N_cos, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&N_sin, &final_N_sin, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&N_cross, &final_N_cross, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0){

			N_cos = final_N_cos;
			N_sin = final_N_sin;
			N_cross = final_N_cross;

			sigma_cos = sqrt(6e0/N_cos);
			sigma_sin = sqrt(6e0/N_sin);

			for(phi=0; phi<npts_phi; phi++){
				results[o][phi+1] = error_bar_with_phase(phase[phi]);
			}

			// Print out the answer together with parameters used
			print_results();

		}
	}

	if(rank == 0){
		write_2D_array(results, npts_o, npts_phi+1, "mpi_multiple_omega_and_phase_result");
	}

	// Clean up
	free_tilde();
	free_2D_array(C);
	free_2D_array(legendre);
	free_2D_array(temp_cos_tilde);
	free_2D_array(temp_sin_tilde);
	free_2D_array(results);
	free_array(gl_nodes);
	free_array(gl_weights);
	free_array(phase);
	free_2D_int_array(i_j);
	free_2D_int_array(i_l);
	free_bessel();
	free_transfer();

	MPI_Finalize();
}


//
//void bispectrum_plot(){
//	// Calculates bispectrum for l1=l2=l3=l and saves the results to a file
//
//	// Initialise
//	initialise();
//
//	precompute_tilde();
//
//	free_bessel();
//	free_transfer();
//
//	double *b_cos = (double *) calloc(npts_l, sizeof(double));
//	double *b_sin = (double *) calloc(npts_l, sizeof(double));
//	
//	double c, s;
//	int l, i;
//	#pragma omp parallel for private(l,i,c,s)
//	for(l=0; l<npts_l; l++){
//		for(i=0; i<npts_x; i++){
//			c = cos_tilde[T][i][l];
//			s = sin_tilde[T][i][l];
//			b_cos[l] += xvec[i] * xvec[i] * step_x[i] * (c*c*c - 3*c*s*s);
//			b_sin[l] += xvec[i] * xvec[i] * step_x[i] * (-s*s*s + 3*s*c*c);
//		}
//	}
//
//	write_1D_array(b_cos, npts_l, "cos_bispectrum_equal_l");
//	write_1D_array(b_sin, npts_l, "sin_bispectrum_equal_l");
//}
//
//void const_model(){
//	// Computes the bispectrum for constant shape function S to compare with analytic solution
//
//	lmax = 200;
//	initialise();
//	double **transfer = create_2D_array(npts_l, npts_k);
//
//	// Large angle approximation for l<<200. Transfer(l,k) ~= (1/3)*bessel(l,tau0*k)
//	gsl_interp_accel *acc = gsl_interp_accel_alloc();
//	gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, bessel_npts_x);
//	int l, k;
//	for(l=0; l<npts_l; l++){
//		gsl_spline_init(spline, bessel_xvec, get_bessel(l+lmin), bessel_npts_x);
//		for(k=0; k<transfer_npts_k; k++){
//			transfer[l][k] = (1e0/3e0) * gsl_spline_eval(spline, transfer_kvec[k] * tau0, acc);
//		}
//		gsl_interp_accel_reset(acc);	
//	}
//
//	// S=cos(wk) for omega = 0 corresponds to the constant model
//	omega = 0;
//
//	// **** THIS NEEDS TO BE CHANGED ****
//	precompute_tilde();	// Here large angle approx. to the transfer functions needs to be used
//	write_2D_array(cos_tilde[T], npts_x, npts_l, "const_tilde");
//
//
//	// Bispectrum
//	double ***b_const = create_3D_array(npts_l, npts_l, npts_l);
//	double ***b_analytic = create_3D_array(npts_l, npts_l, npts_l);
//	double threshold = 1e-3;
//
//	double *xint = create_array(npts_x);
//	int i;
//	for(i=0; i<npts_x; i++){
//		xint[i] = step_x[i] * xvec[i] * xvec[i];	// to save computing time later
//	}
//
//	#pragma omp parallel private(i)
//	{
//		int l1, l2, l3;
//		int ii, jj, kk;
//		double error;
//
//		#pragma omp for
//		for(ii=0; ii<npts_l; ii++){
//			l1 = ii + lmin;
//
//			for(jj=0; jj<=ii; jj++){
//				l2 = jj + lmin;
//
//				// triangle condition and l1>=l2>=l3 enforced here
//				for(kk=abs(l1-l2)-lmin; kk<=jj; kk++){
//					l3 = kk + lmin;
//
//					b_const[ii][jj][kk] = 0;
//					for(i=0; i<npts_x; i++){
//						b_const[ii][jj][kk] += xint[i] * cos_tilde[T][i][ii] * cos_tilde[T][i][jj] * cos_tilde[T][i][kk];
//					}
//
//					b_analytic[ii][jj][kk] = 1e0 / (27e0 * (2*l1+1) * (2*l2+1) * (2*l3+1)) * (1e0 / (l1+l2+l3+3e0) + 1e0 / (l1+l2+l3));
//
//					error = (b_const[ii][jj][kk] - b_analytic[ii][jj][kk]) / b_analytic[ii][jj][kk];
//
//					if(fabs(error) > threshold){
//						printf("b(%d,%d,%d) = %e, analytic %e, error %e\n", l1, l2, l3, b_const[ii][jj][kk], b_analytic[ii][jj][kk], error);
//					}
//				}
//			}
//		}
//		
//	}
//
//}
//
//
//
//void write_qtilde_integrand(){
//
//	initialise();
//
//	///* Save integrands for k the integration */
//
//	int npts_l_save = 10;
//
//	double ***cos_integrand = create_3D_array(npts_l_save, npts_x, npts_k);
//	double ***sin_integrand = create_3D_array(npts_l_save, npts_x, npts_k);
//
//	free_2D_array(cos_tilde[T]);
//	free_2D_array(sin_tilde[T]);
//
//	cos_tilde[T] = create_2D_array(npts_x, npts_l_save);
//	sin_tilde[T] = create_2D_array(npts_x, npts_l_save);
//
//	// Initialise to zero
//	int i, l;
//	for(i=0; i<npts_x; i++){
//		for(l=0; l<npts_l_save; l++){
//			cos_tilde[T][i][l] = sin_tilde[T][i][l] = 0;
//		}
//	}
//
//	// Perform integral over k to evaluate cos_tilde(x,l) and  sin_tilde(x,l)
//	// cos_tilde(x,l) = (2/pi) * integral{dk * cos(omega*k) * j_l(k*x) * Delta_l(k)}
//	// sin_tilde(x,l) = (2/pi) * integral{dk * sin(omega*k) * j_l(k*x) * Delta_l(k)}
//
//	#pragma omp parallel private(i,l)
//	{
//		// Initialise GSL tools for interpolation and integration
//		gsl_spline *bessel_spline = gsl_spline_alloc(gsl_interp_linear, bessel_npts_x);
//		gsl_interp_accel *bessel_acc = gsl_interp_accel_alloc();
//		gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, npts_k);
//		gsl_interp_accel *acc = gsl_interp_accel_alloc();
//		double bes;	// temporary variables
//		double *transfer;
//		double *y1 = create_array(npts_k);
//		double *y2 = create_array(npts_k);
//		int k;
//
//		#pragma omp for
//		for(l=0; l<npts_l_save; l++){
//			for(i=0; i<npts_x; i++){
//
//				gsl_spline_init(bessel_spline, bessel_xvec, get_bessel(l+lmin), bessel_npts_x);
//				transfer = get_transfer(l+lmin);
//
//				for(k=0; k<npts_k; k++){
//					bes = gsl_spline_eval(bessel_spline, xvec[i] * kvec[k], bessel_acc);
//
//					y1[k] = cos(omega * kvec[k]) * transfer[k] * bes;
//					y2[k] = sin(omega * kvec[k]) * transfer[k] * bes;
//
//					cos_integrand[l][i][k] = y1[k];
//					sin_integrand[l][i][k] = y2[k];
//				}
//
//				gsl_interp_accel_reset(bessel_acc);
//
//				gsl_spline_init(spline, kvec, y1, npts_k);
//				cos_tilde[T][i][l] = (2e0/pi) * gsl_spline_eval_integ(spline, kmin, kmax, acc);
//				gsl_interp_accel_reset(acc);
//
//				gsl_spline_init(spline, kvec, y2, npts_k);
//				sin_tilde[T][i][l] = (2e0/pi) * gsl_spline_eval_integ(spline, kmin, kmax, acc);
//				gsl_interp_accel_reset(acc);
//			}
//		}
//
//		// Free the interpolaters
//		gsl_spline_free(bessel_spline);
//		gsl_interp_accel_free(bessel_acc);
//		gsl_spline_free(spline);
//		gsl_interp_accel_free(acc);
//	}
//
//	// Write data
//	for(l=0; l<npts_l_save; l++){
//		char name1[200], name2[200];
//		sprintf(name1, "cos_tilde_integrand_l_%d", l+lmin);
//		sprintf(name2, "sin_tilde_integrand_l_%d", l+lmin);
//		write_2D_array(cos_integrand[l], npts_x, npts_k, name1);
//		write_2D_array(sin_integrand[l], npts_x, npts_k, name2);
//	}	
//}
//
//
//
//
//
//void mpi_bispectrum(int argc, char **argv){
//
//	// MPI parameters
//	int rank, nprocs;
//
//	MPI_Init(&argc, &argv); 
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
//
//	// Initialise parameters and load data
//	initialise();
//
//
//	//	printf("Proc %d finished initialising\n", rank);
//
//	// Distribute the workload by setting start and end points for vectorised indices.
//	// Note that for non-mpi codes they are set to 0 and npts, respectively
//	start_i_l = rank * (npts_x_l / nprocs) + fmin(rank, npts_x_l % nprocs);
//	end_i_l = (rank + 1) * (npts_x_l / nprocs) + fmin(rank + 1, npts_x_l % nprocs);
//
//	//	printf("Proc %d finished loading bessel & transfer\n", rank);
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	// Precompute sin_tilde(x) and cos_tilde(x)
//	precompute_tilde();
//	printf("Proc %d Finished computing sin_tilde, cos_tilde.\n", rank);
//
//	//	printf("Proc %d finished precomputing\n", rank);
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	// In order to put the data together, we define temporary arrays
//	double **temp_cos_tilde[T] = create_2D_array(npts_x, npts_l);
//	double **temp_sin_tilde[T] = create_2D_array(npts_x, npts_l);
//
//	// MPI_Allreduce to incorporate calculated results
//	MPI_Allreduce(&cos_tilde[T][0][0], &temp_cos_tilde[T][0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//	MPI_Allreduce(&sin_tilde[T][0][0], &temp_sin_tilde[T][0][0], npts_x * npts_l, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//
//	int i, l;	
//	for(i=0; i<npts_x; i++){
//		for(l=0; l<npts_l; l++){
//			cos_tilde[T][i][l] = temp_cos_tilde[T][i][l];
//			sin_tilde[T][i][l] = temp_sin_tilde[T][i][l];
//		}
//	}
//	free_2D_array(temp_cos_tilde[T]);
//	free_2D_array(temp_sin_tilde[T]);
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	// We no longer need bessel and transfer functions data
//	free_bessel();
//	free_transfer();
//
//	// Bispectrum
//	int lstep = 10;
//	npts_l = npts_l / lstep;
//	lmax = lmin + (npts_l - 1) * lstep;
//	double ***bispectrum = create_3D_array(npts_l, npts_l, npts_l);
//
//	// Initialise to zero
//	int l1, l2, l3;
//	for(l1=0; l1<npts_l; l1++){
//		for(l2=0; l2<npts_l; l2++){
//			for(l3=0; l3<npts_l; l3++){
//				bispectrum[l1][l2][l3]= 0;
//			}
//		}
//	}
//
//	// Distribute workload
//	int **l_l = vectorise_two_indices(npts_l, npts_l);	
//	int npts_l_l = npts_l * npts_l;
//	int start_l_l, end_l_l;
//	start_l_l = rank * (npts_l_l / nprocs) + fmin(rank, npts_l_l % nprocs);
//	end_l_l = (rank + 1) * (npts_l_l / nprocs) + fmin(rank + 1, npts_l_l % nprocs);
//
//
//#pragma omp parallel private(i,l1,l2,l3)
//	{
//		gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, npts_x);
//		gsl_interp_accel *acc = gsl_interp_accel_alloc();
//		double *y = create_array(npts_x);
//		int n;
//
//#pragma omp for
//		for(n=start_l_l; n<end_l_l; n++){
//			l1 = l_l[n][0];
//			l2 = l_l[n][1];
//
//			for (l3=0; l3<npts_l; l3++){
//				if (lstep*(l1+l2)+lmin >= lstep*l3 && lstep*(l2+l3)+lmin >= lstep*l1 && lstep*(l3+l1)+lmin >= lstep*l2){ 
//					for(i=0; i<npts_x; i++){
//						y[i] = xvec[i] * xvec[i] * (cos_tilde[T][i][lstep*l1]*cos_tilde[T][i][lstep*l2]*cos_tilde[T][i][lstep*l3] - cos_tilde[T][i][lstep*l1]*sin_tilde[T][i][lstep*l2]*sin_tilde[T][i][lstep*l3] - sin_tilde[T][i][lstep*l1]*cos_tilde[T][i][lstep*l2]*sin_tilde[T][i][lstep*l3] - sin_tilde[T][i][lstep*l1]*sin_tilde[T][i][lstep*l2]*cos_tilde[T][i][lstep*l3]);
//
//					}
//
//					gsl_spline_init(spline, xvec, y, npts_x);
//					bispectrum[l1][l2][l3] = gsl_spline_eval_integ(spline, xvec[0], xvec[npts_x-1], acc);	 
//					gsl_interp_accel_reset(acc);
//				}
//			}
//		}
//
//		gsl_spline_free(spline);
//		gsl_interp_accel_free(acc);
//		free_array(y);
//
//	}
//
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	// Again need temporary variables
//	double ***final_bispectrum = create_3D_array(npts_l, npts_l, npts_l);
//
//	// MPI reduce the bisepctrum in root. 
//	MPI_Reduce(&bispectrum[0][0][0], &final_bispectrum[0][0][0], npts_l*npts_l*npts_l, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//
//	if (rank == 0){
//		char *filename = "cos_bispectrum_3D_2000";
//		write_bispectrum(final_bispectrum, npts_l, lstep, lstep, filename);
//
//		printf("Computing bispectrum for feature models\n");
//		printf("*** Parameters ***\n");
//		printf("omega = %e\n", omega);
//		printf("lmin = %d, lmax = %d\n", lmin, lmax);
//		printf("xmax = %e, npts_x = %d\n", xmax, npts_x);
//		printf("kmax = %e, npts_k = %d\n", kmax, npts_k);
//		printf("Used bessel_4000 data set and transfer functions,  power spectrum from Planck\n");
//		printf("\n");
//		printf("Results saved to file %s\n", filename);
//	}
//
//	free_2D_array(cos_tilde[T]);
//	free_2D_array(sin_tilde[T]);
//
//
//	MPI_Finalize();
//
//
//}
//
//*/

int main(int argc, char **argv){
	mpi_error_bars(argc, argv);
//	mpi_multiple_omega(argc, argv);
//	error_bars();
//	multiple_omega();
//	bispectrum_plot();
//	mpi_bispectrum(argc, argv);
//	const_model();	
//	write_qtilde_integrand();

	return 0;
}

