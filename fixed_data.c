#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "arrays.h"

#ifndef T
#define T 0
#define E 1
#define TT 0
#define TE 1
#define EE 2
#endif

double pi = M_PI;

double deltaphi = 1.5714e-8;
double fskyT = 0.76;
double fskyE = 0.74;

int do_polarisation = 0;

//char *data_dir = "/fast/space/projects/dp002/jf334/data/";
char *data_dir = "./data/";

char *bessel_data_filename = "bessel_4000";
char *bessel_size_filename = "bessel_4000_size";
char *transfer_T_data_filename = "transfer_planck_DDX9_3000_T";
char *transfer_E_data_filename = "transfer_planck_DDX9_3000_E";
char *transfer_size_filename = "transfer_planck_DDX9_3000_size";
char *C_data_filename = "cls_planck_DDX9_lensed.txt";
char *BN_TT_data_filename = "BN_DX11d_smica_case1_HP_T.txt";
char *BN_EE_data_filename = "BN_DX11d_smica_case1_HP_E.txt";

int bessel_npts_l, bessel_npts_x;
int bessel_lmin, bessel_lmax;
double *bessel_xvec;
static double **bessel;
double tau0; 	// distance to the surface of last scattering 

int transfer_npts_l, transfer_npts_k;
int transfer_lmin, transfer_lmax;
double *transfer_kvec;
static double **transfer_T, **transfer_E;

int cl_lmin, cl_lmax, cl_npts_l;
static double *C_TT, *C_TE, *C_EE;

int BN_lmin, BN_lmax, BN_npts_l;
static double *beam_TT, *beam_TE, *beam_EE;
static double *noise_TT, *noise_TE, *noise_EE;

void load_bessel(){

	FILE *size_file, *data_file;
	char size_filename[100], data_filename[100];
	int sizes[2], size;

	strcat(strcpy(size_filename, data_dir), bessel_size_filename);
	strcat(strcpy(data_filename, data_dir), bessel_data_filename);
	size_file = fopen(size_filename, "r");
	data_file = fopen(data_filename, "r");
	if (size_file == NULL || data_file == NULL){
		printf("Error loading bessel size and data files. Try again, you can do it!\n");
		exit(1);
	}

	fread(sizes, sizeof(int), 2, size_file);
	fclose(size_file);

	size = sizes[0] * sizes[1];	// total size
	// printf("Bessel data file has sizes %d X %d = %d\n", sizes[0], sizes[1], size);

	bessel = create_2D_array(sizes[0], sizes[1]);
	fread(bessel[0], sizeof(double), size, data_file);
	fclose(data_file);

	// The first row of bessel data contains the tau0 and x values (linear)
	// while the first column contains the l values
	// It is implied that the bessel function vanishes at x = 0

	bessel_npts_l = sizes[0] - 1;
	bessel_npts_x = sizes[1];

	tau0 = bessel[0][0];
	bessel[0][0] = 0;
	bessel_xvec = &bessel[0][0];	// A vector containing x values of the bessel data
	bessel_lmin = (int) bessel[1][0];
	bessel_lmax = (int) bessel[bessel_npts_l][0];
	int l;
	for(l=1; l<=bessel_npts_l; l++){
		if(bessel[l][0] != l + bessel_lmin - 1){
			printf("Fault in the bessel data file. %dth row contains l=%f, not %d\n", l+1, bessel[l][0], l+bessel_lmin-1);
			exit(1);
		}else{
			bessel[l][0] = 0;
		}
	}

//	 printf("tau0 = %e, xmin = %e, xmax = %e, lmin = %d, lmax = %d\n", tau0, bessel_xvec[0], bessel_xvec[bessel_npts_x-1], bessel_lmin, bessel_lmax);
}

void load_transfer(){

	FILE *size_file, *T_data_file;
	char size_filename[100], T_data_filename[100];
	int sizes[2], size;

	strcat(strcpy(size_filename, data_dir), transfer_size_filename);
	strcat(strcpy(T_data_filename, data_dir), transfer_T_data_filename);
	size_file = fopen(size_filename, "r");
	T_data_file = fopen(T_data_filename, "r");
	if (size_file == NULL || T_data_file == NULL){
		printf("Error loading transfer size and data files. Try again, you can do it!\n");
		exit(1);
	}

	fread(sizes, sizeof(int), 2, size_file);
	fclose(size_file);

	size = sizes[0] * sizes[1]; // total size
	// printf("Transfer function data file has size %d X %d = %d\n", sizes[0], sizes[1], size);

	transfer_T = create_2D_array(sizes[0], sizes[1]);
	fread(transfer_T[0], sizeof(double), size, T_data_file);
	fclose(T_data_file);

	// The first row of transfer functions data contains values of k (exponential),
	// while the first column of the data contains l values
	// It is implied that the transfer function vanishes at k = 0

	transfer_npts_l = sizes[0] - 1;
	transfer_npts_k = sizes[1];

	transfer_kvec = &transfer_T[0][0];	// A vector containing k values of the transfer functions data
	transfer_lmin = (int) transfer_T[1][0];
	transfer_lmax = (int) transfer_T[transfer_npts_l][0];
	int l;
	for(l=1; l<=transfer_npts_l; l++){
		if(transfer_T[l][0] != transfer_lmin + l - 1){
			printf("Fault in the transfer function T data file. Row #%d has l=%e, not %d\n", l, transfer_T[l][0], transfer_lmin+l-1);
			exit(1);
		}else{
			transfer_T[l][0] = 0;
		}
	}

//	printf("Transfer. kmin = %e, kmax = %e, lmin = %d, lmax = %d\n", transfer_kvec[0], transfer_kvec[transfer_npts_k-1], transfer_lmin, transfer_lmax);

	if(do_polarisation == 1){

		FILE *E_data_file;
		char E_data_filename[100];
		strcat(strcpy(E_data_filename, data_dir), transfer_E_data_filename);
		E_data_file = fopen(E_data_filename, "r");
		if(E_data_filename == NULL){
			printf("Error loading transfer E data file. Try again, you can do it!\n");
			exit(1);
		}

		transfer_E = create_2D_array(sizes[0], sizes[1]);
		fread(transfer_E[0], sizeof(double), size, E_data_file);
		fclose(E_data_file);

		//Check if it is consistent with the T data file
		int k;
		for(k=0; k<transfer_npts_k; k++){
			if(fabs(transfer_E[0][k]-transfer_T[0][k])>1e-6){
				printf("Warning: transfer T and E data files have different kvec: T %e E %e\n", transfer_T[0][k], transfer_E[0][k]);
				break;
			}
		}
		for(l=1; l<=transfer_npts_l; l++){
			if(transfer_E[l][0] != transfer_lmin + l - 1){
				printf("Fault in the transfer function E data file. Row #%d has l=%e, not %d\n", l, transfer_E[l][0], transfer_lmin+l-1);
				exit(1);
			}else{
				transfer_E[l][0] = 0;
			}
		}
	}
}

void load_C(){

	FILE *angular_power_spectrum_file;
	char angular_power_spectrum_filename[100];
	char line[200], **cptr = (char **) malloc(1 * sizeof(char *));
	int size = 3499;
	C_TT = create_array(size);
	if(do_polarisation == 1){
		C_TE = create_array(size);
		C_EE = create_array(size);
	}

	strcat(strcpy(angular_power_spectrum_filename, data_dir), C_data_filename);
	angular_power_spectrum_file = fopen(angular_power_spectrum_filename, "r");
	if (angular_power_spectrum_file == NULL){
		printf("Error loading angular power spectrum file. Try again, you can do it! \n");
		exit(1);
	}

	int i, l;
	i = 0;
	while (fgets(line, 200, angular_power_spectrum_file)){

		l = (int) strtod(line, cptr);
		if (i == 0){
			cl_lmin = l;
			if(cl_lmin != 0){
				printf("Cls data file starts with l = %d instead of 0, which is assumed. Try again, you can do it! \n", cl_lmin);
				exit(1);
			}
		}
		C_TT[i] = strtod(*cptr, cptr);
		C_TT[i] = 2e0 * pi * C_TT[i] / (l * (l + 1e0));

		if(do_polarisation == 1){
			C_TE[i] = strtod(*cptr, cptr);
			C_EE[i] = strtod(*cptr, cptr);
			C_TE[i] = 2e0 * pi * C_TE[i] / (l * (l + 1e0));
			C_EE[i] = 2e0 * pi * C_EE[i] / (l * (l + 1e0));
		}
	
		i++;
	}

	

	cl_lmax = l;
	cl_npts_l = i;
	fclose(angular_power_spectrum_file);
	free(cptr);

//	printf("Angular power spectrum data file has %d points in l\n", cl_npts_l);
	//printf("lmin = %d, lmax = %d\n", cl_lmin, cl_lmax);
	//for (i=0; i<100; i++) printf("%e ", C[i]);
	//printf("\n");

}

void load_BN(){

	FILE *BN_TT_file;	// Beam and Noises
	char BN_TT_filename[100];
	char line[200], **cptr = (char **) malloc(1 * sizeof(char *));
	int size = 2501;

	beam_TT = create_array(size);
	noise_TT = create_array(size);

	strcat(strcpy(BN_TT_filename, data_dir), BN_TT_data_filename);
	BN_TT_file = fopen(BN_TT_filename, "r");

	if (BN_TT_file == NULL){
		printf("Error loading beam and noise file. Try again, you can do it!\n");
		exit(1);
	}

	int i, l;
	i = 0;
	while (fgets(line, 200, BN_TT_file)){
		
		l = (int) strtod(line, cptr);
		if (l == 0) BN_lmin = l;
		beam_TT[i] = strtod(*cptr, cptr);
		noise_TT[i] = strtod(*cptr, cptr);
		i++;
	}

	BN_lmax = l;
	BN_npts_l = i;
	fclose(BN_TT_file);

	if(do_polarisation == 1){
		FILE *BN_EE_file;
		char BN_EE_filename[100];
		strcat(strcpy(BN_EE_filename, data_dir), BN_EE_data_filename);
		BN_EE_file = fopen(BN_EE_filename, "r");

		beam_TE = create_array(size);
		beam_EE = create_array(size);
		noise_TE = create_array(size);
		noise_EE = create_array(size);
			
		if (BN_EE_file == NULL){
			printf("Error loading beam and noise file. Try again, you can do it!\n");
			exit(1);
		}

		i = 0;
		while (fgets(line, 200, BN_EE_file)){
			beam_EE[i] = strtod(*cptr, cptr);
			noise_EE[i] = strtod(*cptr, cptr);
			i++;
		}
		fclose(BN_EE_file);

		for(l=0; l<BN_npts_l; l++){
			beam_TE[l] = sqrt(beam_TT[l] * beam_EE[l]);
			noise_TE[l] = 0;
		}
	}
	
	free(cptr);

}

double *get_bessel(int l){
	// bessel[1][0] has j(l=lmin, x=0)
	return bessel[l - bessel_lmin + 1];
}

double *get_transfer(int pol, int l){
	/* pol = T or E */
	// transfer[1][0] has transfer(l=lmin, k=0)
	if(pol == T){
		return transfer_T[l - transfer_lmin + 1];
	}else{
		return transfer_E[l - transfer_lmin + 1];
	}
}

double *get_C(int pol){
	/* pol = TT, TE or EE */
	double *cp = create_array(cl_npts_l);
	int l;
	if(pol == TT){
		for(l=0; l<cl_npts_l; l++) cp[l] = C_TT[l];	
	}else if(pol == TE){
		for(l=0; l<cl_npts_l; l++) cp[l] = C_TE[l];	
	}else if(pol == EE){
		for(l=0; l<cl_npts_l; l++) cp[l] = C_EE[l];
	}
	return cp;
}

double *get_beam(int pol){
	/* pol = TT, TE or EE */
	if(pol == TT){
		return beam_TT;
	}else if(pol == TE){
		return beam_TE;
	}else{
		return beam_EE;
	}
}

double *get_noise(int pol){
	/* pol = TT, TE or EE */
	if(pol == TT){
		return noise_TT;
	}else if(pol == TE){
		return noise_TE;
	}else{
		return noise_EE;
	}
}
