#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "arrays.h"

double pi = M_PI;

double deltaphi = 1.5714e-8;
double fskyT = 0.76;
double fskyE = 0.74;

//char *data_dir = "/fast/space/projects/dp002/jf334/data/";
char *data_dir = "./data/";

char *bessel_data_filename = "bessel_4000";
char *bessel_size_filename = "bessel_4000_size";
char *transfer_T_data_filename = "transfer_planck_DDX9_3000_T";
char *transfer_E_data_filename = "transfer_planck_DDX9_3000_E";
char *transfer_size_filename = "transfer_planck_DDX9_3000_size";

char *plot_dir = "./plot_data/";

int bessel_npts_l, bessel_npts_x;
int bessel_lmin, bessel_lmax;
double *bessel_xvec;
double **bessel;
double tau0; 	// distance to the surface of last scattering 

int transfer_npts_l, transfer_npts_k;
int transfer_lmin, transfer_lmax;
double *transfer_kvec;
double **transfer;

int cl_lmin, cl_lmax, cl_npts_l;
double *C;

int BN_lmin, BN_lmax, BN_npts_l;
double *beam;
double *noise;

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

	// Now shift the array so that the first row contains actual data
	for(l=0; l<bessel_npts_l; l++){
		bessel[l] = bessel[l+1];
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

	transfer = create_2D_array(sizes[0], sizes[1]);
	fread(transfer[0], sizeof(double), size, T_data_file);
	fclose(T_data_file);

	// The first row of transfer functions data contains values of k (exponential),
	// while the first column of the data contains l values
	// It is implied that the transfer function vanishes at k = 0

	transfer_npts_l = sizes[0] - 1;
	transfer_npts_k = sizes[1];

	transfer_kvec = &transfer[0][0];	// A vector containing k values of the transfer functions data
	transfer_lmin = (int) transfer[1][0];
	transfer_lmax = (int) transfer[transfer_npts_l][0];
	int l;
	for(l=1; l<=transfer_npts_l; l++){
		if(transfer[l][0] != transfer_lmin + l - 1){
			printf("Fault in the transfer function data file. Row #%d has l=%e, not %d\n", l, transfer[l][0], transfer_lmin+l-1);
			exit(1);
		}else{
			transfer[l][0] = 0;
		}
	}

	// Now shift the array so that the first row contains actual data
	for(l=0; l<transfer_npts_l; l++){
		transfer[l] = transfer[l+1];
	}

//	printf("Transfer. kmin = %e, kmax = %e, lmin = %d, lmax = %d\n", transfer_kvec[0], transfer_kvec[transfer_npts_k-1], transfer_lmin, transfer_lmax);

}

void load_cls(){

	FILE *angular_power_spectrum_file;
	char angular_power_spectrum_filename[100];
	char line[200], **cptr = (char **) malloc(1 * sizeof(char *));
	C = create_array(3499);

	strcat(strcpy(angular_power_spectrum_filename, data_dir), "cls_planck_DDX9_lensed.txt");
	angular_power_spectrum_file = fopen(angular_power_spectrum_filename, "r");
	if (angular_power_spectrum_file == NULL){
		printf("Error loading angular power spectrum file. Try again, you can do it! \n");
		exit(1);
	}

	int i, l;
	i = 0;
	while (fgets(line, 200, angular_power_spectrum_file)){

		l = (int) strtod(line, cptr);
		if (i == 0) cl_lmin = l;
		C[i] = strtod(*cptr, cptr);
		C[i] = 2e0 * pi * C[i] / (l * (l + 1e0));
		i++;
	}

	cl_lmax = l;
	cl_npts_l = i;
	fclose(angular_power_spectrum_file);

//	printf("Angular power spectrum data file has %d points in l\n", cl_npts_l);
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
		printf("Error loading beam and noise file. Try again, you can do it!\n");
		exit(1);
	}

	int i, l;
	i = 0;
	while (fgets(line, 200, BN_file)){
		
		l = (int) strtod(line, cptr);
		if (l == 0) BN_lmin = l;
		beam[i] = strtod(*cptr, cptr);
		noise[i] = strtod(*cptr, cptr);
		i++;
	}

	BN_lmax = l;
	BN_npts_l = i;
	fclose(BN_file);

}
