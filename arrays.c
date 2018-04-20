#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

const int use_icc = 0;

double *create_array(int d1){
	double *arr = (double *) malloc(d1 * sizeof(double));
	if (arr == NULL){
		printf("Insuffcient memory for creating a %d array. Try again, you can do it!\n", d1);
		exit(1);
	}
	return arr;
}

double **create_2D_array(int d1, int d2){
	double **arr = (double **) malloc(d1 * sizeof(double *));
	arr[0] = (double *) malloc(d1 * d2 * sizeof(double));
	if (arr[0] == NULL){
		printf("Insuffcient memory for creating a %d X %d array. Try again, you can do it!\n", d1, d2);
		exit(1);
	}
	int i;
	for(i=1; i<d1; i++) arr[i] = arr[i-1] + d2;
	return arr;
}

double ***create_3D_array(int d1, int d2, int d3){
	double ***arr = (double ***) malloc(d1 * sizeof(double **));
	arr[0] = (double **) malloc(d1 * d2 * sizeof(double *));
	arr[0][0] = (double *) malloc(d1 * d2 * d3 * sizeof(double));
	if (arr[0][0] == NULL){
		printf("Insufficient memory for creating a %d X %d X %d array. Try again, you can do it!\n", d1, d2, d3);
		exit(1);
	}
	int i, j, n=0;
	for(i=1; i<d1; i++) arr[i] = arr[i-1] + d2;
	for(i=0; i<d1; i++){
		for(j=0; j<d2; j++){
			arr[i][j] = arr[0][0] + d3 * n;
			n++;
		}
	}
	return arr;
}


double *create_aligned_array(int d1){
	assert(d1 % 8 == 0);
	double *arr = (double *) aligned_alloc(64, d1 * sizeof(double));
	if (arr == NULL){
		printf("Insuffcient memory for creating a %d array. Try again, you can do it!\n", d1);
		exit(1);
	}
	return arr;
}


double **create_aligned_2D_array(int d1, int d2){
	assert(d1 % 8 == 0);
	assert(d2 % 8 == 0);
	double **arr = (double **) aligned_alloc(64, d1 * sizeof(double *));
	arr[0] = (double *) aligned_alloc(64, d1 * d2 * sizeof(double));
	if (arr[0] == NULL){
		printf("Insuffcient memory for creating a %d X %d array. Try again, you can do it!\n", d1, d2);
		exit(1);

	}
}


void free_array(double *arr){
	free(arr);
}

void free_2D_array(double **arr){
	free(arr[0]);
	free(arr);
}

void free_aligned_array(double *arr){
	free(arr);
}

void free_aligned_2D_array(double **arr){
	free(arr[0]);
	free(arr);
}

void free_int_array(int * arr){
	free(arr);
}

void free_2D_int_array(int **arr){
	free(arr[0]);
	free(arr);
}

int **vectorise_two_indices(int size1, int size2){
	int **vec_ind = (int **) malloc(size1 * size2 * sizeof(int *));
	vec_ind[0] = (int *) malloc(2 * size1 * size2 * sizeof(int));
	int i, j, n = 0;
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
	int i, j, k, n = 0;
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

