#ifndef ARRAYS_H
#define ARRAYS_H

double *create_array(int d1);
double **create_2D_array(int d1, int d2);
double ***create_3D_array(int d1, int d2, int d3);

double *create_aligned_array(int d1);
double **create_aligned_2D_array(int d1, int d2);

void free_array(double *arr);
void free_2D_array(double **arr);
void free_int_array(int *arr);
void free_2D_int_array(int **arr);

int **vectorise_two_indices(int size1, int size2);
int **vectorise_three_indices(int size1, int size2, int size3);


#endif /* ARRAYS_H */
