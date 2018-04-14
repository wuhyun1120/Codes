#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *output_dir = "/home/wuhyun/Treasure/plot_data/";
//char *output_dir = "./plot_data/";

void write_1D_array(double *data, int size, char *output_filename){

	char filename[100];
	strcat(strcpy(filename, output_dir), output_filename);
	FILE *output_file = fopen(filename, "w");
	if (output_file == NULL){
		printf("Error opening file %s to write data. Try again, you can do it!\n", output_filename);
		exit(1);
	}
	
	int i;
	for(i=0; i<size; i++){
		fprintf(output_file, "%e\n", data[i]);
	}
	fclose(output_file);

}

void write_2D_array(double **data, int size1, int size2, char *output_filename){

	char filename[100];
	strcat(strcpy(filename, output_dir), output_filename);
	FILE *output_file = fopen(filename, "w");
	if (output_file == NULL){
		printf("Error opening file %s to write data. Try again, you can do it!\n", output_filename);
		exit(1);
	}
	
	int i,j;
	for(i=0; i<size1; i++){
		for(j=0; j<size2; j++){
			fprintf(output_file, "%e ", data[i][j]);
		}
		fprintf(output_file, "\n");
	}
	fclose(output_file);

}

void write_bispectrum(double ***data, int size, int lmin, int lstep, char *output_filename){

	char filename[100];
	strcat(strcpy(filename, output_dir), output_filename);
	FILE *output_file = fopen(filename, "w");
	if (output_file == NULL){
		printf("Error opening file %s to write data. Try again, you can do it!\n", output_filename);
		exit(1);
	}
	
	int l1, l2, l3;
	for(l1=0; l1<size; l1++){
		for(l2=0; l2<=l1; l2++){
			for(l3=0; l3<=l2; l3++){
				if (lstep*(l1+l2)+lmin>=lstep*l3 && lstep*(l2+l3)+lmin>=lstep*l1 && lstep*(l3+l1)+lmin>=lstep*l2){
					fprintf(output_file, "%d %d %d %e\n", lstep*l1+lmin, lstep*l2+lmin, lstep*l3+lmin, data[l1][l2][l3]);
				}

			}
		}
	}
	fclose(output_file);
}
