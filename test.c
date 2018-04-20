#include <stdio.h>
#include <stdlib.h>

static double *xvec;
static double *step_x;
static double xmax = 1.68e4;		
static int npts_x;

void make_xvec(){
	///* Makes a vector of x which is denser around the last scattering surface ~13900 */  
	
//	double step3 = 200.0;
//	double step2 = 50.0;
//	double step1 = 10.0;
	double step3 = 60.0;
	double step2 = 30.0;
	double step1 = 10.0;
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

	step_x = (double *) malloc(npts_x * sizeof(double));
//	for(i=1; i<npts_x-1; i++) step_x[i] = 0.5 * (xvec[i+1] - xvec[i-1]);
//	step_x[0] = 0.5 * (xvec[1] - xvec[0]);
//	step_x[npts_x-1] = 0.5 * (xvec[npts_x-1] - xvec[npts_x-2]);
	for(i=0; i<npts_x-1; i++) step_x[i] = xvec[i+1] - xvec[i];
	step_x[npts_x-1] = step3;

}

int main(){

	make_xvec();
	printf("npts_x = %d\n", npts_x);


	free(xvec);
	free(step_x);
	return 0;
}
