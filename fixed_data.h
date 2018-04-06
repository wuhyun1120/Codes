#ifndef FIXED_DATA_H
#define FIXED_DATA_H

#define MAXLEN 200 // max length of extern char arrays

extern double pi;

extern double deltaphi;
extern double fskyT;
extern double fskyE;

extern char data_dir[MAXLEN];

extern char bessel_data_filename[MAXLEN];
extern char bessel_size_filename[MAXLEN];
extern char transfer_T_data_filename[MAXLEN];
extern char transfer_E_data_filename[MAXLEN];
extern char transfer_size_filename[MAXLEN];


int bessel_npts_l, bessel_npts_x;
int bessel_lmin, bessel_lmax;
double *bessel_xvec;
double **bessel;
double tau0;

int transfer_npts_l, transfer_npts_k;
int transfer_lmin, transfer_lmax;
double *transfer_kvec;
double **transfer;

int cl_lmin, cl_lmax, cl_npts_l;
double *C;

int BN_lmin, BN_lmax, BN_npts_l;
double *beam;
double *noise;

void load_bessel();
void load_transfer();
void load_cls();
void load_BN();

#endif /* FIXED_DATA_H */
