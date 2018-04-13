#ifndef FIXED_DATA_H
#define FIXED_DATA_H

#define MAXLEN 200 // max length of extern char arrays

#ifndef T
#define T 0
#define E 1
#define TT 0
#define TE 1
#define EE 2
#endif

extern const double pi;

extern const double deltaphi;
extern const double fskyT;
extern const double fskyE;

extern const int do_polarisation;

extern char data_dir[MAXLEN];

extern char *bessel_data_filename;
extern char *bessel_size_filename;
extern char *transfer_T_data_filename;
extern char *transfer_E_data_filename;
extern char *transfer_size_filename;


extern int bessel_npts_l, bessel_npts_x;
extern int bessel_lmin, bessel_lmax;
extern double *bessel_xvec;
extern double tau0;

extern int transfer_npts_l, transfer_npts_k;
extern int transfer_lmin, transfer_lmax;
double *transfer_kvec;

extern int cl_lmin, cl_lmax, cl_npts_l;
extern int BN_lmin, BN_lmax, BN_npts_l;

extern void load_bessel();
extern void load_transfer();
extern void load_C();
extern void load_BN();

extern double *get_bessel(int l);
extern double *get_transfer(int pol, int l);
extern double *get_C(int pol);
extern double *get_beam(int pol);
extern double *get_noise(int pol);

extern void free_bessel();
extern void free_transfer();
extern void free_C();
extern void free_BN();


#endif /* FIXED_DATA_H */
