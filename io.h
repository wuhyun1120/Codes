#ifndef IO_H
#define IO_H

#define MAXLEN 200 // max length of char arrays

extern char output_dir[MAXLEN];
extern void write_1D_array(double *data, int size, char *output_filename);
extern void write_2D_array(double **data, int size1, int size2, char *output_filename);
extern void write_bispectrum(double ***data, int size, int lmin, int lstep, char *output_filename);


#endif /* IO_H */
