CC = mpicc
LIBS = -lgsl -lgslcblas -lm -fopenmp
CFLAGS = -I. -g -O2
FILES = *.c

mpi_error_bar : $(FILES)
	$(CC) $(FILES) -o $@ $(CFLAGS) $(LIBS)
