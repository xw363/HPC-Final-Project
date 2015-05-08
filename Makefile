CC = gcc
CFLAGS = -Wall -g -fopenmp
LDFLAGS = -lfftw3 -lfftw3_omp -lcudart -lcufft -lm
OBJ = welch.o fftw.o cufft.o utility.o

.PHONY: clean

all: welch-fftw welch-fftw-openmp welch-cufft welch-cufft-openmp
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
welch-fftw: welch-fftw.o $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
welch-fftw-openmp: welch-fftw-openmp.o $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
welch-cufft: welch-cufft.o $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
welch-cufft-openmp: welch-cufft-openmp.o $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm *.o
	rm welch-fftw
	rm welch-fftw-openmp
	rm welch-cufft
	rm welch-cufft-openmp
