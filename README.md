#Final Project, High Performance Computing (NYU Spring 2015)

## Compile the programs
Check the followings before compiling the programs:
- fftw3 and CUDA toolkit are installed on your system.
- CUDA header directory is added to include path, and CUDA library directory is
added to linking and library path.
- The compilation flag `-lfftw3_omp` works. If not, change it to
`-lfftw3_threads` in `Makefile`.
Enter `make` in terminal to compile the programs.
