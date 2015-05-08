/**
 * File: fftw.c
 * Description: Use the fftw3 library to compute discrete Fourier transform
                from real to real. x is truncated if n > nfft, and zero-padded
                if n < nfft. It is assumed that x points to an allocated memory
                block of size n, and xfft points to an allocated memory block
                of size nfft.
 * Author: Xiaojun Wu <xiaojun.wu@nyu.edu>
 */
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <omp.h>
#include "welch.h"

welchStatus_t fftw(double *x, int n, double *xfft, int nfft, int useOpenMP)
{
    fftw_plan plan;              /* The FFT plan */
    double *xPadded;            /* Zero-padded x, if necessary */
    fftw_complex *xfftComplex;  /* xfft in complex numbers */
    welchStatus_t status;
    int i;                      /* Index of for loops */

    /* Pad x with 0 if n < nfft */
    status = padZero(x, n, &xPadded, nfft);
    if (status != WELCH_SUCCESS) {
        fprintf(stderr, "Failed to allocate memory in fftw()\n");

        return WELCH_FAILURE;
    }

    /* Initialize the complex version of xfft */
    xfftComplex = (fftw_complex*) fftw_malloc(nfft * sizeof(fftw_complex));
    if (xfftComplex == NULL) {
        fprintf(stderr, "Failed to allocate memory in fftw()\n");

        free(xPadded);

        return WELCH_FAILURE;
    }

    /* Set and run a fftw plan */
    if (useOpenMP) {
        fftw_init_threads();
        fftw_plan_with_nthreads(omp_get_max_threads());
    }
    plan = fftw_plan_dft_r2c_1d(nfft, xPadded, xfftComplex, FFTW_ESTIMATE);
    fftw_execute(plan);

    /* Convert complex result to real */
    xfft[0] = xfftComplex[0][0];
    for (i = 1; i <= (nfft - 1) / 2; ++i) {
        xfft[2 * i - 1] = xfftComplex[i][0];
        xfft[2 * i] = xfftComplex[i][1];
    }
    if (nfft % 2 == 0) {
        xfft[nfft - 1] = xfftComplex[nfft / 2][0];
    }

    fftw_destroy_plan(plan);
    free(xPadded);
    fftw_free(xfftComplex);

    return WELCH_SUCCESS;
}
