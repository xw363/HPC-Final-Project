/**
 * File: cufft.c
 * Description: Use CUDA's cufft library to compute discrete Fourier transform
                from real to real. x is truncated if n > nfft, and zero-padded
                if n < nfft. It is assumed that x points to an allocated memory
                block of size n, and xfft points to an allocated memory block
                of size nfft.
 * Author: Xiaojun Wu <xiaojun.wu@nyu.edu>
 */
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>
#include "welch.h"

welchStatus_t cufft(double *x, int n, double *xfft, int nfft)
{
    cufftHandle plan;                      /* The cufft plan */
    cufftDoubleReal *xPadded, *d_xPadded;  /* Zero-padded x on host and GPU */
    cufftDoubleComplex *xfftComplex;       /* xfft in complex numbers on host */
    cufftDoubleComplex *d_xfftComplex;     /* xfftComplex on GPU */
    int i;                                 /* For loop index */
    cudaError_t cudaStatus;
    cufftResult cufftStatus;
    welchStatus_t status;

    /* Pad x with 0 if n < nfft */
    status = padZero(x, n, &xPadded, nfft);
    if (status != WELCH_SUCCESS) {
        fprintf(stderr, "Error in cufft(): Failed to allocate memory "
                "on host.\n");

        return WELCH_FAILURE;
    }

    /* Initilize complex version of xfft */
    xfftComplex = (cufftDoubleComplex*) malloc((nfft / 2 + 1) *
                                               sizeof(cufftDoubleComplex));
    if (xfftComplex == NULL) {
        fprintf(stderr, "Error in cufft(): Failed to allocate memory "
                "on host.\n");

        free(xPadded);

        return WELCH_FAILURE;
    }

    /* Initialize arrays on GPU */
    cudaStatus = cudaMalloc((void**) &d_xPadded, nfft * sizeof(cufftDoubleReal));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Error in cufft(): Failed to allocate memory "
                "on GPU.\n");

        free(xPadded);
        free(xfftComplex);

        return WELCH_FAILURE;
    }

    cudaStatus = cudaMalloc((void**) &d_xfftComplex,
                            (nfft / 2 + 1) * sizeof(cufftDoubleComplex));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Error in cufft(): Failed to allocate memory "
                "on GPU.\n");

        free(xPadded);
        free(xfftComplex);
        cudaFree(d_xPadded);

        return WELCH_FAILURE;
    }

    /* Copy padded x to GPU memory */
    cudaStatus = cudaMemcpy(d_xPadded, xPadded, nfft * sizeof(cufftDoubleReal),
                            cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Error in cufft(): Failed to copy data to GPU.\n");

        free(xPadded);
        free(xfftComplex);
        cudaFree(d_xPadded);
        cudaFree(d_xfftComplex);

        return WELCH_FAILURE;
    }

    /* Set cufft plan */
    cufftStatus = cufftPlan1d(&plan, nfft, CUFFT_D2Z, 1);
    if (cufftStatus != CUFFT_SUCCESS) {
        fprintf(stderr, "Error in cufft(): Failed to get a CUFFT plan.\n");

        free(xPadded);
        free(xfftComplex);
        cudaFree(d_xPadded);
        cudaFree(d_xfftComplex);

        return WELCH_FAILURE;
    }

    /* Run cufft plan */
    cufftStatus = cufftExecD2Z(plan, d_xPadded, d_xfftComplex);
    if (cufftStatus != CUFFT_SUCCESS) {
        fprintf(stderr, "Error in cufft(): Failed to execute a CUFFT plan.\n");

        free(xPadded);
        free(xfftComplex);
        cudaFree(d_xPadded);
        cudaFree(d_xfftComplex);

        return WELCH_FAILURE;
    }

    /* Retrieve result from GPU */
    cudaStatus = cudaMemcpy(xfftComplex, d_xfftComplex,
                            nfft * sizeof(cufftDoubleReal),
                            cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess)  {
        fprintf(stderr, "Error in cufft(): Failed to copy data from GPU.\n");

        free(xPadded);
        free(xfftComplex);
        cudaFree(d_xPadded);
        cudaFree(d_xfftComplex);

        return WELCH_FAILURE;
    }

    /* Convert complex result to real */
    xfft[0] = xfftComplex[0].x;
    for (i = 1; i <= (nfft - 1) / 2; ++i) {
        xfft[2 * i - 1] = xfftComplex[i].x;
        xfft[2 * i] = xfftComplex[i].y;
    }
    if (nfft % 2 == 0) {
        xfft[nfft - 1] = xfftComplex[nfft / 2].x;
    }

    cufftDestroy(plan);
    free(xPadded);
    free(xfftComplex);
    cudaFree(d_xPadded);
    cudaFree(d_xfftComplex);

    return WELCH_SUCCESS;
}
