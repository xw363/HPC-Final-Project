/**
 * File: welch-cufft-openmp.c
 * Description: Test the welch() function with cufft and openmp library.
 *
 * Author: Xiaojun Wu <xiaojun.wu@nyu.edu>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include "welch.h"

#define PI 3.1415926535897932384626
#define N 16384

int main(int argc, char *argv[])
{
    double *signal, *Pxx, *frequency;
    int lenSignal, lenSegment, lenOverlap, lenPxx, nfft, samplingFrequency;
    int tid, i;
    welchStatus_t status;
    struct timeval tic, toc;  /* Start and finish time */
    double total_time;

    /* Set up variables */
    lenSignal = N;
    lenSegment = N / 4;
    lenOverlap = N / 8;
    samplingFrequency = 1000;
    nfft = N / 2;

    /* Generate input signal */
#pragma omp parallel private(tid, signal, Pxx, frequency, lenPxx, tic, toc)
{
    tid = omp_get_thread_num();
    signal = malloc(lenSignal * sizeof(double));

    if (signal != NULL) {
        for (i = 0; i < lenSignal; ++i) {
            signal[i] = 5 * sin(2 * PI * i / N);
        }

        /* Run the algorithm */
        gettimeofday(&tic, NULL);
        status = welch(signal, &Pxx, &frequency, samplingFrequency, lenSignal,
                       lenSegment, lenOverlap, &lenPxx, "rectangular", "cufft",
                       nfft);
        gettimeofday(&toc, NULL);

        /* Print time spent on the Welch method */
        if (status == WELCH_SUCCESS) {
            total_time = toc.tv_sec - tic.tv_sec
                         + (toc.tv_usec - tic.tv_usec) / 1e6;
            printf("Thread %d: Welch method completed in %.8f seconds.\n",
                   tid, total_time);
        } else {
            printf("Thread %d: Welch method failed.\n", tid);
        }

        free(signal);
        free(Pxx);
        free(frequency);
    } else {
        printf("Thread %d: Failed to allocate memory for signals.\n", tid);
    }
}

    return EXIT_SUCCESS;
}
