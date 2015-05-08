/**
 * File: welch.c
 * Description: Implements the Welch method for signals of real numbers, which
 *              calls user-specified FFT routines.
 *
 * Author: Xiaojun Wu <xiaojun.wu@nyu.edu>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "welch.h"

#define FFTW 0
#define FFTW_OPENMP 1
#define CUFFT 2

welchStatus_t welch(double *signal, double **Pxx, double **frequency,
                    double samplingFrequency, int lenSignal, int lenSegment,
                    int lenOverlap, int *lenPxx, char *windowType,
                    char *fftType, int nfft)
{
    double *signalfft;          /* FFT of windowed signal */
    double *PxxInternal;        /* All computation of Pxx is done to this
                                   variable so that Pxx is not touched if some
                                   error occurs. */
    double *frequencyInternal;  /* Similar purpose, but for frequency */
    double lenPxxInternal;      /* Similar purpose, but for lenPxx */
    int numSegment;             /* Number of segments */
    double scale;               /* Scale for Pxx */
    double *window;             /* Array representing the window function */
    double *windowedSignal;     /* Signal convolved with window */
    double normSquared;         /* Placeholder for squared norm of an array */
    int fftCall;                /* Type of FFT implementation to call */
    int i, j;                   /* Loop indices */
    int status;                 /* Function status */

    /* Check inputs */
    if (samplingFrequency <= 0) {
        fprintf(stderr, "Sampling frequency of signal must be positive.\n");

        return WELCH_FAILURE;
    }

    if (lenSignal <= 0) {
        fprintf(stderr, "Length of signal must be positive.\n");

        return WELCH_FAILURE;
    }

    if (lenSegment <= 0) {
        fprintf(stderr, "Length of segment must be positive.\n");

        return WELCH_FAILURE;
    }

    if (lenOverlap < 0) {
        fprintf(stderr, "Number of overlapping points must be non-negative\n");

        return WELCH_FAILURE;
    }

    if (lenSignal < lenSegment) {
        fprintf(stderr, "Length of segment must be smaller than length "
                "of signal.\n");

        return WELCH_FAILURE;
    }

    if (lenSegment < lenOverlap) {
        fprintf(stderr, "Length of overlap must be smaller than length "
                "of segment.\n");

        return WELCH_FAILURE;
    }

    if ((lenSignal - lenOverlap) % (lenSegment - lenOverlap) != 0) {
        fprintf(stderr, "Unable to determine integral number of segments.\n");

        return WELCH_FAILURE;
    }

    /* Get window function */
    if (strcmp(windowType, "rectangular") == 0) {
        window = (double*) malloc(lenSignal * sizeof(double));

        if (window == NULL) {
            fprintf(stderr, "Failed to allocate memory in Welch method.\n");

            return WELCH_FAILURE;
        }

        for (i = 0; i < lenSignal; ++i) {
            window[i] = 1.0;
        }
    } else {
        fprintf(stderr, "Unrecoginzed type of window function.\n");

        return WELCH_FAILURE;
    }

    if (strcmp(fftType, "fftw") == 0) {
        fftCall = FFTW;
    } else if (strcmp(fftType, "fftw_openmp") == 0) {
        fftCall = FFTW_OPENMP;
    } else if (strcmp(fftType, "cufft") == 0) {
        fftCall = CUFFT;
    } else {
        fprintf(stderr, "Error in welch(): Unrecoginzed FFT implementation.\n");

        free(window);

        return WELCH_FAILURE;
    }

    if (nfft <= 0) {
        fprintf(stderr, "Number of FFT points must be positive.\n");

        free(window);

        return WELCH_FAILURE;
    }

    /* Initialize variables */
    signalfft = (double*) malloc(nfft * sizeof(double));
    if (signalfft == NULL) {
        fprintf(stderr, "Failed to allocate memory in welch().\n");

        free(window);

        return WELCH_FAILURE;
    }

    if (nfft % 2 == 0) {
        lenPxxInternal = nfft / 2 + 1;
    } else {
        lenPxxInternal = (nfft + 1) / 2;
    }

    PxxInternal = NULL;
    PxxInternal = (double*) malloc(lenPxxInternal * sizeof(double));
    if (PxxInternal == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for Pxx in "
                        "welch(). Pxx is not modified.\n");

        free(signalfft);
        free(window);

        return WELCH_FAILURE;
    }

    for (i = 0; i < lenPxxInternal; ++i) {
        PxxInternal[i] = 0.0;
    }

    frequencyInternal = NULL;
    frequencyInternal = (double*) malloc(lenPxxInternal * sizeof(double));
    if (frequencyInternal == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for frequencies in "
                        "welch().\n");

        free(signalfft);
        free(window);

        return WELCH_FAILURE;
    }

    windowedSignal = (double*) malloc(lenSegment * sizeof(double));
    if (windowedSignal == NULL) {
        fprintf(stderr, "Failed to allocate memory in welch().\n");

        free(signalfft);
        free(PxxInternal);
        free(frequencyInternal);
        free(window);

        return WELCH_FAILURE;
    }

    /* Compute scale */
    normSquared = 0.0;
    for (i = 0; i < lenSegment; ++i) {
        normSquared += window[i] * window[i];
    }
    scale = 1.0 / (samplingFrequency * normSquared);

    /* Compute FFT of each segment, and add squared sums to Pxx */
    for (i = 0; i + lenSegment <= lenSignal; i += lenSegment - lenOverlap) {
        /* Convolve the current segment with window */
        for (j = 0; j < lenSegment; ++j) {
            windowedSignal[j] = signal[i + j] * window[j];
        }

        /* Compute FFT of the convolved signal */
        if (fftCall == FFTW) {
            status = fftw(windowedSignal, lenSegment, signalfft, nfft, 0);
        } else if (fftCall == FFTW_OPENMP) {
            status = fftw(windowedSignal, lenSegment, signalfft, nfft, 1);
        } else {
            status = cufft(windowedSignal, lenSegment, signalfft, nfft);
        }

        if (status == WELCH_FAILURE) {
            free(signalfft);
            free(PxxInternal);
            free(frequencyInternal);
            free(window);
            free(windowedSignal);

            return WELCH_FAILURE;
        }

        /* Add the squared magnitude of FFT to Pxx */
        PxxInternal[0] += signalfft[0] * signalfft[0];
        for (j = 1; j < lenPxxInternal; ++j) {
            if (nfft % 2 == 0 && j == lenPxxInternal - 1) {
                /* If nfft is even, the last term is real does not have a
                 * corresponding imaginary term */
                PxxInternal[j] += signalfft[nfft - 1] * signalfft[nfft - 1];
            } else {
                PxxInternal[j] += signalfft[2 * j - 1] * signalfft[2 * j - 1]
                                  + signalfft[2 * j] * signalfft[2 * j];
            }
        }
    }

    /* Scale Pxx and average it over number of segments */
    numSegment = (lenSignal - lenOverlap) / (lenSegment - lenOverlap);
    for (i = 0; i < lenPxxInternal; ++i) {
        if (i == 0 || i == lenPxxInternal - 1) {
            PxxInternal[i] *= scale / numSegment;
        } else {
            PxxInternal[i] *= scale * 2 / numSegment;
        }
    }

    /* Get frequencies */
    for (i = 0; i < lenPxxInternal; ++i) {
        frequencyInternal[i] = i * samplingFrequency / nfft;
    }

    free(signalfft);
    free(window);
    free(windowedSignal);

    /* Return Pxx and its length */
    *Pxx = PxxInternal;
    *frequency = frequencyInternal;
    *lenPxx = lenPxxInternal;

    return WELCH_SUCCESS;
}
