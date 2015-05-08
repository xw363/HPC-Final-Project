/**
 * File: welch.h
 * Description: Defines macros and functions for Welch an FFT implementations
 *
 * Author: Xiaojun Wu <xiaojun.wu@nyu.edu>
 */
#ifndef WELCH_H
#define WELCH_H

#include <stddef.h>

/**
 * Function return status
 */
typedef enum {
    WELCH_SUCCESS = 0,
    WELCH_FAILURE = 1
} welchStatus_t;

/**
 * The Welch method for real signals
 * signal - input signal
 * Pxx - spectral density estimate.
 * frequency - frequencies where the spectral density is estimated.
 * samplingFrequency - sampling frequency of the signal
 * lenSignal - length of the signal array / number of samples
 * lenSegment - length of a single segment of signals
 * lenOverlap - length of overlap for two consecutive segments
 * lenPxx - length of spectral density estimate, determined by this function
 * windowType - type of window function to apply
 *              (only rectangular window can be used at this time)
 * fftType - type of FFT implementation to use
 * nfft - number of points to do FFT
 *
 * Returns a welchStatus_t
 */
welchStatus_t welch(double *signal, double **Pxx, double **frequency,
                    double samplingFrequency, int lenSignal, int lenSegment,
                    int lenOverlap, int *lenPxx, char *windowType,
                    char *fftType, int nfft);

/**
 * FFT routine wrappers
 * x - input data
 * n - length of input array
 * xfft - Fourier transformed data
 * nfft - length of FFT
 * useOpenMP - enable OpenMP in fftw. 1 for yes, 0 for no
 *
 * Returns a welchStatus_t
 */
welchStatus_t fftw(double *x, int n, double *xfft, int nfft, int useOpenMP);
welchStatus_t cufft(double *x, int n, double *xfft, int nfft);

/* Utility functions */

/**
 * Get window function
 * windowType - type of the desired window funtion
 * window - returned array representing the window function
 * lenWindow - length of the window
 *
 * Returns a welchStatus_t
 */
welchStatus_t getWindow(char *windowType, double *window, int lenWindow);

/**
 * Pad an array with 0. If n = nPadded, xPadded is simply a copy of x.
 * x - array to be padded
 * n - size of x
 * xPadded - zero-padded x
 * nPadded - size of xPadded
 *
 * Returns a welchStatus_t
 */
welchStatus_t padZero(double *x, int n, double **xPadded, int nPadded);

#endif
