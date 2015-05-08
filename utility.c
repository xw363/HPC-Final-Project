/**
 * File: utility.c
 * Description: Utility functions defined in welch.h
 *
 * Author: Xiaojun Wu <xiaojun.wu@nyu.edu>
 */
#include <stdio.h>
#include <stdlib.h>
#include "welch.h"

welchStatus_t getWindow(char *windowType, double *window, int lenWindow)
{
    /* TODO: Implement this function */

    return WELCH_SUCCESS;
}

welchStatus_t padZero(double *x, int n, double **xPadded, int nPadded)
{
    int i;

    if (n > nPadded) {
        fprintf(stderr, "Error in padZero(): The original array has larger"
                "size than the array to pad.\n");

        return WELCH_FAILURE;
    }

    *xPadded = (double*) malloc(nPadded * sizeof(double));
    if (*xPadded == NULL) {
        fprintf(stderr, "Error in padZero(): Failed to allocate memory\n");

        return WELCH_FAILURE;
    }

    /* Copy x */
    for (i = 0; i < n; ++i) {
        (*xPadded)[i] = x[i];
    }

    /* Pad zeros */
    for (; i < nPadded; ++i) {
        (*xPadded)[i] = 0.0;
    }

    return WELCH_SUCCESS;
}
