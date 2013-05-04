#ifndef __FFTW_H__
#define __FFTW_H__

#include<fftw3.h>

/* real to complex forward */
void fft_r2c_(fftw_plan* p, int* n, double* data, 
                double complex* res);

void fft_c2r_(fftw_plan* p, int* n, double complex* data, 
             double* res);

void fft_c2c_(fftw_plan* p, int* n, double complex* data, 
             double complex* res, int* inverse);
#endif
