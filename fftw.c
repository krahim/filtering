#include<complex.h>
#include<fftw3.h>
#include"fftw.h"

void fft_r2c_(fftw_plan* p, int* n, double* data, 
             double complex* res) {
   
   int i, nc = *n/2 +1;

   *p = fftw_plan_dft_r2c_1d(*n, data, res, FFTW_ESTIMATE);

   fftw_execute(*p);

   for(i=nc; i < *n; i++) {
         res[i] = conj(res[*n - i]); 
   }
         
}
 

void fft_c2r_(fftw_plan* p, int* n, double complex* data, 
                double* res) {
   
   *p = fftw_plan_dft_c2r_1d(*n, data, res, FFTW_ESTIMATE);

   fftw_execute(*p);

}


void fft_c2c_(fftw_plan* p, int* n, double complex* data, 
             double complex* res, int* inverse) {
   int sign;
   
   if(*inverse == 1) {
      sign = FFTW_BACKWARD;
   } else {
      sign = FFTW_FORWARD;
   }
   
   *p = fftw_plan_dft_1d(*n, data, res, sign, FFTW_ESTIMATE);
   
   fftw_execute(*p);

}
