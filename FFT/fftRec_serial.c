#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <sys/time.h>
#define PI 3.14159265358979323846

/* Standard NON-Parallel FFT Implementation. For comparison purposes. */

double complex * fftRec(double complex *, int, int);
int main(int argc, char ** argv)
{
  int i;
  double complex * arr;
  double complex * out;

 int numel=16384;
 double z;

 arr = (double complex *) malloc(sizeof(double complex)*numel);
 out = (double complex *) malloc(sizeof(double complex)*numel);
   
  /* Random Number Generation. */
  for (i = 0; i < numel; i++)
    {
       z=rand() % 40 - 20;
       arr[i]=z;
    }

struct timeval  tv1, tv2;
gettimeofday(&tv1, NULL);

out = fftRec(arr, 1, numel);

gettimeofday(&tv2, NULL);


free(arr);
free(out);

printf ("Total Serial time = %f seconds\n",
(double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec));

return 0;
}


double complex * fftRec(double complex * x, int strd, int n)
{
  /* Recursion Condition. */
  if(n == 1)
  {
    double complex * out = malloc(sizeof(double complex)*n);
    out[0] = x[0];
    return out;
  }
  else
  {
    double complex * out1, * out2;
    int i;
    /* Recursive Calls. */
    out1 = fftRec(x, 2*strd, n/2);
    out2 = fftRec(x+strd, 2*strd, n/2);

    double complex * out = malloc(sizeof(double complex)*n);

    for(i = 0; i < n/2; i++)
    {
    /* FFT Computations. */
      out[i] = out1[i] + cexp((2*PI*I*i)/n)*out2[i];
      out[i+n/2] = out1[i] + cexp((2*PI*I*(i+n/2))/n)*out2[i];
    }
    free(out1);
    free(out2);
    return out;
  }
}