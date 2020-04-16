#include "signal.hpp"

static void FindMinMax(float A[],int size,  float *Pmin, float *Pmax)
{
    float min=A[0];
    float max=A[0];
    int  i ;
    for(i=1;i<size;i++)
    {
        if(A[i]>max)
            max=A[i];
        if(A[i]<min)
            min=A[i];
    }
    *Pmin = min;
    *Pmax = max;
}

static void filter(int ord, float *a, float *b, int np, float *x, float *y)
{
    int i, j;
    y[0] = b[0] * x[0];
    for (i = 1; i<ord + 1; i++)
    {
        y[i] = 0.0;
        for (j = 0; j<i + 1; j++)
            y[i] = y[i] + b[j] * x[i - j];
        for (j = 0; j<i; j++)
            y[i] = y[i] - a[j + 1] * y[i - j - 1];
    }
    /* end of initial part */
    for (i = ord + 1; i<np + 1; i++)
    {
        y[i] = 0.0;
        for (j = 0; j<ord + 1; j++)
            y[i] = y[i] + b[j] * x[i - j];
        for (j = 0; j<ord; j++)
            y[i] = y[i] - a[j + 1] * y[i - j - 1];
    }
    return;
} /* end of filter */


