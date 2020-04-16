#ifndef FFT_H
#define FFT_H

#include <math.h>

//#include "macro.h"
//#include "../hardware.h"
//#include "../parameters.h"
#include "dynamicMemory.h"

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 * Structure                                                                   *
 ******************************************************************************/

struct objFFT
{

    // +-------------------------------------------------------------------+
    // | Parameters                                                        |
    // +-------------------------------------------------------------------+

    // Number of points in the FFT (must be a power of 2)
    unsigned int FFT_SIZE;

    // Number of levels in the FFT (number of poinst = 2^(number of levels)
    unsigned int FFT_NBLEVELS;

    // Half the number of points in the FFT
    unsigned int FFT_HALFSIZE;

    // Inverse of the number of points in the FFT
    float FFT_SIZE_INV;

    // Number of elements in the array for SIMD group optimization
    unsigned int FFT_SIMD_GROUP;

    // Number of elements in the array for SIMD individual optimization
    unsigned int FFT_SIMD_INDIVIDUAL;

    // +-------------------------------------------------------------------+
    // | Variables                                                         |
    // +-------------------------------------------------------------------+

    // Define the FFT factor
    // The FFT factor is Wn(r) = exp(-j*2*pi*r/N) for r = 0 ... (N/2 - 1)
    // It is precomputed when the FFT module is initialized. Once this is done,
    // the results will be used for each FFT to be computed
    //
    // Size: [FFT_HALFSIZE]
    float* WnReal;
    float* WnImag;
    // Size: [FFT_SIMD_GROUP*4]
    float* simdWnReal;
    float* simdWnImag;

    // Define working array
    // This array is used during the FFT computations. It is set as a global
    // variable to avoid stacking and unstacking all these values when there is
    // a function call. This will speed up the computations.
    //
    // Size: [FFT_SIZE]
    float* workingArrayReal;
    float* workingArrayImag;

    // Define arrays to compute 2 FFTs with only one FFT for real signals
    //
    // Size: [FFT_SIZE]
    float* fftTwiceReal;
    float* fftTwiceRealFlipped;
    float* fftTwiceImag;
    float* fftTwiceImagFlipped;

    // Define array to compute 1 FFT for real signals
    //
    // Size: [FFT_SIZE]
    float* emptyArray;

    // Define array to trash unused result from IFFT
    //
    // Size: [FFT_SIZE]
    float* trashArray;

    // Define array to store the reverse bit order addressing
    //
    // Size: [FFT_SIZE]
    unsigned int* revBitOrderArray;

    // Define arrays to store the a indexes for SIMD optimization
    //
    // Size: [FFT_SIMD_GROUP]
    float** simdARealGroups;
    float** simdAImagGroups;
    float** simdBRealGroups;
    float** simdBImagGroups;
    float** simdRRealGroups;
    float** simdRImagGroups;
    // Size: [FFT_SIMD_INDIVIDUAL]
    float** simdAIndividual;
    float** simdBIndividual;

};

/*******************************************************************************
 * Prototypes                                                                  *
 ******************************************************************************/

void fftInit(struct objFFT* myFFT, unsigned int size);

void fftTerminate(struct objFFT* myFFT);

void fftCompute(struct objFFT* myFFT, float* sourceArrayReal, float* sourceArrayImag, float* destArrayReal, float* destArrayImag);

void fftComputeOnce(struct objFFT* myFFT, float *sourceArray, float *destArrayReal, float *destArrayImag);

void fftComputeTwice(struct objFFT* myFFT, float *sourceArray1, float *sourceArray2, float *destArray1Real, float *destArray1Imag, float *destArray2Real, float *destArray2Imag);

void ifftCompute(struct objFFT* myFFT, float *sourceArrayReal, float *sourceArrayImag, float *destArrayReal, float *destArrayImag);

void ifftComputeOnce(struct objFFT* myFFT, float *sourceArray1Real, float *sourceArray1Imag, float *destArray1);

void ifftComputeTwice(struct objFFT* myFFT, float *sourceArray1Real, float *sourceArray1Imag, float *sourceArray2Real, float *sourceArray2Imag, float *destArray1, float *destArray2);

#ifdef __cplusplus
}
#endif

#endif
