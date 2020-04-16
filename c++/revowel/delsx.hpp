//
// Created by zyc on 2020/3/31.
//

#ifndef _DELSX_HPP_
#define _DELSX_HPP_
#include "fft.h"
#include<algorithm>
//#include "signal.hpp"
#include <vector>
using namespace std ;
typedef struct objDelsx
{
    int fs;
    int inc;
    int wlen;
    int nfft;

    int decinc ;
    int NN ;
    int MM ;

    struct objFFT* myFFT;

}objDelsx;


void delsxlInit();
void delsxTerminate();
//void delsxprocess(short *in ,  int L,  short *out);
//void delsxprocess( vector<short >in ,  int L,  vector<short > &out);
void delsxprocess( short *in ,  int L,  vector<short > &out);




#endif


