#include "delsx.hpp"
#include "signal.hpp"
//#include "fft.h"
#include <cmath>
#include <iostream>
#include <vector>
#define pi 3.1415926

static  struct objDelsx *mysx;

void delsxlInit()
{
    mysx = (struct objDelsx *) malloc (sizeof( struct objDelsx) );
    mysx->fs = 16000;
    mysx->wlen = 400;
    mysx->inc = 160 ;
    mysx->nfft = 1024 ;

    mysx->decinc = 3;
    mysx->NN = 4 ;
    mysx->MM = 2 ;

    mysx->myFFT = (struct objFFT*) malloc(sizeof(struct objFFT));
    fftInit(mysx->myFFT ,  mysx->nfft );

}

void delsxTerminate()
{
    fftTerminate(mysx->myFFT);
    free(mysx->myFFT);
    free(mysx);
}

static int crosszeros(float* y , int L)
{
    int zcr =0 ;
    for (int i = 0; i < L-1; ++i) {
        if (y[i]*y[i+1] < 0)
            zcr = zcr +1 ;
    }
    return zcr ;
}

static void win1fun( float* wnd )
{
    for (int i = 0; i <mysx->wlen ; ++i) {
        wnd[i] = 0.54 - 0.46* cos(2*i* pi/( mysx->wlen-1) ) ;
    }
    return ;
}

static void win2fun( vector<float > &w1 , vector<float > &w2)
{
    for (int j = 0; j < mysx->inc* mysx->MM ; ++j) {
        float tmp1 = (float) j / (mysx->inc* mysx->MM) ;
        w1[j] = (float) j / (mysx->inc* mysx->MM) ;
//        w2[j] =  1 -  w1[j] ;
    }


    for (int k = 0 ; k < mysx->inc* mysx->MM; ++k) {
        float tmp1 = (float) (mysx->inc* mysx->MM -k)/(mysx->inc* mysx->MM) ;
        w2[k] = (float) (mysx->inc* mysx->MM -k)/(mysx->inc* mysx->MM) ;
    }

    return ;
}


static void DetFeature( vector<float > x, vector<float >xlow , int num , vector<int >&Am , vector<float >&cz ,vector<float >&Ylow)
{
    std::vector<float>wnd(mysx->wlen);
    std::vector<float>y(mysx->nfft);
    std::vector<float>Yreal(mysx->nfft);
    std::vector<float>Yimag(mysx->nfft);
    std::vector<float>ylow(mysx->nfft);
    std::vector<float>Ylowreal(mysx->nfft);
    std::vector<float>Ylowimag(mysx->nfft);
    std::vector<float>A(mysx->nfft/2);
    std::vector<float>Alow(mysx->nfft/2);
    win1fun(wnd.data());

    for (int i = 0; i < num ; ++i) {
        int ind1 = mysx->inc *i;
        int ind2 = mysx->inc *i + mysx->wlen;

        for (int j = 0; j < mysx->wlen ; ++j) {
            y[j] = x[ind1 + j ]* wnd[j];
            ylow[j] = xlow[ind1 + j ] * wnd[j] ;
//            y[j] = x[ind1 + j ];
//            ylow[j] = xlow[ind1 + j ]  ;
        }
        for (int j = mysx->wlen ; j < mysx->nfft; ++j) {
            y[j] = 0 ;
            ylow[j] = 0;
        }

        fftComputeOnce(mysx->myFFT, y.data(), Yreal.data() , Yimag.data());
        fftComputeOnce(mysx->myFFT, ylow.data(), Ylowreal.data() , Ylowimag.data());
        for (int k = 0; k < mysx->nfft/2; ++k) {
            A[k] = sqrt( Yreal[k]* Yreal[k] + Yimag[k]*Yimag[k]);
            Alow[ k] = sqrt( Ylowreal[k]* Ylowreal[k] + Ylowimag[k]*Ylowimag[k]);
        }

        for (int l = 0; l < mysx->nfft/2; ++l) {
            int ind3 = i* mysx->nfft/2 + l;
            Ylow[ind3] =  Alow[l];
        }

        int numfs1 = floor(1800.0/mysx->fs*mysx->nfft);
        int numfs2 = floor(3000.0/mysx->fs*mysx->nfft);
        float tmpA1 =0.0;
        float tmpA2 = 0.0;
        for (int k1 =  numfs2 ; k1 <512 ; ++k1) {
             tmpA1 = tmpA1 + A[k1];
        }

        for (int k2 = 0; k2 < numfs1 ; ++k2) {
             tmpA2 =  tmpA2 + A[k2];
        }

        if (tmpA2 < tmpA1)
            Am[i] = 1;
        else
            Am[i] = 0;

        cz[i] = crosszeros(y.data() , mysx->wlen);
    }
//    int tmp1 =1;
    return ;
}

static void DetTime( vector<int>Am, vector<float >cz, int num , vector<int>&t1 , vector<int>&t2, int &numseg)
{
    int maxcz = 150 ;
    int minlen = 9 ;
    int status = 0;
    int segn = 0 ;
    int count = 0 ;

    for (int j = 0; j < num; ++j) {
        t1[j] = 0;
        t2[j] = 0;
    }

    for (int i = 0; i < num ; ++i) {
        switch (status){
            case 0 :
            case 1 :
                if (Am[i]  && cz[i] > maxcz)
                {
                    t1[segn] = i ;
                    status = 2 ;
                    count = count +1 ;
                }
                else if (Am[i]  || cz[i] > maxcz){
                    status =1;
                    count = count +1 ;
                }
                else{
                    status =0 ;
                    count = 0 ;
                    t1[segn] = 0 ;
                    t2[segn] = 0 ;
                }
                break ;
            case 2:
                if (Am[i]  && cz[i] > maxcz){
                    count = count +1 ;
                }
                else{
                    if  (count < minlen){
                        status = 0 ;
                        count = 0 ;
                    }
                    else{
                        status =3 ;
                        t2[segn] = t1[segn] + count ;
                    }
                }
                break ;
            case 3:
                segn =segn + 1;
                status = 0 ;
                count = 0 ;
                t1[segn] = 0;
                t2[segn] = 0;
                break ;
        }
    }

    for (int j = 0; j < num; ++j) {
        if (t1[j] != 0 ){
            numseg = j+1 ;
        }

    }
    if (t2[numseg-1] == 0){
        t2[numseg-1] = num;
    }

    for (int k = 0; k < numseg; ++k) {
        t2[k] = t2[k]- mysx->decinc;
    }

    return ;

}

static void DetFre(vector<float>Ylow , vector<int >t1 , vector<int >t2 , int numseg ,  vector<float >&f0)  //Ylow: 512* numseg
{
    int ind1 = floor( 2000.0/ mysx->fs *mysx->nfft );
    int ind2 = floor( 7500.0 / mysx->fs *mysx->nfft );
    std::vector<float>EA(mysx->nfft/2) ;
    std::vector<float>tmpyy(512) ;

    for (int i = 0; i < numseg; ++i) {
         int T1 = t1[i] ;
         int T2 = t2[i] ;

        for (int l = 0; l < mysx->nfft/2; ++l) {
            EA[l] = 0 ;
        }
         // 2
        for (int n = T1; n < T2+1; ++n) {
            int Tseg = n* mysx->nfft/2 ;
            for (int j = 0; j <  mysx->nfft/2; ++j) {
                tmpyy[j] =  Ylow[Tseg + j] ;
                EA[j] = EA[j] + Ylow[Tseg + j] ;
            }
        }

        // 3
        float tmpAE = 0;
        int indmax = ind1 ;
        for (int k =  ind1; k < ind2; ++k) {
            if(EA[k] > tmpAE)
            {
                tmpAE = EA[k] ;
                indmax = k +1 ;
            }
        }

        // 4
        f0[i]  =  (float)indmax / mysx->nfft * mysx->fs;

        }

    return ;
}

void BandStop(float fs ,float fc,  vector<float > &a ,vector<float > &b )
{
    float beta =0.5 ;
    float fb = 600 ;
    float Q = beta *fc / fb ;
    float K = tan(pi *fc /fs);

    b[0] = Q*(1+ K*K) / (K*K*Q + K + Q);
    b[1] = 2*Q*(K*K  - 1)/  (K*K*Q + K + Q);
    b[2] = Q*(1+ K*K) / (K*K*Q + K + Q);
    a[0] = 1;
    a[1] = 2*Q*(K*K - 1)/ (K*K*Q + K + Q);
    a[2] = (K*K*Q - K + Q) / (K*K*Q + K + Q);

    return ;
}


static void segfilter( vector<float>x1 ,int fs , float f0 , vector<float > w1 , vector<float >w2 , vector<float>&x3 )
{
    int Lin = x1.size();
    vector<float>a(3);
    vector<float>b(3);
    std::vector<float >xfilter(Lin);
    std::vector<float >x2(Lin);

    BandStop( fs , f0,  a,  b );
//     a[0] = 1.0000 ;
//     a[1] = -0.7561 ;
//     a[2] =  0.4757 ;
//     b[0] = 0.7405 ;
//     b[1] = -0.7561 ;
//     b[2] =  0.7352 ;

    filter(2 , a.data(), b.data(), Lin-1 ,x1.data(), xfilter.data() );

    int Lout =  Lin - mysx->inc * mysx->MM * 2;

    for (int i = 0; i < Lout; ++i) {
        x2[i] = xfilter[i+ mysx->inc*mysx->MM];
        x3[i] = x2[i];
    }

    for (int i = 0; i < mysx->inc*2 ; ++i) {
        x3[i] =x2[i] * w1[i];
    }

    for (int i = Lout - mysx->inc*2 ; i < Lout; ++i) {
        int indw2 =i-(Lout - mysx->inc*2);
        x3[i] = x2[i] * w2[indw2];
    }

//    float testw2 =x2[5600 - Lout + mysx->inc*2];
//    float test1 =x2[5600];
//    float test2 =x3[5600];

    return ;
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

    for (i = ord + 1; i<np + 1; i++)
    {
        y[i] = 0.0;
        for (j = 0; j<ord + 1; j++)
            y[i] = y[i] + b[j] * x[i - j];
        for (j = 0; j<ord; j++)
            y[i] = y[i] - a[j + 1] * y[i - j - 1];
    }

    return;
}


void delsxprocess(short *in ,  int L,  vector<short > &out)
{
    int num =  floor((L- mysx->wlen)/mysx->inc) -1 ;
    std::vector<float> x(L) ;
    std::vector<float> xlow(L) ;
    std::vector<int >Am(num);
    std::vector<float>cz(num);
    std::vector<float>Ylow(mysx->nfft/2 *num);

    for (int i = 0; i < L; ++i) {
        float tmp = (float) in[i];
        x[i] = tmp /32768 ;
    }
    float maxValue = *max_element(x.begin(),x.end());
    for (int i = 0; i < L; ++i) {
        x[i] = x[i] / maxValue ;
    }

    // lowfilter
    float alow[3] = {1.0000 ,  -0.0336  ,  0.3166};
    float blow[3] = { 0.5224 ,   0.4517  ,  0.1694};
    filter(2 , alow, blow, L-1 ,x.data(), xlow.data() );

    DetFeature( x, xlow , num ,Am , cz , Ylow) ;
    std::vector<int>t1(num);
    std::vector<int>t2(num);
    std::vector<float >f0(num);
    int numseg =0 ;
    DetTime(Am , cz , num , t1 ,t2 , numseg);
    DetFre(Ylow, t1, t2 , numseg ,f0) ;

    std::vector<float >w1(mysx->inc * mysx->MM );
    std::vector<float >w2(mysx->inc * mysx->MM );
    win2fun( w1 , w2);
    for (int j = 0; j < numseg ; ++j) {
        t1[j] = t1[j] +1 ;
        t2[j] = t2[j]+ 1 ;
        int indx1 = t1[j]* mysx->inc  ;
        int indx2 = t2[j]* mysx->inc +mysx->wlen -1 ;
        int Lx =  indx2 -indx1 +1 ;

        int indNN1 = (t1[j] - mysx->NN)* mysx->inc  ;
        int indNN2 = (t2[j] + mysx->NN)* mysx->inc + mysx->wlen -1 ;
        int LNN = indNN2 - indNN1+1 ;
        std::vector<float >segframe(LNN);
        std::vector<float >Segframe(LNN - mysx->inc * mysx->MM*2 );
        for (int i = 0; i < LNN ; ++i) {
            segframe[i] = xlow[i+  indNN1];
        }
        segfilter(segframe, mysx->fs , f0[j] ,w1 ,w2,  Segframe);

        int indMM1 = (t1[j] - mysx->MM)* mysx->inc  ;
        int indMM2 = (t2[j] + mysx->MM)* mysx->inc + mysx->wlen -1 ;
        int LMM = indMM2 - indMM1 +1 ;
        std::vector<float >xmid( LMM );
        std::vector<float >Xmid( LMM );
        std::vector<float >Seg( LMM );
        for (int i = 0; i < mysx->inc * mysx->MM ; ++i) {
            xmid[i] = x[i + indMM1] ;
            Xmid[i] = x[i + indMM1] *w2[i];
        }
        for (int i = mysx->inc * mysx->MM; i < mysx->inc * mysx->MM+ Lx; ++i) {
            xmid[i] = 0;
            Xmid[i] = 0 ;
        }
        for (int i = mysx->inc * mysx->MM+ Lx ;  i < LMM ; ++i) {
            int tmpw1 = i - ( mysx->inc * mysx->MM+ Lx )  ;
            int tmpxlow = indx2 + tmpw1 + 1;
            xmid[i] = x[tmpxlow] ;
            Xmid[i] = x[tmpxlow] *w1[tmpw1];
        }
        float test1 = Xmid[5429];
        for (int k = 0; k < LMM; ++k) {
            Seg[k] = Xmid[k] + Segframe[k];
            x[indMM1+k] =  Seg[k] ;
//            out[k] = (short)( Seg[k]*32768*maxValue );
        }

        int Tind =5350 ;
//        float testXmid= Xmid[Tind];
//        float testSeg = Segframe[Tind];
//        float test2 = Seg[Tind];
//        float test3= x[indMM1+ Tind];


    }

    for (int j = 0; j <L ; ++j)
    {
        out[j] = (short)( x[j]*32768*maxValue );
    }

    return ;

}
