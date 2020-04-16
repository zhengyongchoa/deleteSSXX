import numpy as np
import librosa
import math
import scipy
from librosa.output import write_wav
from scipy.signal import  butter, lfilter
from scipy.signal import get_window
import matplotlib.pyplot as plt

def crosszeros(x):
    wlen = len(x)
    zcr = 0
    for j in range(wlen-1):
         if x[j] *x[j+1] < 0 :
             zcr = zcr+1
    return  zcr

def DetFeature(Xframe,xsfft , fs=16000 ,nfft=1024):
    num  = Xframe.shape[1] -2
    numfs1 = np.floor(1800/fs*nfft).astype(int)
    numfs2 = np.floor(3000/fs*nfft ).astype(int)
    Am = np.zeros(num)
    cz = np.zeros(num)
    for i in range(num):
        A = np.abs(xsfft[:,i])
        x = Xframe[:,i]
        wnd = get_window('hamming', int(400))
        x= x*wnd
        X = np.fft.fft(x,n=1024)
        Ax = np.abs(X)
        # plt.plot(Ax[0:513])
        # plt.show()
        tmp1 = np.sum(Ax[0:numfs1])
        tmp2 = np.sum( Ax[numfs2: 513 ] )
        if tmp2 > tmp1 :
            Am[i] = 1
        else:
            Am[i] = 0
        cz[i] = crosszeros(Xframe[:,i])

    return cz ,Am


def DetTime(cz ,Am):
    maxcz = 150
    minlen = 9
    status = 0
    segn = 0
    count = 0
    num = cz.size
    t1 = np.zeros(num).astype(int)
    t2 = np.zeros(num).astype(int)
    for i in range(num):
        if status==0 or status ==1:
            if Am[i] and cz[i] > maxcz :
                t1[segn] = i
                status = 2
                count = count+1
            elif  Am[i] or cz[i] > maxcz :
                status = 1
                count = count + 1

            else:
                status =0
                count = 0
                t1[segn] = 0
                t2[segn] = 0

        elif status==2:
            if Am[i] and cz[i] > maxcz :
                count = count + 1
            else:

                if count < minlen :
                    status = 0
                    count = 0
                else:
                    status = 3
                    t2[segn] = t1[segn] + count

        else:
            segn = segn +1
            status = 0
            count =0
            t1[segn] = 0
            t2[segn] = 0

    numseg = 0
    for j in range(num):
        if t1[j] != 0:
            numseg =j +1
    if t2[numseg-1] ==0:
        t2[numseg-1] = num
    t1 = t1[0: numseg]
    t2 = t2[0: numseg]
    decinc = 3
    for j  in range(numseg):
        t2[j] = t2[j] - decinc
    return t1 ,t2

def DetFre( x ,t1 ,t2 ,fs1=2000 ,fs2 =7500 , fs =16000, nfft = 1024):
    xframe = librosa.util.frame(x, frame_length=400, hop_length=160)
    num = xframe.shape[1] - 2
    wnd = get_window('hamming', int(400))
    Xframe = np.zeros((1024,num))
    for i in range(num) :
        x = xframe[:, i]
        x = x * wnd
        X = np.fft.fft(x, n=1024)
        Ax = np.abs(X)
        Xframe[:,i] = Ax

    numseg = t1.size
    ind1 = np.floor(fs1/fs*nfft).astype(int)
    ind2 = np.floor(fs2/fs*nfft).astype(int)
    f0 = np.zeros(t1.size)
    for j in  range(numseg):
        EA = Xframe[:, t1[j]:t2[j]+1 ].sum(axis =1)
        indfs = EA[ind1:ind2+1].argmax(0)
        f0[j] = (ind1+indfs+ 1)/nfft*fs
    return f0


def BandStop(fc ,fs):
    a = np.zeros(3)
    b = np.zeros(3)
    beta = 0.5
    fb = 600
    Q = beta * fc / fb
    K = math.tan (math.pi * fc / fs)

    b[0] = Q * (1 + K * K) / (K * K * Q + K + Q)
    b[1] = 2 * Q * (K * K - 1) / (K * K * Q + K + Q)
    b[2] = Q * (1 + K * K) / (K * K * Q + K + Q)
    a[0] = 1
    a[1] = 2 * Q * (K * K - 1) / (K * K * Q + K + Q)
    a[2] = (K * K * Q - K + Q) / (K * K * Q + K + Q)
    return  a ,b

def segfilter(x1 ,f0 , fs=16000 ,inc = 160 ,MM =2 ):
    Lin = x1.size
    a ,b = BandStop(f0 , fs)
    y1 = lfilter(b ,a , x1)
    y2 = y1[inc*2 : Lin - inc*2]
    Ly2 = y2.size
    w1 =  np.arange(inc*MM )/(inc*MM)
    w2 =  np.arange(inc*MM , 0 , -1) / (inc*MM)
    y3 = y2
    y3[0:inc*MM] = y3[0:inc*MM] * w1
    y3[ Ly2- inc*MM  :Ly2] =  y3[ Ly2- inc*MM  :Ly2] * w2

    return y3


if __name__ =="__main__" :
    # filename = '/Users/momo/Desktop/DSX/audio_16k/21.wav'
    filename = '../audio_16k/21.wav'
    outfilename = '../audio_out/0410/p_out_21.wav'
    x ,fs  = librosa.load(filename,sr =16000)
    xmax = x.max(0)
    x = x /xmax
    xframe = librosa.util.frame(x, frame_length=400, hop_length=160)
    xfft = librosa.stft(x, n_fft=1024, hop_length=160, win_length=400, window='hamming')
    cz , Am  = DetFeature(xframe,xfft)
    t1 ,t2 = DetTime(cz ,Am)

    a1 = [1.0000  , -0.0336  ,  0.3166]
    b1 = [0.5224  ,  0.4517  ,  0.1694]
    xlow = lfilter(b1 ,a1 ,x)
    f0 = DetFre(xlow , t1 ,t2)

    numseg =t1.size

    inc = 160
    wlen = 400
    NN = 4
    MM = 2
    w1 =  np.arange(inc*MM )/(inc*MM)
    w2 =  np.arange(inc*MM , 0 , -1) / (inc*MM)

    for j in range(numseg):
        t1[j]=t1[j]  + 1
        t2[j]=t2[j]  + 1
        indx1 = t1[j]*inc
        indx2 = t2[j]*inc + wlen

        sf1 = (t1[j] - NN)*inc
        sf2 = (t2[j] + NN)*inc +wlen
        segframe = xlow[sf1 : sf2]
        SegFrame = segfilter(segframe , f0[j] )

        xlowind1 = (t1[j] - MM)*inc
        xlowind2 = (t2[j] + MM)*inc +wlen
        xmid = x[xlowind1 :xlowind2]
        Weight =  np.concatenate((w2 , np.zeros(indx2-indx1) ,w1))
        Xmid = xmid* Weight
        Seg = Xmid + SegFrame
        x[ xlowind1 :xlowind2 ] =  Seg

    output = x*xmax

    write_wav(outfilename , output, 16000)



