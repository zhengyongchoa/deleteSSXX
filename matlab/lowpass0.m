% function [y, N] = LowPassFilter( x, Fs, Fc)
clear ; clc; close all;
 filename ='/Users/momo/Desktop/DSX/enhance/intest/21.wav';
 [x  ,Fs] = audioread(filename);

Fc =5500;

m=exp(-1/(Fs/Fc));
b = [1-m , 0];
a = [1 , -m] ;

[H ,f] = freqz(b ,a ,1000 , Fs );
figure(1)
plot(f , 20*log10(abs(H)));



y=filter(b,a,x);

audiowrite('./out/0326/lowpass_fc_21out.wav',y,Fs);
