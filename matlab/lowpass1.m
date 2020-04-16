clear ;close all;
% ¿¿Æ×¡¢
% filename ='/Users/momo/Desktop/out/0312/psbÐÂÈÎÎñ/old_male.wav';
filename ='./enhance/intest/25.wav';
[x  ,fs] = audioread(filename);

parametricEQ = fdesign.parameq('N,Flow,Fhigh,Gref,G0,GBW,Gst', 2, 4700, 8000 ,0,-15,-5,-1,fs);
Hd = design( parametricEQ,'cheby2');
[b1 ,a1] =tf(Hd);
% fvtool(Hd)

% a1=[ 1.0000   -0.3018]; 
% b1 = [ 0.3491    0.3491];
y1 = filter(b1,a1,x);
% y1 =y1*0.25;
% sound(y1,fs);
audiowrite('/Users/momo/Desktop/out/0324/0408/lowpass_2jie_out_25.wav',y1,fs);

