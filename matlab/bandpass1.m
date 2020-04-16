clear ;close all;
% filename ='/Users/momo/Desktop/out/0312/psbÐÂÈÎÎñ/old_male.wav';
filename ='./enhance/intest/25.wav';
[x  ,fs] = audioread(filename);

%
parametricEQ = fdesign.parameq('N,Flow,Fhigh,Gref,G0,GBW,Gst', 2, 1700, 2300 ,0,-40,-10,-1,fs);
Hd = design( parametricEQ,'cheby2');
[b1 ,a1] =tf(Hd);
fvtool(Hd)


y1 = filter(b1,a1,x);
% y1 =y1*0.25;
% sound(y1,fs);
audiowrite('./out/0325afternoon/bandpassN_2_25out.wav',y1,fs);


