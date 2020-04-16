clear  ; clc; close all;
filename ='./enhance/intest/23.wav';
[x  ,fs] = audioread(filename);
% fs =16000;
fs2=fs/2;                       % 降采样后采样频率的一半
fp1= 2500;                     % 通带频率
fs1= 5000;                     % 阻带频率
wp1=fp1/fs2;                    % 归一化通带频率
ws1=fs1/fs2;                    % 归一化阻带频率
Ap=3; As=30;                    % 通带波纹和阻带衰减
[n,Wn]=buttord(wp1,ws1,Ap,As);  % 求滤波器原型阶数和带宽
n =1;
[b,a]=butter(n,Wn);         % 求数字滤波器系数
% a=[ 1.0000   -0.3018]; 
% b = [ 0.3491    0.3491];
[H,f]=freqz(b,a,1000,fs);   % 求数字滤波器幅频曲线

y1=filter(b,a,x);           % 对降采样后的数据进行滤波
% figure(1)
% plot(f,(abs(H)))
figure(2)
plot(f,20*log10(abs(H)),'k');
grid; 
xlabel('频率/Hz'); ylabel('幅值')
title('巴特沃斯滤波器的幅值响应')
set(gcf,'color','w');

audiowrite('./out/0325/lowpass1_1_23out.wav',y1,fs);

