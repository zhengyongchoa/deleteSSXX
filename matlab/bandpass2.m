clear ;close all
filename ='./enhance/intest/24.wav';
[x  ,fs] = audioread(filename);
% fs = 16000 ;
fc = 3000 ;
fb = 600 ;

beta = 0.2; % 0.2 ; 0.5 ; 主要决定陷波的深度
Q =beta* fc / fb ;
K = tan(pi* fc/fs);

b0 =   Q*(1+K^2)/(K^2*Q + K + Q);
b1 =   2*Q*(K^2 - 1)/(K^2*Q + K + Q) ;
b2 =   Q*(1+ K^2) / (K^2*Q + K + Q) ;
a0 =   1 ;
a1 =   2*Q*(K^2 - 1)/(K^2*Q + K + Q) ;
a2 =  (K^2*Q - K + Q) /(K^2*Q + K + Q) ;

a=[a0 a1  a2 ];
b=[b0 b1  b2 ];
[H,f]=freqz(b,a,1000,fs);   % 求数字滤波器幅频曲线
y1=filter(b,a,x);           % 对降采样后的数据进行滤波
figure(1)
plot(f,20*log10(abs(H)),'k');
grid; 
xlabel('频率/Hz'); ylabel('幅值')
title('巴特沃斯滤波器的幅值响应')

audiowrite('/Users/momo/Desktop/out/0324/0407/0.2_out24.wav',y1,fs);




