clear ;close all
filename ='./enhance/intest/24.wav';
[x  ,fs] = audioread(filename);
% fs = 16000 ;
fc = 3000 ;
fb = 600 ;

beta = 0.2; % 0.2 ; 0.5 ; ��Ҫ�����ݲ������
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
[H,f]=freqz(b,a,1000,fs);   % �������˲�����Ƶ����
y1=filter(b,a,x);           % �Խ�����������ݽ����˲�
figure(1)
plot(f,20*log10(abs(H)),'k');
grid; 
xlabel('Ƶ��/Hz'); ylabel('��ֵ')
title('������˹�˲����ķ�ֵ��Ӧ')

audiowrite('/Users/momo/Desktop/out/0324/0407/0.2_out24.wav',y1,fs);




