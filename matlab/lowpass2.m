clear  ; clc; close all;
filename ='./enhance/intest/23.wav';
[x  ,fs] = audioread(filename);
% fs =16000;
fs2=fs/2;                       % �����������Ƶ�ʵ�һ��
fp1= 2500;                     % ͨ��Ƶ��
fs1= 5000;                     % ���Ƶ��
wp1=fp1/fs2;                    % ��һ��ͨ��Ƶ��
ws1=fs1/fs2;                    % ��һ�����Ƶ��
Ap=3; As=30;                    % ͨ�����ƺ����˥��
[n,Wn]=buttord(wp1,ws1,Ap,As);  % ���˲���ԭ�ͽ����ʹ���
n =1;
[b,a]=butter(n,Wn);         % �������˲���ϵ��
% a=[ 1.0000   -0.3018]; 
% b = [ 0.3491    0.3491];
[H,f]=freqz(b,a,1000,fs);   % �������˲�����Ƶ����

y1=filter(b,a,x);           % �Խ�����������ݽ����˲�
% figure(1)
% plot(f,(abs(H)))
figure(2)
plot(f,20*log10(abs(H)),'k');
grid; 
xlabel('Ƶ��/Hz'); ylabel('��ֵ')
title('������˹�˲����ķ�ֵ��Ӧ')
set(gcf,'color','w');

audiowrite('./out/0325/lowpass1_1_23out.wav',y1,fs);

