clear ;close all;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
% log： 0401 c版本
% log:  0409 c版本
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
filename = '../case_good/g3.wav';
outfile ='../audio_out/out_good/out_m_g2.wav';
% filename = '../audio_16k/female.wav';
% outfile ='../audio_out/0410/out_m_female.wav';
[x  ,fs] = audioread(filename);
xmax=max(x) ;
x = x./max(x); 


%% load lowpass filter
% a1=[ 1.0000   -0.3018]; 
% b1 = [ 0.3491    0.3491];
% a1=[ 1.0000   -0.5252    0.2208];
% b1 = [ 0.1739    0.3478    0.1739];
a1=[1.0000   -0.0336    0.3166]; %  比较弱的低通滤波器
b1=[0.5224    0.4517    0.1694];
xlow = filter(b1, a1 , x);

%%
L1 = length(x);
time=(0:L1-1)/fs;
inc =160;
nfft = 1024;
wlen =400;
num = floor((L1-wlen)/inc )-1;
frameTime=frame2time(num, wlen, inc, fs);% 计算各帧对应的时间
sf =zeros(num,1);
maxcz =150;
Afs1 = 1 ;   Afs2 = 1800;
Afs3 = 3000; Afs4 = fs/2 ;
numfs1 = ceil(Afs1/fs*nfft) ;
numfs2 = floor(Afs2/fs*nfft) ;
numfs3 = floor(Afs3/fs*nfft);
numfs4 =  floor(Afs4/fs*nfft);
for i = 1:num
   ind1 =  inc* (i-1) ;
   ind2  = ind1 +wlen ;
   y = x(ind1+1 : ind2);
   ylow = xlow(ind1+1 : ind2);
    wnd=hamming(wlen);                        % 设置窗函数
    y = y.* wnd;
    ylow = ylow.* wnd;
    A=abs(fft(y,nfft));                      % 取来一帧数据FFT后取幅值
    Alow = abs(fft(ylow,nfft)); 

    Am(i) = sum(A(numfs3:512)) - sum(A(1:numfs2)) ;
    if Am(i) >0   
        Am(i) =1;
    else     
        Am(i) =0 ;     
    end
    
    cz(i) =  crosszeros(y);
    
    if  cz(i)>maxcz  &&  Am(i)==1
        sf(i) = 1;
    end
%     AA(:,i) = A;
    AAlow(:,i) = Alow ;
end

% figure(2)
% plot(Am);
% title('子带能量');
% figure(5)
% plot(cz)
% title('过0率');
% figure(6)
% plot( frameTime ,sf);
% title('简单求时间节点');


%%                     
minlen  = 9;    % 初始化
status  = 0;
count   = 0;
segn= 1;

for i = 1:num
    
    switch status
        case {0,1}                      %  0 =  静音， 1 = 开始
            if Am(i)  &&  cz(i)>maxcz
                 x1(segn) = i;
                 status = 2;
                 count = count + 1;
                
            elseif Am(i)  ||  cz(i)>maxcz         % 可能处于语音段
                 status = 1;
                 count = count + 1;
            else
                 status = 0;                    % 静音段
                 count = 0;
                 x1(segn) =0 ;
                 x2(segn) = 0;
            end


        case 2                              % 语音段
            if Am(i)  &&  cz(i)>maxcz
                 count = count + 1;
            else
                if count < minlen            % 噪音段    
                     status = 0;
                     count = 0;  
                else               
                    status =3; 
                    x2(segn) = x1(segn) + count;
                end 
           end

        case 3  
             segn = segn+1;
             status = 0;
             count = 0;
             x1(segn) = 0;
             x2(segn) = 0;
    
    end
end


numseg = length(x1);
if x1(1)==0 
   disp('没有检测到SSXX ！');
   return;
end
if x1(numseg)==0, numseg = numseg -1; end
if x2(numseg) ==0
    x2(numseg) = num;
end

decinc =ceil( wlen/inc);
for i= 1:numseg
    x2(i) = x2(i) - decinc ;
    
end

SF=zeros(1,num);                         % 按x1和x2，对SF和NF赋值
EAind1 =  floor(2000/fs*1024);
EAind2 = floor(7500/fs*1024);
for i=1 : numseg
    SF(x1(i):x2(i))=1;
%     EA = sum(AA(:,x1(i):x2(i)),2);
%     tmp1 = AAlow(:,x1(i):x2(i)) ;
    EA = sum(AAlow(:,x1(i):x2(i)),2);
    [~ , indi] = max( EA(EAind1 :EAind2));
    x3(i) =  (indi+ EAind1-1) / 1024 *fs;
end
figure(16)
plot(frameTime,SF)
% plot(SF)
title('第2种方法')

% 
NN =4 ; MM=2;
w1 = (0:inc*MM -1)'./(inc*MM);
w2 = (inc*MM : -1 :1)'./(inc*MM);

for i=1 : numseg

    indx1 = x1(i)*inc +1;
    indx2 = x2(i)*inc + wlen;
%     indx2 = x2(i)*inc-inc ;
    % xlow进行滤波
    indNN1 = (x1(i) -NN)*inc +1;
    indNN2 = (x2(i) +NN)*inc +wlen;
    segframe = xlow(indNN1 : indNN2);
    tmpXX = x(indNN1 : indNN2);
    Segframe = segfilter(segframe , inc ,x3(i) ,fs);
    testSeg = Segframe;testSeg(inc*2+1:end-inc*2) = 0 ;
    %xlow 头尾操作
    indMM1 = (x1(i) - MM)*inc +1;
    indMM2 = (x2(i) + MM)*inc +wlen;
%     xmid = xlow( indMM1: indMM2);
    xmid = x( indMM1: indMM2);
    Weight = [w2 ;zeros((indx2-indx1+1) ,1) ;w1];
    Xmid = xmid .* Weight;
%     figure(1) ;plot(x(indx1 : indx2));title('only_x');

     Seg =  Xmid+Segframe;
     x(indMM1: indMM2) = Seg;

%     figure(2) ;plot(Segframe);title('Segframe');
%     figure(3);plot(Xmid,'r');title('xlow中部归零 ');hold on;
%     plot(testSeg);
    
%     figure(4);plot(Seg);title('合并');
%     figure(5);plot(x(segframe_inx1: segframe_inx2) ,'r' );title('合并后的4+N+4');
%     hold on ;plot(tmpXX,'b');
end
output =x.*xmax;
% figure(1)
% plot(x);
audiowrite(outfile,output,fs);












