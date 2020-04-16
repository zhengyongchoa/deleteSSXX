function Y = segfilter(x , inc ,fs ,FS)
winL1 =2 ;
winL2 =2 ;

Lx = length(x);

fs1 = fs-300;
fs2 = fs+300;
parametricEQ = fdesign.parameq('N,Flow,Fhigh,Gref,G0,GBW,Gst', 2, fs1, fs2 ,0,-40,-10,-1,FS);
Hd = design( parametricEQ,'cheby2');
[b ,a] =tf(Hd);

y1 =filter(b,a, x);
y2 = y1(inc*winL1+1: Lx-winL1*inc);
Ly2 =length(y2);
% figure(9)
% plot(x);
% hold on
% plot(y1 ,'r');

w1 = (0:inc*winL2 -1)'./(inc*winL2) ;
w2 = (inc*winL2 : -1 :1)'./(inc*winL2);

y3 =y2;
y3(1:inc*winL2) = y2(1:inc*winL2).* w1 ;
y3(Ly2 -inc*winL2+1: Ly2) = y2( Ly2 -inc*winL2+1: Ly2) .* w2;
% figure(10)
% plot(y2,'k');
% hold on
% plot(y3 ,'r');


Y = y3;
end

