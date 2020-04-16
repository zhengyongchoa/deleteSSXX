function  y = crosszeros(x)
wlen =length(x);
zcr1 =0;
for j=1: (wlen- 1)          % 在一帧内寻找过零点
         if x(j)* x(j+1)< 0       % 判断是否为过零点
             zcr1=zcr1+1;   % 是过零点，记录1次
         end
end

y = zcr1 ;
end