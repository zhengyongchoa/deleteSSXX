function  y = crosszeros(x)
wlen =length(x);
zcr1 =0;
for j=1: (wlen- 1)          % ��һ֡��Ѱ�ҹ����
         if x(j)* x(j+1)< 0       % �ж��Ƿ�Ϊ�����
             zcr1=zcr1+1;   % �ǹ���㣬��¼1��
         end
end

y = zcr1 ;
end