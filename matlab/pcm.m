clear;

[X, Fs]= audioread('../case_good/g1.wav');

x = X * 32768;

fileId = fopen( '../case_good/g1.pcm','w');

fwrite(fileId, x,'int16') ;

fclose(fileId);



f1= fopen('/Users/momo/Desktop/չʾ/songbai.pcm','r');

y= fread(f1 ,inf,'int16');

fclose(f1);

Y =y / 32768 ;



plot(Y);

figure

plot(X);

