function y=lowPassFilterByFFT(y,Fs,bw,guard)

Y=fft(y);
N=numel(Y);
Nb=floor(bw/Fs*N);
Ng=floor(guard/Fs*N);

filter=zeros(1,N);
filter(1:Nb+1)=1;
filter(N-Nb+1:N)=1;
filter(Nb+1:Nb+Ng)=cos((0:Ng-1)/Ng*pi/2);
filter(N-Nb-Ng+2:N-Nb+1)=cos((Ng:-1:1)/Ng*pi/2);
plot(filter);
Y=Y.*filter;
y=ifft(Y);
