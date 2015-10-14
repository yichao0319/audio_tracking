clear all;

config=1;
if config==1
    load 'preambleC1.mat';
elseif config==2
    load 'preambleC2.mat';
elseif config==3
    load 'preambleC3.mat';
end

% generating 80ms/160ms symbols
Fs=44100;
Ts=1/Fs;
K=0.16/Ts;

symbol=zeros(1,K);
Np=numel(preamble);
symbol(1:Np)=preamble;

Ns=5000;
analogData=zeros(1,Ns*K);
for i=1:Ns
    analogData((1:K)+(i-1)*K)=symbol;
end

% upconverting
fc=18000;
T=numel(analogData);
analogData=real(analogData.*exp(1i*2*pi*fc*(1:T)*Ts));
plot(analogData);

% scaling
analogData=analogData/max(abs(analogData))*0.98;
plot(analogData);

% generate audio file
signal_d=zeros(T,2);
signal_d(:,2)=analogData.';
audiowrite('./ofdmN128CP32BW2kSeq.wav',signal_d,Fs,'BitsPerSample',16);
    