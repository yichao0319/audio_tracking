clear all;

% parameters
Nfft=128;
Ncp=32;
symbolRate=1/2205;
Ts=1/44100;
Fs=1/Ts;
fc=18000;
rollOff=0.25;
nSamp=20;
filterSpan=10;
sampleDelay=filterSpan*nSamp/2+nSamp*Ncp;

% load signals
load 'tx_signal.mat';

% downconvert
T=numel(analogData);
analogData=analogData.*exp(-1i*2*pi*fc*(1:T)*Ts);
% analogDataFFT=fft(analogData);
% plot(abs(analogDataFFT));
analogData=lowPassFilterByFFT(analogData,Fs,2000,500);
analogDataFFT=fft(analogData);
plot(abs(analogDataFFT));

% get the signals for the preamble
load 'preamble.mat';

% split the received signals
dataindex=Fs/20+1;
dataSpan=(Nfft+Ncp-1)*nSamp+nSamp*filterSpan+1+Fs/100*2;
dataInterval=Fs/10;
Ns=100;
analogData=split(analogData,dataindex,dataSpan,dataInterval,Ns);

Nd=numel(analogData(1,:));
Np=numel(preamble);
h=zeros(Ns,Nfft);
for j=1:Ns
    
    % timing sychronization    
    corr=zeros(1,Nd-Np+1);
    maxCorr=0;
    maxCorrIndex=0;
    for i=0:Nd-Np
        corr(i+1)=abs(analogData(j,i+1:i+Np)*preamble');
        if corr(i+1)>maxCorr
            maxCorr=corr(i+1);
            maxCorrIndex=i+1;
        end
    end
    disp(maxCorrIndex);
    % plot(corr);
    startIndex=maxCorrIndex+sampleDelay;
    
    % downsampling
    data_fft=analogData(j,startIndex:nSamp:startIndex+nSamp*(Nfft-1));
    
    % convert to frequency domain
    data=ifft(data_fft);
    
    % calculate channel gain
    data_orgn=2*m_sequence([1 0 0 0 0 1 1 1])-1;
    data_orgn=data_orgn(1:Nfft);
    Hf=data./data_orgn;
    ht=ifft(Hf);
    % stem(abs(ht));
    h(j,:)=ht;
end
