clear all;

% parameters
config=1;
if config==1
    Nfft=128;
    Ncp=32;
    symbolRate=1/2205;
    Ts=1/44100;
    Fs=1/Ts;
    fc=18000;
    rollOff=0.25;
    nSamp=20;
    filterSpan=10;
    preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
    sampleInterval=0.16;
elseif config==2
    Nfft=256;
    Ncp=128;
    symbolRate=1/4410;
    Ts=1/44100;
    Fs=1/Ts;
    fc=18000;
    rollOff=0.25;
    nSamp=10;
    filterSpan=10;
    preambleSeq=m_sequence([0 0 0 0 1 0 0 0 1]);
elseif config==3
    Nfft=64;
    Ncp=16;
    symbolRate=1/2205;
    Ts=1/44100;
    Fs=1/Ts;
    fc=18000;
    rollOff=0.25;
    nSamp=20;
    filterSpan=10;
    preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
    sampleInterval=0.08;
end
sampleDelay=filterSpan*nSamp/2+nSamp*Ncp;

% load signals
% load 'txch_signal.mat';
fileIndex='1436306148';
fileSuffix='_point_0_audio_audio.wav';
fileName=strcat(fileIndex,fileSuffix);
[analogData,~]=audioread(fileName);
analogData=analogData.';

% downconvert
T=numel(analogData);
analogData=analogData.*exp(-1i*2*pi*fc*(1:T)*Ts);
% analogDataFFT=fft(analogData);
% plot(abs(analogDataFFT));
analogData=lowPassFilterByFFT(analogData,Fs,2000,500);
analogDataFFT=fft(analogData);
plot(abs(analogDataFFT));

% get the signals for the preamble
if config==1
    load 'preambleC1.mat';
elseif config==2
    load 'preambleC2.mat';
elseif config==3
    load 'preambleC3.mat';
end

% course timing synchornizing
windowSize=300;
detectLength=sampleInterval/Ts*2;
courseStartIndex=findStartIndexByDoubleWin(analogData,windowSize,detectLength);

% split the received signals
dataindex=courseStartIndex;
dataSpan=(Nfft+Ncp-1)*nSamp+nSamp*filterSpan+1+Fs/100*2;
dataInterval=sampleInterval/Ts-dataSpan;
Ns=400;
analogData=split(analogData,dataindex,dataSpan,dataInterval,Ns);

Nd=numel(analogData(1,:));
Np=numel(preamble);
h=zeros(Ns,Nfft);
startIndex=zeros(1,Ns);
startIndexCorr=zeros(1,Ns);
for j=1:Ns
    
    % refined timing sychronization  
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
    
    startIndex(j)=maxCorrIndex+sampleDelay;
    startIndexCorr(j)=0;
    minh=100;
    data_orgn=2*preambleSeq-1;
    data_orgn=data_orgn(1:Nfft);
    for i=startIndex(j)-10:startIndex(j)+10
        % downsampling
        data_fft=analogData(j,i:nSamp:i+nSamp*(Nfft-1));
        
        % convert to frequency domain
        data=fft(data_fft);
        
        % calculate channel gain        
        Hf=data./data_orgn;
        htemp=ifft(Hf);
        if abs(htemp(Nfft))<minh
            minh=abs(htemp(Nfft));
            ht=htemp;
            startIndexCorr(j)=i;
        end
    end
    disp(startIndexCorr(j)-sampleDelay);
    % stem(abs(ht));
    h(j,:)=ht;
end
h1Enrg=abs(h(:,1)).^2.';

% get distance
fileSuffix='_point_0_audio_truth.txt';
gtFileName=strcat(fileIndex,fileSuffix);
fileID=fopen(gtFileName);
dataMat=textscan(fileID,'%d64 %d %d');
fclose(fileID);
gtTime=dataMat{1};
gtTime=gtTime-gtTime(1);
pixelX=dataMat{3};
pixelY=dataMat{2};

% translate pixels to positions
initPosition=[0.98;0.58];
initPixel=[112;330];
pixelScale=1/970;
posiX=(double(pixelX)-initPixel(1))*pixelScale+initPosition(1);
posiY=(double(pixelY)-initPixel(2))*pixelScale+initPosition(2);
p=[posiX;posiY];

% calculate the recording time
recordingTime=((courseStartIndex+startIndex)*Ts+(0:Ns-1)*0.08)*1000;

% calculate the distance
currentTimeIndex=1;
recordingPosi=zeros(1,Ns);
dist=zeros(1,Ns);
C1=[0;0];
CHeight=0.085;
for i=1:Ns
    currentTimeIndex=findNearestGroundTruthTime(recordingTime(i),gtTime,currentTimeIndex);
    recordingPosi(i)=p(currentTimeIndex);
    dist(i)=normCorr(CHeight,recordingPosi(i)-C1);
end
res=[dist;h1Enrg];

 scatter(log(h1Enrg),log(dist));
 ylabel('distance (log scale)');
 xlabel('Estimated Direct Path Energy (log scale)');
 