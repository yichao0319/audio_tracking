clear all;

output_dir = './gen_data/';

% paramters
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
elseif config==4
    Nfft=16;
    Ncp=8;
    symbolRate=1/2940;
    Ts=1/44100;
    Fs=1/Ts;
    fc=18250;
    rollOff=0.25;
    nSamp=15;
    filterSpan=10;
    preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
end



% generate data
data=preambleSeq;
data=data(1:Nfft);

% modulate with BPSK
data=data*2-1;

% generate OFDM symbol
data_fft=ifft(data,Nfft);
% plot(abs(data_fft));
ofdmSymbol=[data_fft(Nfft-Ncp+1:Nfft), data_fft];
plot(abs(ofdmSymbol));

% pulse shaping
firConfig=fdesign.pulseshaping(nSamp,'Raised Cosine','Nsym,Beta',filterSpan,rollOff);
fir=design(firConfig);
firValue=fir.Numerator;
preamble=upfirdn(ofdmSymbol, firValue, nSamp);
analogData=[zeros(1,Fs/100) preamble zeros(1,Fs/100)];
plot(abs(preamble));
% disp(numel(analogData));
% return

fprintf('  one packet size = %dx%d\n', size(analogData));

% generating multiple symbols
intervalSignal=zeros(1,Fs/1);
Ns=60;
multiSymbol=[];
for i=1:Ns
    multiSymbol=[multiSymbol, intervalSignal, analogData, intervalSignal]; 
    fprintf('  one packet size with interval = %dx%d\n', size([intervalSignal, analogData, intervalSignal]));
end
analogData=multiSymbol;

% upconverting
T=numel(analogData);
analogData=real(analogData.*exp(1i*2*pi*fc*(1:T)*Ts));
% dataSpectrum=fft(analogData);
% dataSpectrum=circshift(dataSpectrum,[0 floor(T/2)]);
% Th=(T-1)/2;
% plot((-Th:Th)/Th*Fs/2,abs(dataSpectrum));

% scaling
analogData=analogData/max(abs(analogData))*0.9;
plot(abs(analogData));

% % model channel effects
% htap=5;
% h=zeros(1,htap*nSamp);
% h(1:nSamp:(htap-1)*nSamp+1)=[1 5 2 1 1];
% analogData=conv(analogData,h);
% plot(abs(analogData));

% write the results
if config==1
    save([output_dir 'preambleC1.mat'],'preamble');
elseif config==2
    save([output_dir 'preambleC2.mat'],'preamble');
elseif config==3
    save([output_dir 'preambleC3.mat'],'preamble');
elseif config==4
    save([output_dir 'preambleC4.mat'],'preamble');
end
save([output_dir 'txch_signal.' num2str(fc) '.mat'],'analogData');


% generate audio file
signal_d=zeros(T,2);
signal_d(:,2)=analogData.';
audiowrite([output_dir 'ofdm.' num2str(fc) '.wav'], signal_d, Fs, 'BitsPerSample', 16);
