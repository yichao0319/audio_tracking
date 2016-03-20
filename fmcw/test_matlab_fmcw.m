close all;
clear all;

fc = 17000;
c = 300;
lambda = c/fc;
range_max = 10;
% tm = 5.5*range2time(range_max,c);
tm = 0.1;
range_res = 0.01;
% bw = range2bw(range_res, c);
bw = 3000;
sweep_slope = bw/tm;
fr_max = range2beat(range_max,sweep_slope,c);

fd_max = speed2dop(5,lambda);
fb_max = fr_max + fd_max;

% fs = max(2*fb_max,bw);
fs = 44100;


hw = phased.FMCWWaveform('SweepBandwidth', bw,...
   'SampleRate',fs,'SweepDirection','Up','SweepTime',tm, ...
   'NumSweeps',1);
x_base = step(hw);
% x = x .* exp(1i*2*pi*fc*(1:(tm*fs))/fs)';
x = x_base .* cos(2*pi*fc*[1:length(x_base)]/fs)';
x = lowPassFilterByFFT(x', fs, 22000, 0)';

% analogData = real(x);
analogData = x;
analogData = analogData .* cos(2*pi*fc*[1:length(analogData)]/fs)';
analogData = lowPassFilterByFFT(analogData', fs, 2*bw, 0)';


fh = figure(1); clf;
% spectrogram(x,32,16,128,hw.SampleRate,'yaxis');
% size(x)
subplot(211); plot(real(x));
xlabel('Time (s)'); ylabel('Amplitude (v)');
title('FMCW signal'); axis tight;
subplot(212); spectrogram(x,32,16,32,fs,'yaxis');
title('FMCW signal spectrogram');
% waitforbuttonpress


ranges = [1:100:3000];
for ri = 1:length(ranges)
    N = round(ranges(ri) / 30000 * fs);
    fprintf('range=%fcm, N=%d, len-x=%d', ranges(ri), N, length(analogData));
    xr = [zeros(N,1); analogData(1:end-N,1)];
    % xr = [x(end-N+1:end,1); x(1:end-N,1)];
    xd = dechirp(xr,x_base);

    fb_rng = rootmusic(pulsint(xd,'coherent'), 1, fs);
    rng_est(ri) = beat2range(fb_rng, sweep_slope, c)*2;
    fprintf('  est range=%fcm\n', rng_est(ri)*100);
end


return

fb=17000;
Fs=44100;
Ts=1/Fs;
interval=0.04;
K=Fs*interval;
round=1000;
T=round*interval;
N=Fs*T;
L=1000;

signal=zeros(1,N);
signal_left=zeros(1,N);
signal_right=zeros(1,N);


% generate frequency modulated signals
fminL=14000;
fminR=10500;
B=2500;
w=0.5-0.5*cos(2*pi/K*(0:K-1)); % hann window

for i=1:round
    t=(0:K-1)*Ts;
    signal_left((1:K)+(i-1)*K)=cos(2*pi*(1/2*t.^2*B/interval+fminL*t));
    signal_right((1:K)+(i-1)*K)=cos(2*pi*(1/2*t.^2*B/interval+fminR*t));
end

signal_d=[transpose(signal_left), transpose(signal_right)];

% spectrogram(signal_left(1:K), 32, 16, 128, Fs,'yaxis');
% waitforbuttonpress

signal=signal_left+signal_right;
plot(abs(signal))
ylim([-2,2])

Fsig=fft(signal);
plot(abs(Fsig));
N=numel(Fsig)/2;
Fsig2=[Fsig(N+1:2*N) Fsig(1:N)];
plot((-N:N-1)/N*Fs/2,abs(Fsig2));

audiowrite('./toneChirpSoundLR.wav',signal_d,Fs,'BitsPerSample', 16);
% audiowrite('./toneMixedSin.wav',signal,Fs,'BitsPerSample', 16);