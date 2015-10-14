%nBits_ofdm = 48;
clear; clc;
load 'info.mat';
filename = 'audio_ofdm.wav';

fftlen = 64;
nSyms = 30;
cplen = fftlen;
modulator = modem.pskmod(2);
demodulator = modem.pskdemod(2);

mod_sym = modulate(modulator,preamble_bit);
tx_sym = ifft(mod_sym);

tx_sym = [tx_sym; tx_sym]; %Cyclic Prefix
tx_signals = repmat(tx_sym, nSyms,1);

% Upsample the transmitted signals (44.1 times)
%tx_upsampled = upsample(tx_signals, 44);
tx_upsampled = interp(tx_signals,44);
Fs = 44100;
ts = linspace(0, 50, Fs*50)';
ts = ts(1:length(tx_upsampled));
Fc = 16000;
%carrier_signal = cos(2*pi*Fc*ts);
carrier_signal = exp(1j*2*pi*Fc*ts);
carrier_sin = -sin(2*pi*Fc*ts);
carrier_cosin = cos(2*pi*Fc*ts);

tx_cos = real(tx_upsampled) .* carrier_cosin;
tx_sin = imag(tx_upsampled) .* carrier_sin;
tx_audio = tx_cos + tx_sin;

rx_cos_1 = tx_cos .* carrier_cosin;
rx_cos = tx_audio .* carrier_cosin;
[B, A] = butter(2, [0.0001 2000/Fs]);
rx_cos_1 = filter(B, A, rx_cos_1);
rx_cos = filter(B, A, rx_cos);

rx_sin_1 = tx_sin .* carrier_sin;
rx_sin = tx_audio .* carrier_sin;
rx_sin_1 = filter(B, A, rx_sin_1);
rx_sin = filter(B, A, rx_sin);

tx_passband = tx_upsampled .* carrier_signal;
%rx = tx_passband .* carrier_signal;

%rx = rx_cos + 1j*rx_sin;
rx = rx_cos_1 + 1j*rx_sin_1;

%rx = tx_upsampled;
rx_downsampled = downsample(rx,44);
rx_1 = fft(rx_downsampled(1:fftlen));
check = max(abs(rx_1 - mod_sym))
%scatterplot(rx_1);
%{
plot(abs(fftshift(fft(rx_cos))))
figure;
idx_len = 3000;
%plot([rx_cos_1(1:idx_len), rx_cos_2(1:idx_len), real(tx_upsampled(1:idx_len))])
plot([rx_sin_1(1:idx_len), rx_sin(1:idx_len), imag(tx_upsampled(1:idx_len))])
%}

% Add 1 second idle period before data transmission
idle = zeros(Fs,1);
tx_audio = [idle; tx_audio];
audiowrite(filename, tx_audio, Fs, 'BitsPerSample', 24);
