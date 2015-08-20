clear; clc;
filename = 'chirp_record.wav';
Fs = 44100;
t = 10;
interval = 0.01;
%interval = 1;
fft_size = Fs*interval*10;
fft_width = Fs / fft_size;
%fft_idx  = linspace(-Fs/2, fft_width, Fs/2)';
fft_idx  = (-22040:10:22050)';
%fft_idx  = (-22049:1:22050)';
%fft_idx  = (-22050:1:22049)';

ts = linspace(0, interval, Fs*interval)';
f0 = 16000;
f1 = 20000;
%f0 = 0;
%f1 = 10000;
audio_rx = audioread(filename);
audio_rx = audio_rx(5000:end);
rx_tmp = audio_rx(1:fft_size);
fft_out = fft(rx_tmp);
plot(fft_idx, fftshift(abs(fft_out)));
rx_fd = [fft_idx, abs(fft_out)];



%{
Nfft = 2048*2;
window = floor(Nfft/4);
noverlap = floor(window/4); % 75% overlap
%[S,F,T,P] = spectrogram(audio_signal(1:4410), window, noverlap, Nfft, Fs, 'yaxis');
%spectrogram(chirp_signal, window, noverlap, Nfft, Fs, 'yaxis');
spectrogram(audio_rx, window, noverlap, Nfft, Fs, 'yaxis');
%imagesc(T, F, real(S));
%colorbar;
%ylim([17000 21000]);
ylim([f0 f1]);

xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Time-Frequency plot of a Audio signal');
%[S,F,T,P] = spectrogram(audio_signal(1:4410), window, noverlap, Nfft, Fs, 'yaxis');
%}