clear; clc;
filename = 'chirp_audio.wav';
Fs = 44100;
t = 10;
%interval = 1;
interval = 0.01;
total_time = 1000;
nRepeats = round(total_time / interval);
ts = linspace(0, interval, Fs*interval)';
%f0 = 0;
f0 = 16000;
f1 = 20000;
%f0 = 0;
%f1 = 10000;

%t = 0:1/1e3:2;
chirp_signal = chirp(ts,f0,interval,f1);
size(chirp_signal)
audio_signal = repmat(chirp_signal, nRepeats, 1);
audiowrite(filename, audio_signal, Fs, 'BitsPerSample', 24);
%sound(audio_signal);
fft_idx  = (-22050:1:22049)';
audio_fd = fftshift(abs(fft(audio_signal(1:44100))));
plot(fft_idx, audio_fd);
rx_fd = [fft_idx, audio_fd];

return;
Nfft = 2048*2;
window = floor(Nfft/4);
noverlap = floor(window/4); % 75% overlap
%[S,F,T,P] = spectrogram(audio_signal(1:4410), window, noverlap, Nfft, Fs, 'yaxis');
%spectrogram(chirp_signal, window, noverlap, Nfft, Fs, 'yaxis');
spectrogram(audio_signal, window, noverlap, Nfft, Fs, 'yaxis');
%imagesc(T, F, real(S));
%colorbar;
ylim([16500 20500]);
%ylim([0 11000]);

xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Time-Frequency plot of a Audio signal');
[S,F,T,P] = spectrogram(audio_signal(1:4410), window, noverlap, Nfft, Fs, 'yaxis');
