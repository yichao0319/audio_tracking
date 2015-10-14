%nBits_ofdm = 48;
function dec_audio_ofdm
    clear; clc;
    load 'info.mat';
    filename = 'RecordAudio.wav';    
    %filename = 'audio_ofdm.wav';
    
    rx_audio = audioread(filename);
    rx_audio = mean(rx_audio,2);
    %fftlen = 64;
    nSyms = 100;
    demodulator = modem.pskdemod(2);
    Fs = 44100;
    Fc = 16000;
    data_idx = [3:28, 45:55];
    
    %tx_sym = [tx_sym; tx_sym]; %Cyclic Prefix
    %tx_signals = repmat(tx_sym, nSyms,1);

    rx_signal = pass2base(rx_audio, Fc, Fs);
    idx = synchronization(rx_signal, fftlen, preamble_td)
    rx_preamble = fft(rx_signal(idx:idx+fftlen-1));
    %scatterplot(rx_sym);
    h = rx_preamble.*conj(preamble_fd);

    for x=1:20
        close all;
        rx_sym = fft(rx_signal(idx + fftlen*2*x:fftlen*2*x+idx+fftlen-1));
        rx_sym = rx_sym ./ h;
        rx_sym = rx_sym(data_idx);
        %evm = abs(rx_sym-preamble_fd);
        %plot(evm)
        scatterplot(rx_sym);        
        %axis([-2, 2, -2, 2]);
    end
        
    %{
    for x=-5:5
        rx_1 = fft(rx_signal(x+idx:x+idx+fftlen-1));
        scatterplot(rx_1);
    end
    %}
    return;
end

function idx = synchronization(rx_signal, fftlen, preamble_td)
    corr_result = xcorr(rx_signal, preamble_td);
    corr_result = corr_result(length(rx_signal):end);
    corr_result = abs(corr_result);
    %plot(abs(corr_result),'r');
    %hold on;
    %plot(abs(rx_signal),'g');    
    
    
    %% Synchronization
    %1. Initial finding
    start_idx = 4000;
    thr = 0.001;
    while(corr_result(start_idx) < thr)
        start_idx = start_idx+1;
    end
    start_idx
    
    %2. Find the peak
    %Use the cross-correlation results of 20 consecutive symbols
    nSyms = 10;
    corr_total = zeros(fftlen,1);
    corr_result = corr_result(start_idx: start_idx+ fftlen*nSyms-1);
    %plot(corr_result)
    corr_result = reshape(corr_result,fftlen, nSyms);
    mean_corr = mean(corr_result,2);
    [~, peak_idx] = max(mean_corr);
    idx = start_idx + peak_idx-1;
end

%pass_signal is audio signal that are real numbers
function base_signal = pass2base(pass_signal, Fc, Fs)
    rx_len = length(pass_signal);
    ts = linspace(0, 50, Fs*50)';
    ts = ts(1:rx_len);
    carrier_sin = -sin(2*pi*Fc*ts);
    carrier_cosin = cos(2*pi*Fc*ts);
    base_real = pass_signal .* carrier_cosin;
    base_imag = pass_signal .* carrier_sin;

    [B, A] = butter(2, [0.0001 2000/Fs]);
    base_real = filter(B, A, base_real);
    base_imag = filter(B, A, base_imag);
    base_signal = base_real + 1j*base_imag;
    base_signal = downsample(base_signal,44);
end