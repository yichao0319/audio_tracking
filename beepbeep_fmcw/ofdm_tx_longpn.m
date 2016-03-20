
%% ofdm_tx_longpn(128, 1)
function ofdm_tx_longpn(pn_len, config)
    % clear all;

    if nargin < 1, pn_len = 128; end
    if nargin < 2, config = 1; end


    output_dir = './tx_sound/';

    % paramters
    % config=1;
    if config==1
        Nfft=128;
        Ncp=32;
        symbolRate=1/2205;
        Ts=1/44100;
        Fs=1/Ts;
        fc1=12000;
        fc2=5000;
        rollOff=0.25;
        nSamp=20;
        filterSpan=10;
        % preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
        preambleSeq=m_sequence([1 0 0 0 0 1 1 1 0 1 0 1 0 0 0 1 0]);
    elseif config==2
        Nfft=128;
        Ncp=32;
        symbolRate=1/2205;
        Ts=1/44100;
        Fs=1/Ts;
        fc1=5000;
        fc2=12000;
        rollOff=0.25;
        nSamp=20;
        filterSpan=10;
        % preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
        preambleSeq=m_sequence([1 0 0 0 0 1 1 1 0 1 0 1 0 0 0 1 0]);
    elseif config==3
        Nfft = 128;
        Ncp = 32;
        symbolRate = 1/2205;
        Ts = 1/44100;
        Fs = 1/Ts;
        fc1 = 5000;
        fc2 = 12000;
        rollOff = 0.25;
        nSamp = 20;
        filterSpan = 10;
        % preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
        preambleSeq=m_sequence([1 0 0 0 0 1 1 1 0 1 0 1 0 0 0 1 0]);
    else
        error('wrong config');
    end



    % generate data
    data=preambleSeq;
    % data=data(1:Nfft);
    data = data(1:min(end,pn_len));
    data = reshape(data(1:Nfft*floor(length(data)/Nfft)), Nfft, []);
    % modulate with BPSK
    data=data*2-1;

    % generate OFDM symbol
    data_fft = ifft(data,Nfft);
    data_fft = reshape(data_fft, 1, []);
    preamble_size = length(data_fft);
    fprintf(' Size of Preamble: %d (~%fs)\n', preamble_size, preamble_size*Nfft*Ts);
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generating multiple symbols
    % intervalSignal=zeros(1, ceil(Fs*2));
    % fprintf('  one packet size with interval = %dx%d\n', size([intervalSignal, analogData, intervalSignal]));
    % Ns=60;
    % oneSymbol = [intervalSignal, analogData, intervalSignal];
    % multiSymbol=[];
    % for i=1:Ns
    %     multiSymbol=[multiSymbol, intervalSignal, analogData, intervalSignal];
    %     % fprintf('  one packet size with interval = %dx%d\n', size([intervalSignal, analogData, intervalSignal]));
    % end
    % analogData=multiSymbol;
    %%%%%%%%
    Ns = 30;
    oneSymbol = zeros(1, Fs*1);
    oneSymbol(1:length(analogData)) = analogData;
    analogData = repmat(oneSymbol, 1, Ns);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % upconverting
    T=numel(analogData);
    analogData=real(analogData.*exp(1i*2*pi*fc1*(1:T)*Ts));
    analog_one_symbol = real(oneSymbol.*exp(1i*2*pi*fc1*(1:numel(oneSymbol))*Ts));
    % dataSpectrum=fft(analogData);
    % dataSpectrum=circshift(dataSpectrum,[0 floor(T/2)]);
    % Th=(T-1)/2;
    % plot((-Th:Th)/Th*Fs/2,abs(dataSpectrum));

    % scaling
    analogData=analogData/max(abs(analogData))*0.9;
    analog_one_symbol = analog_one_symbol / max(abs(analog_one_symbol)) * 0.9;
    plot(abs(analogData));



    % % model channel effects
    % htap=5;
    % h=zeros(1,htap*nSamp);
    % h(1:nSamp:(htap-1)*nSamp+1)=[1 5 2 1 1];
    % analogData=conv(analogData,h);
    % plot(abs(analogData));

    % write the results
    save([output_dir 'preamble.' num2str(config) '.' num2str(fc1) '.' num2str(preamble_size) '.mat'],'preamble');
    save([output_dir 'txch_signal.' num2str(config) '.' num2str(fc1) '.' num2str(preamble_size) '.mat'],'analogData');
    save([output_dir 'txch_1signal.' num2str(config) '.' num2str(fc1) '.' num2str(preamble_size) '.mat'],'analog_one_symbol');



    %% Chirp Signal
    time_len = length(analogData) / Fs;
    [chirp_analogData] = gen_chirp(fc2, config, time_len, output_dir);


    % generate audio file
    signal_d = zeros(T,2);
    signal_d(:,2)= analogData.' + chirp_analogData;
    signal_d(:,2) = signal_d(:,2) / max(abs(signal_d(:,2))) * 0.9;
    audiowrite([output_dir 'ofdm.' num2str(config) '.' num2str(fc1) '.' num2str(preamble_size) '.wav'], signal_d, Fs, 'BitsPerSample', 16);


    fh = figure(1); clf;
    spectrogram(signal_d(1+10*Fs:13*Fs,2), 256, 64, 256, Fs, 'yaxis');
    ylim([5 15]);
end






function [chirp_analogData] = gen_chirp(fc, config, time_len, output_dir)

    %% --------------------
    %% DEBUG
    %% --------------------
    DEBUG0 = 0;
    DEBUG1 = 1;
    DEBUG2 = 1;  %% progress
    DEBUG3 = 1;  %% verbose
    DEBUG4 = 1;  %% results


    %% --------------------
    %% Constant
    %% --------------------
    % output_dir = './tx_sound/';


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, fc = 12000; end
    if nargin < 2, time_len = 60; end


    %% --------------------
    %% Variable
    %% --------------------
    fig_idx = 0;

    fs = 44100;
    c = 300;

    if config == 1
        tm = 0.1;
        bw = 2500;
    elseif config == 2
        tm = 0.1;
        bw = 2500;
    elseif config == 3
        tm = 0.04;
        bw = 2500;
    else
        error('unknown config');
    end



    %% --------------------
    %% Main starts
    %% --------------------

    times = [0:1/fs:tm-1/fs]';
    %% --------------------
    %% Generate Chirp -- Method 1
    %% --------------------
    % hw = phased.FMCWWaveform('SweepBandwidth', bw,...
    %    'SampleRate',fs,'SweepDirection','Up','SweepTime',tm, ...
    %    'NumSweeps',1);
    % x_base = step(hw);
    % x = x_base .* exp(1i*2*pi*fc*(1:(tm*fs))/fs)';

    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % subplot(211); spectrogram(x_base,256,64,256,fs,'yaxis');
    % set(gca, 'ylim', [0, 5])
    % title('Base signal spectrogram');

    % subplot(212); spectrogram(x,256,64,256,fs,'yaxis');
    % title('FMCW signal spectrogram');
    % return


    %% --------------------
    %% Generate Chirp -- Method 2
    %% --------------------
    hw = phased.FMCWWaveform('SweepBandwidth', bw,...
       'SampleRate',fs,'SweepDirection','Up','SweepTime',tm, ...
       'NumSweeps',1);
    x_base = step(hw);
    x = x_base .* (cos(2*pi*fc*[1:length(x_base)]/fs).');
    x = lowPassFilterByFFT(x.', fs, 20000, 0).';

    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % subplot(211); spectrogram(x_base,256,64,256,fs,'yaxis');
    % set(gca, 'ylim', [0, 5])
    % title('Base signal spectrogram');

    % subplot(212); spectrogram(x,256,64,256,fs,'yaxis');
    % % set(gca, 'ylim', [10, 16])
    % title('FMCW signal spectrogram');


    %% --------------------
    %% Generate Chirp -- Method 3
    %% --------------------
    % x = chirp(times,fc,tm,fc+bw);

    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % subplot(211); plot(times,real(x));
    % xlabel('Time (s)'); ylabel('Amplitude (v)');
    % title('FMCW signal'); axis tight;

    % subplot(212); spectrogram(x,256,64,256,fs,'yaxis');
    % % set(gca, 'ylim', [0, 5000])
    % % set(gca, 'ylim', [0, 5000]+fc-1000)
    % title('FMCW signal spectrogram');


    num_chirp = floor(time_len / tm);
    time_len = num_chirp * tm;
    fprintf('  #chirps=%d\n', num_chirp);

    %% scaling
    chirp_signal = x/max(abs(x))*0.9;
    chirp_signal_base = x_base/max(abs(x_base))*0.9;

    %% repeating
    analogData = repmat(chirp_signal, num_chirp, 1);


    %% save file
    filename = sprintf('tx_chirp.%d.B%d.T%.2f', fc, bw, tm);
    save([output_dir filename '.mat'], 'chirp_signal_base');

    % signal_d = zeros(time_len*fs, 2);
    % signal_d(:,2) = analogData;
    % audiowrite([output_dir filename '.wav'], signal_d, fs, 'BitsPerSample', 16);

    chirp_analogData = analogData;
end

