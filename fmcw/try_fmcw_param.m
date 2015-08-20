function try_fmcw_param()
    fc = 18e3;
    c = 300;
    lambda = c/fc;

    range_max = 1;
    % tm = 5*range2time(range_max,c);
    tm = 10e-3;

    range_res = 0.2;
    % bw = range2bw(range_res,c);
    bw = 2e3;
    sweep_slope = bw/tm;

    fr_max = range2beat(range_max,sweep_slope,c);

    v_max = 4;
    fd_max = speed2dop(2*v_max,lambda);

    fb_max = fr_max+fd_max;

    % fs = max(2*fb_max,bw);
    % fs = 48e3;
    fs = 2*bw;

    fprintf('freq: %.2fKHz\n', fc / 1e3);
    fprintf('max target range: %.2fm\n', range_max);
    fprintf('range resolution: %fm\n', range_res);
    fprintf('max target speed: %.2fm/s\n', v_max);
    fprintf('sweep time: %.2fms\n', tm*1e3);
    fprintf('sweep bandwidth: %.2fHz\n', bw);
    fprintf('max beat freq: %.2fHz\n', fb_max);
    fprintf('sample rate: %.2fHz\n', fs);

    hwav = phased.FMCWWaveform( ...
        'SweepTime', tm, ...
        'SweepBandwidth', bw,...
        'SampleRate', fs, ...
        'NumSweeps', 2);

    s = step(hwav);
    font_size = 20;
    % subplot(211); 
    % plot(0:1/fs:tm-1/fs,real(s));
    % xlabel('Time (s)', 'FontSize', font_size); 
    % ylabel('Amplitude (v)', 'FontSize', font_size);
    % set(gca, 'FontSize', font_size);
    % title('FMCW signal'); 
    % axis tight;

    % subplot(212); 
    figure(1);
    size(s)
    window = floor(fs/1024);
    noverlap = floor(window/4); % 75% overlap
    Nfft = fs;
    spectrogram(s, window, noverlap, Nfft, fs, 'yaxis');
    set(gca, 'FontSize', font_size);
    h = get(gca, 'xlabel');
    set(h, 'FontSize', font_size);
    h = get(gca, 'ylabel');
    set(h, 'FontSize', font_size);
    set(gca, 'YLim', [0 1*bw]);
    title('FMCW signal spectrogram');

    x = [zeros(10,1); s(1:end-10)];
    y = dechirp(x,s);

    figure(2);
    [Pxx,F] = periodogram(x,[],1024,fs,'centered');
    plot(F/1000,10*log10(Pxx)); grid;
    xlabel('Frequency (kHz)', 'FontSize', font_size);
    ylabel('Power/Frequency (dB/Hz)', 'FontSize', font_size);
    title('Periodogram Power Spectral Density Estimate Before Dechirping');

    figure(3);
    [Pxx,F] = periodogram(y,[],1024,fs,'centered');
    plot(F/1000,10*log10(Pxx)); grid;
    xlabel('Frequency (kHz)', 'FontSize', font_size);
    ylabel('Power/Frequency (dB/Hz)', 'FontSize', font_size);
    title('Periodogram Power Spectral Density Estimate Before Dechirping');


    hchannel = phased.FreeSpace('PropagationSpeed',c,...
        'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);


    
end
