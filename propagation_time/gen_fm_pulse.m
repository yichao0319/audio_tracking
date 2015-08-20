%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%%
%% - Input:
%%
%%
%% - Output:
%%
%%
%% example:
%%
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gen_fm_pulse()
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
    

    %% --------------------
    %% Variable
    %% --------------------
    input_dir  = '';
    output_dir = '';


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, arg = 1; end
    if nargin < 1, arg = 1; end


    %% --------------------
    %% Main starts
    %% --------------------
    fs = 44100;                             % Sampling rate
    prop_speed = 331;                       % Propagation speed
    % max_range = 5;                          % 5m
    pulse_bw = fs / 2;                      % Pulse bandwidth
    pulse_width = 1/pulse_bw;               % Pulse width
    range_res = prop_speed / (2*pulse_bw);
    % prf = prop_speed/(2*max_range);         % Pulse repetition frequency
    prf = fs/1000;
    max_range = prop_speed / (2*prf);
    max_range


    fc = 18000;
    num_pulse_int = 10;
    lambda = prop_speed/fc;
    fast_time_grid = unigrid(0,1/fs,1/prf,'[)');
    slow_time_grid = (0:num_pulse_int-1)/prf;


    % hfm = phased.LinearFMWaveform('SampleRate', fs, ...
    %     'PulseWidth', pulse_width, ...
    %     'PRF', prf, ...
    %     'SweepBandwidth', pulse_bw, 'SweepDirection', 'Up',...
    %     'Envelope','Rectangular',...
    %     'OutputFormat', 'Pulses', 'NumPulses', 1);
    % hfm = phased.LinearFMWaveform('PulseWidth',100e-6,...
    %     'SweepBandwidth',2e5,'PRF',1e3);
    hfm = phased.LinearFMWaveform(...
        'SampleRate', fs, ...
        'PulseWidth', pulse_width, ...
        'PRF', prf, ...
        'SweepBandwidth', pulse_bw);
    hfm
    plot(hfm)

    y = step(hfm);
    size(y)

    %% up-convertion
    tx = real(y) .* cos(2*pi*fc/fs * [1:length(y)])';
    size(tx)

    %% receive
    % rx = [zeros(0.5*fs) tx zeros(1*fs) tx];
    rx = [zeros(0.5*fs, 1); tx];
    size(rx)

    %% down-convertion
    rx_base_real = rx .* cos(2*pi*fc/fs * [1:length(rx)])';
    rx_base_imag = rx .* sin(2*pi*fc/fs * [1:length(rx)])';

    [b,a] = butter(6, 1000/fs);
    freqz(b,a)
    rx_base_real = filter(b, a, rx_base_real);
    rx_base_imag = filter(b, a, rx_base_imag);
    rx_base = complex(rx_base_real, rx_base_imag);

    figure
    subplot(4,1,1)
    plot(real(y))
    title('tx base');

    subplot(4,1,2)
    plot(tx)
    title('tx up-convert');

    subplot(4,1,3)
    plot(rx)
    title('rx');

    subplot(4,1,4)
    plot(rx_base_real)
    title('rx base');



    matchingcoeff = getMatchedFilter(hfm);
    hmf = phased.MatchedFilter(...
        'Coefficients',matchingcoeff);
    [rx_pulses] = step(hmf, rx_base);

    matchingcoeff = getMatchedFilter(hfm);
    hmf = phased.MatchedFilter(...
        'Coefficients',matchingcoeff);
    [tx_pulses] = step(hmf, y);

    figure
    subplot(311);
    plot(real(rx));
    xlabel('Samples'); ylabel('Amplitude');
    title('Input Signal');

    subplot(312),plot(real(rx_pulses));
    xlabel('Samples'); ylabel('Amplitude');
    title('Matched Filter Output');

    % subplot(313),plot(real(rx_pulses_int));
    % xlabel('Samples'); ylabel('Amplitude');
    % title('Matched Filter Output');

    [~, idx] = max(tx_pulses)
    tx_time = idx / fs
    [~, idx] = max(rx_pulses)
    rx_time = idx / fs
    delay = rx_time - tx_time
    
    return
    
    % t = unigrid(0, 1/hfm.SampleRate, 1/hfm.PRF, '[)');
    % figure;
    % subplot(2,1,1)
    % plot(t,real(y))
    % axis tight;
    % title('Real Part');
    % subplot(2,1,2);
    % plot(t,imag(y)); xlabel('Seconds');
    % title('Imaginary Part');
    % axis tight;


    % x = [zeros(fs*0.5, 1); y];
    % x = y;
    % hlfm = hfm;
    % [afmag_lfm,delay_lfm,doppler_lfm] = ambgfun(x,...
    %     hlfm.SampleRate,hlfm.PRF);
    % [afmag_lfm,delay_lfm] = ambgfun(x,...
    %     hlfm.SampleRate, hlfm.PRF,...
    %     'Cut','Doppler');

    % figure;
    % imagesc(delay_lfm*1e6, doppler_lfm/1e3, afmag_lfm);
    % surf(delay_lfm*1e6,doppler_lfm/1e3,afmag_lfm,...
    %     'LineStyle','none');
    % axis tight; grid on; view([140,35]); colorbar;
    % xlabel('Delay \tau (\mus)');
    % ylabel('Doppler f_d (kHz)');
    % title('Linear FM Pulse Waveform Ambiguity Function');

    % figure
    % stem(delay_lfm, afmag_lfm)
    % xlabel('Delay (seconds)');
    % title('Autocorrelation of Linear FM Pulse');
    % axis([-2*range_res/prop_speed 2*range_res/prop_speed 0 1]); 
    % set(gca,'XTick',1e-5 * (-5:5))


    matchingcoeff = getMatchedFilter(hfm);
    hmf = phased.MatchedFilter(...
        'Coefficients',matchingcoeff);
    rx = [zeros(fs*0.5, 1); y];
    [rx_pulses] = step(hmf, rx);

    matchingcoeff = getMatchedFilter(hfm);
    hmf = phased.MatchedFilter(...
        'Coefficients',matchingcoeff);
    [tx_pulses] = step(hmf, y);

    rx_pulses_int = pulsint(real(rx_pulses));

    figure
    subplot(311);
    plot(real(rx));
    xlabel('Samples'); ylabel('Amplitude');
    title('Input Signal');

    subplot(312),plot(real(rx_pulses));
    xlabel('Samples'); ylabel('Amplitude');
    title('Matched Filter Output');

    subplot(313),plot(real(rx_pulses_int));
    xlabel('Samples'); ylabel('Amplitude');
    title('Matched Filter Output');

    % [~,range_detect] = findpeaks(real(rx_pulses));
    % range_detect / fs
    [~, idx] = max(tx_pulses)
    tx_time = idx / fs
    [~, idx] = max(rx_pulses)
    rx_time = idx / fs
    delay = rx_time - tx_time



    % matchingdelay = size(matchingcoeff,1)-1;
    % matchingdelay
    % size(x)
    % rx_pulses = buffer(rx_pulses(matchingdelay+1:end),size(rx_pulses,1));
    % figure
    % plot(rx_pulses)
    % % threshold = threshold * db2pow(mfgain);

    % range_gates = prop_speed*fast_time_grid/2;
    % htvg = phased.TimeVaryingGain(...
    %     'RangeLoss',2*fspl(range_gates,lambda),...
    %     'ReferenceLoss',2*fspl(max_range,lambda));

    % rx_pulses = step(htvg,rx_pulses);


    % rx_pulses = pulsint(rx_pulses,'noncoherent');
    % % helperRadarPulsePlot(rx_pulses,threshold,...
    % %     fast_time_grid,slow_time_grid,1);

    % threshold = 0.1
    % [~,range_detect] = findpeaks(rx_pulses);
    % % true_range = round(tgt_rng)
    % range_estimates = round(range_gates(range_detect))
end