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
%%   process_rx_chirp()
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function process_rx_chirp(filename, fc)

    %% --------------------
    %% DEBUG
    %% --------------------
    DEBUG0 = 0;
    DEBUG1 = 1;
    DEBUG2 = 1;  %% progress
    DEBUG3 = 1;  %% verbose
    DEBUG4 = 1;  %% results

    IF_TEST = 0;


    %% --------------------
    %% Constant
    %% --------------------
    input_dir  = './rx_sound/';
    tx_dir     = './tx_sound/';
    output_dir = '';


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, filename = 'tx_chirp.15000.B3000.T0.10'; end
    if nargin < 2, fc = 15000; end

    tx_chirp_filename = sprintf('tx_chirp.%d.B3000.T0.10', fc);


    %% --------------------
    %% Variable
    %% --------------------
    fig_idx = 0;

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

    [B,A] = butter(4, bw/fs);


    %% --------------------
    %% Main starts
    %% --------------------

    %% ============================================================
    %% Load data
    %% ============================================================
    if DEBUG2, fprintf('Load data\n'); end

    %% TX Chirp
    tmp = load([tx_dir tx_chirp_filename '.mat']);
    tx_chirp = tmp.chirp_signal_base;
    chirp_len = length(tx_chirp);
    fprintf('  chirp size = %dx%d\n', size(tx_chirp));

    %% RX Audio
    audio_filename1 = [input_dir filename '.wav'];
    audio_filename2 = [input_dir filename '.aac'];
    audio_filename3 = [tx_dir filename '.wav'];
    if exist(audio_filename1, 'file') == 2,
        audio_filename = audio_filename1;
    elseif exist(audio_filename2, 'file') == 2,
        audio_filename = audio_filename2;
    elseif exist(audio_filename3, 'file') == 2,
        audio_filename = audio_filename3;
        IF_TEST = 1;
    else
        error('Cannot find the file');
    end
    fprintf('  audio file: %s\n', audio_filename);
    [analogData,~] = audioread(audio_filename);

    %% ----------------
    %% XXX: need to modify
    % len = length(analogData) / Fs
    if size(analogData,2) > 1
        analogData = analogData(:,2);
    end
    %% ----------------
    fprintf('  data size = %dx%d\n', size(analogData));


    %% ====================
    %% downconvert
    %% ====================
    if DEBUG2, fprintf('Downconvert\n'); end

    T = numel(analogData);
    analogData = analogData .* exp(-1i*2*pi*fc*(1:T)/fs)';

    if DEBUG2, fprintf('Low Pass Filter\n'); end
    % analogData = lowPassFilterByFFT(analogData', fs, 2000, 500)';
    % analogData = lowPassFilterByFFT(analogData', fs, 2000, 500)';
    analogData = filter(B, A, analogData);


    %% ====================
    %% Plot Tx Chirp Signal
    %% ====================
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(211);
    plot(real(tx_chirp));
    xlabel('Time (s)'); ylabel('Amplitude (v)');
    title('FMCW signal'); axis tight;

    subplot(212);
    spectrogram(tx_chirp,256,64,256,fs,'yaxis');
    set(gca, 'YLim', [0 5]);
    title('TX Chirp');
    % return


    %% ============================================================
    %% Roughly Sync Rx and Tx Chirt
    %% ============================================================
    if DEBUG2, fprintf('Roughly Sync Rx and Tx Chirt\n'); end

    if IF_TEST == 0
        for idx = 1:chirp_len
            rx_chirp = analogData(idx:idx+chirp_len-1);
            r = corrcoef(rx_chirp, tx_chirp);
            corr(idx) = r(1,2);
        end
        [~,idx] = max(corr);
        analogData = analogData(idx+100:end);
    end


    %% ====================
    %% Plot Rx Chirp Signal
    %% ====================
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(211);
    plot(real(analogData(1:chirp_len)));
    xlabel('Time (s)'); ylabel('Amplitude (v)');
    title('FMCW signal'); axis tight;

    subplot(212);
    spectrogram(analogData(1:chirp_len),256,64,256,fs,'yaxis');
    set(gca, 'YLim', [0 5]);
    title('RX signal');



    %% ============================================================
    %% DEBUG: simulate the shift
    %% ============================================================
    if DEBUG2, fprintf('Simulate the Shift\n'); end

    if IF_TEST
        gt_dist_dif  = 1.4; %% m
        shift_time   = gt_dist_dif / c;
        shift_sample = round(shift_time * fs);

        analogData = [zeros(shift_sample,1); analogData];
    end


    %% ============================================================
    %% Dechirp the Signal
    %% ============================================================
    if DEBUG2, fprintf('Dechirp the Signal\n'); end

    num_rx_chirps = floor(length(analogData) / length(tx_chirp));
    for ci = 1:num_rx_chirps-10
        std_idx = (ci-1) * length(tx_chirp) + 1;
        end_idx = ci * length(tx_chirp);
        rx_chirp = analogData(std_idx:end_idx);

        xd = dechirp(rx_chirp, tx_chirp);
        xr(:,ci) = xd;

        fb_rng = rootmusic(pulsint(xd,'coherent'),1,fs);
        range1(ci) = beat2range(fb_rng, sweep_slope, c)*2;


        %% ------------
        %% my dechirp1
        %% ------------
        ir = real(ifft(fft(rx_chirp)./fft(tx_chirp)));
        ir = abs(ir);
        ir(1) = 0;
        ir(end) = 0;
        ir = ir(1:round(fs*tm*3/5));
        ir = filter(B,A,ir);

        [~,idx] = max(ir);
        range2(ci) = c * idx / fs;

        fh = figure(4); clf;
        plot(ir, '-b.');
        set(gca, 'xlim', [0, 2000])

        fprintf('  %d/%d: peak freq = %f\n', ci, num_rx_chirps, idx)


        %% ------------
        %% my dechirp2
        %% ------------
        ir = fft(rx_chirp.*tx_chirp);
        % ir = filter(B,A,ir);
        ir = abs(ir);

        [~,idx] = max(ir);
        range3(ci) = idx/chirp_len * tm * c * 2;
        % range3(ci) =

        fh = figure(5); clf;
        plot(ir, '-b.');
        set(gca, 'xlim', [0, 2000])


        fprintf('  %d/%d: peak freq = %f\n', ci, num_rx_chirps, idx)
        % fprintf('  real=%.2f, est=%.2f, my_est=%.2f\n', gt_dist_dif, range(ci), R);
        fprintf('    est1=%.2f, est2=%.2f, est3=%.2f\n', range1(ci), range2(ci), range3(ci));

        % return
        waitforbuttonpress
    end

    % range1 = range1 - range1(1);
    % range2 = range2 - range2(1);

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    plot(range1, '-r')
    hold on;
    plot(range2, '-b')

    fb_rng = rootmusic(pulsint(xr,'coherent'),1,fs);
    rng_est = beat2range(fb_rng,sweep_slope,c) * 2;

    % fprintf('  real=%.2f, est=%.2f\n', gt_dist_dif, rng_est)
    fprintf('  est=%.2f\n', rng_est)


    % hrdresp = phased.RangeDopplerResponse('PropagationSpeed',c,...
    %     'DopplerOutput','Speed', 'OperatingFrequency',fc,...
    %     'SampleRate',fs, 'RangeMethod','FFT', 'SweepSlope',sweep_slope,...
    %     'RangeFFTLengthSource','Property','RangeFFTLength',2048,...
    %     'DopplerFFTLengthSource','Property','DopplerFFTLength',256);


end

