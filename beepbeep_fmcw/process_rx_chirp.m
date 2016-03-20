%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%%
%% example:
%%   process_rx_chirp('rx1', 5000)
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
    if nargin < 1, filename = 'tx_chirp.17000.B2500.T0.10'; end
    if nargin < 2, fc = 17000; end

    tx_chirp_filename = sprintf('tx_chirp.%d.B2500.T0.10', fc);


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
    bw = 2500;
    sweep_slope = bw/tm;
    fr_max = range2beat(range_max,sweep_slope,c);

    fd_max = speed2dop(5,lambda);
    fb_max = fr_max + fd_max;

    % fs = max(2*fb_max,bw);
    fs = 44100;
    % fs = 44097.9;

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
    fprintf('  data size = %dx%d\n', size(analogData));

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
    analogData_orig = analogData;
    analogData = analogData .* cos(2*pi*fc*[1:length(analogData)]/fs).';


    if DEBUG2, fprintf('Low Pass Filter\n'); end
    analogData = lowPassFilterByFFT(analogData.', fs, 1.1*bw, 0).';
    % analogData = filter(B, A, analogData);


    %% ============================================================
    %% Roughly Sync Rx and Tx Chirt
    %% ============================================================
    if IF_TEST == 0
        if DEBUG2, fprintf('Roughly Sync Rx and Tx Chirt\n'); end

        for idx = 1:chirp_len
            rx_chirp = analogData(idx:idx+chirp_len-1);
            r = corrcoef(rx_chirp, tx_chirp);
            corr(idx) = r(1,2);
        end
        [~,idx] = max(corr);
        analogData = analogData(idx+chirp_len-200:end);
    end


    %% ====================
    %% Plot Chirp Signal
    %% ====================
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(2,2,[1,2]);
    spectrogram(analogData_orig(1:1*fs), 256, 64, 256, fs, 'yaxis');
    ylim([0 5]);
    title('RX Chirp'); axis tight;

    subplot(2,2,3);
    spectrogram(tx_chirp, 256, 64, 256, fs, 'yaxis');
    ylim([0 5]);
    title('TX Chirp');

    subplot(2,2,4);
    spectrogram(analogData(1:4*chirp_len), 256, 64, 256, fs, 'yaxis');
    ylim([0 5]);
    title('Sync RX Chirp'); axis tight;
    % return


    %% ============================================================
    %% DEBUG: simulate the shift
    %% ============================================================
    if IF_TEST
        if DEBUG2, fprintf('Simulate the Shift\n'); end

        gt_dist_dif  = 0.55; %% m
        shift_time   = gt_dist_dif / c;
        shift_sample = round(shift_time * fs);

        analogData = [zeros(shift_sample,1); analogData];
    end


    %% ============================================================
    %% Deal with clock skew
    %% ============================================================
    if DEBUG2, fprintf('Deal with clock skew\n'); end

    %% uncomment this section of codes to estimate true chirp length
    % gsize = 1;
    % num_rx_chirps = floor(length(analogData) / length(tx_chirp));
    % num_groups = floor(num_rx_chirps/gsize);
    % peaks = [];
    % Nc = 100;
    % idx_min = 1;
    % idx_max = chirp_len;
    % for gi = 1:Nc
    %     % gstd = (gi-1) * gsize + 1;
    %     % gend = gi * gsize;
    %     std_idx = (gi-1) * gsize * chirp_len + 1;
    %     end_idx = std_idx + gsize * chirp_len - 1;

    %     rx_chirp = analogData(std_idx:end_idx);

    %     % for idx = 1:chirp_len

    %     for idx = idx_min:idx_max
    %         tmp = analogData(std_idx + [idx:idx+gsize*chirp_len-1] - 1);
    %         r = corrcoef(tmp, repmat(tx_chirp, gsize, 1));
    %         corr(idx) = r(1,2);
    %     end
    %     [~,idx] = max(corr);
    %     if gi > 1
    %         if abs(idx - prev_idx) > 100
    %             continue;
    %         end
    %     end
    %     prev_idx = idx;
    %     sync_idx(gi) = std_idx + idx - 1;
    %     fprintf('  %d: sync idx = %d (%d)\n', gi, idx, sync_idx(gi));
    %     idx_min = max(1,idx - 100);
    %     idx_max = idx + 100;
    % end

    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % plot((sync_idx(2:end) - sync_idx(1:end-1))/gsize, '-b.'); hold on;
    % plot(((sync_idx(2:end) - sync_idx(1)) ./ (1:(Nc-1)))/gsize, '-r*')
    % grid('on');
    % itvl_first = ((sync_idx(2:end) - sync_idx(1)) ./ (1:(Nc-1))) / gsize;
    % itvl_prev  = (sync_idx(2:end) - sync_idx(1:end-1)) / gsize;
    % fprintf('  new chirp len: prev=%.2f, first = %.2f (last=%.2f)\n', mean(itvl_prev), mean(itvl_first), itvl_first(end));
    % new_chirp_len = mean(itvl_prev);
    % new_chirp_len = mean(itvl_first);
    % new_chirp_len = itvl_first(end);
    % new_chirp_len = mean(itvl_first(max(1,end-20):end));
    % new_chirp_len = 4407.94;
    new_chirp_len = 4407.95;
    fprintf('    select: %.2f\n', new_chirp_len);
    % return




    %% ============================================================
    %% Dechirp the Signal -- method 1
    %% ============================================================
    if DEBUG2, fprintf('Dechirp the Signal -- method 1\n'); end

    gsize = 1;
    num_rx_chirps = floor(length(analogData) / new_chirp_len);
    num_groups = floor(num_rx_chirps/gsize);
    for gi = 1:num_groups
        gstd = (gi-1) * gsize + 1;
        gend = gi * gsize;

        std_idx = round((gi-1) * gsize * new_chirp_len + 1);
        end_idx = std_idx + gsize * chirp_len - 1;

        rx_chirp = analogData(std_idx:end_idx);
        n = 2^nextpow2(length(rx_chirp));
        Y = fft(rx_chirp.*repmat(tx_chirp,gsize,1), n);

        P2 = abs(Y/n);
        P1 = P2(1:n/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        [~,idx] = max(P1(1:n/2));
        % freq = fs * (0:(n/2)) / n;
        freq = [0:(fs/n):(fs/2-fs/n)]';
        peak_f = freq(idx);
        range3(gi) = peak_f * 34600 * tm / bw;
        fprintf('  %d: range=%f, idx=%d\n', gi, range3(gi), idx);

        if gi == 1
            %% used to plot figure (fix the max freq to plot)
            sel_idx = 2*idx;
        end


        fh = figure(1); clf;
        subplot(2,1,1)
        times = [1:length(range3)] * tm;
        plot(times, range3, '-b.');
        hold on;
        plot(times(gi), range3(gi), 'ro')
        grid('on')

        subplot(2,1,2)
        spectrogram(rx_chirp,256,64,256,fs,'yaxis');
        set(gca, 'ylim', [0 5])
        % waitforbuttonpress;
    end


    % fix_range = 20;
    % slope = (range3(fix_range) - range3(1)) / fix_range;
    % range3 = range3 - slope * [0:length(range3)-1];
    range3 = range3 - range3(1);


    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    plot(range3, '-b.');
    grid('on')
    ylabel('range (cm)')
    % set(gca, 'ylim', [0 100])
    return




    %% ============================================================
    %% Dechirp the Signal -- method 1 -- plot movie
    %% ============================================================
    if DEBUG2, fprintf('Dechirp the Signal -- method 1\n'); end

    gsize = 1;
    num_rx_chirps = floor(length(analogData) / new_chirp_len);
    num_groups = floor(num_rx_chirps/gsize);
    for gi = 1:num_groups
        gstd = (gi-1) * gsize + 1;
        gend = gi * gsize;

        std_idx = round((gi-1) * gsize * new_chirp_len + 1);
        end_idx = std_idx + gsize * chirp_len - 1;

        rx_chirp = analogData(std_idx:end_idx);
        n = 2^nextpow2(length(rx_chirp));
        Y = fft(rx_chirp.*repmat(tx_chirp,gsize,1), n);

        P2 = abs(Y/n);
        P1 = P2(1:n/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        [~,idx] = max(P1(1:n/2));
        freq = [0:(fs/n):(fs/2-fs/n)]';
        peak_f = freq(idx);


        fh = figure(5); clf;
        subplot(2,1,1)
        times = [1:length(range3)] * tm;
        plot(times, range3, '-b.');
        hold on;
        plot(times(gi), range3(gi), 'ro')
        grid('on')

        subplot(2,1,2)
        plot(freq, P1(1:n/2), '-b.')
        set(gca, 'xlim', [0 freq(sel_idx)])

        % print(fh, '-dpng', sprintf('./movie/%s.%.2d.png', filename, gi));
        % print(fh, '-dpsc', sprintf('./movie/%s.%.2d.eps', filename, gi));
        % waitforbuttonpress

    end

    % fix_range = 20;
    % slope = (range3(fix_range) - range3(1)) / fix_range;
    % range3 = range3 - slope * [0:length(range3)-1];
    range3 = range3 - range3(1);


    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    plot(range3, '-b.');
    grid('on')
    % set(gca, 'ylim', [0 100])
    return


end

