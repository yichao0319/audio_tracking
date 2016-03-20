%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%%
%% example:
%%   retrieve_signals('rx.3.1', 3)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function retrieve_signals(filename, config)
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
    rx_dir     = '../beepbeep_fmcw/rx_sound/';
    tx_dir     = '../beepbeep_fmcw/tx_sound/';
    output_dir = './spikes/';
    fig_dir    = './fig/';

    c = 34600;


    %% --------------------
    %% Variable
    %% --------------------
    fig_idx = 0;

    %% PN Seq
    max_pn_len        = 128;
    new_itvl_samp     = 44079.5;

    [Nfft, Ncp, symbolRate, Ts, fs, fc_pn, fc_fmcw, rollOff, nSamp, filterSpan, preambleSeq] = get_config_param(config);
    pn_len            = Nfft*floor(max_pn_len/Nfft);
    preamble_filename = [tx_dir 'preamble.' num2str(config) '.' num2str(fc_pn) '.' num2str(pn_len) '.mat'];
    tx_filename       = [tx_dir 'txch_1signal.' num2str(config) '.' num2str(fc_pn) '.' num2str(pn_len) '.mat'];

    %% FMCW
    [bw, tm, new_chirp_len] = get_fmcw_config_param(config);
    tx_chirp_filename = sprintf('tx_chirp.%d.B%d.T%.2f', fc_fmcw, bw, tm);


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

    %% TX PN seq One Signal
    fprintf('  tx signal file: %s\n', tx_filename);
    tmp = load(tx_filename);
    tx_analogData = tmp.analog_one_symbol;
    itvl_samp = length(tx_analogData);
    fprintf('  tx signal size = %dx%d\n', size(tx_analogData));

    %% TX Preambles
    fprintf('  preamble file: %s\n', preamble_filename);
    load(preamble_filename);
    fprintf('  preamble size = %dx%d\n', size(preamble));


    %% RX Audio
    audio_filename1 = [rx_dir filename '.wav'];
    audio_filename2 = [rx_dir filename '.aac'];
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
    analogData = analogData.';
    fprintf('  data size = %dx%d\n', size(analogData));

    %% ----------------
    %% XXX: need to modify
    if size(analogData,1) > 1
        analogData = analogData(2,:);
    end
    %% ----------------
    fprintf('  data size = %dx%d\n', size(analogData));

    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % spectrogram(analogData(1+0*fs:3*fs), 256, 64, 256, fs, 'yaxis');
    % ylim([5 15]);
    % return



    %% ====================
    %% downconvert
    %% ====================
    if DEBUG2, fprintf('Downconvert\n'); end

    T = numel(analogData);

    %% PN
    analogData_pn = analogData .* exp(-1i*2*pi*fc_pn*[1:T]/fs);
    analogData_pn = lowPassFilterByFFT(analogData_pn, fs, 1500, 0);

    %% FMCW
    analogData_fmcw = analogData .* exp(-1i*2*pi*fc_fmcw*[1:T]/fs);
    % analogData_fmcw = analogData .* cos(2*pi*fc_fmcw*[1:T]/fs);
    analogData_fmcw = lowPassFilterByFFT(analogData_fmcw, fs, 1.1*bw, 0);



    %% ======================================================
    %% PN: course timing synchornizing
    %% ======================================================
    if DEBUG2, fprintf('PN: course timing synchornizing\n'); end

    windowSize = 3000;
    detectLength = ceil(itvl_samp * 1.5);
    [courseStartIndex, mag_corr] = findStartIndexByDoubleWin(analogData_pn, windowSize, detectLength);
    fprintf('  Course Start Index = %d (%fs)\n', courseStartIndex, courseStartIndex/fs);

    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % subplot(3,1,1);
    % plot([1:T]/fs, abs(analogData_pn), '-b');
    % hold on;
    % plot(ones(1,2)*courseStartIndex/fs, get(gca, 'ylim'), '-r');
    % set(gca, 'xlim', [1 detectLength]/fs);
    % subplot(3,1,2);
    % plot([1:length(mag_corr)+windowSize]/fs, [zeros(1, windowSize), mag_corr]);
    % hold on;
    % plot(ones(1,2)*courseStartIndex/fs, get(gca, 'ylim'), '-r');
    % % set(gca, 'xlim', [1 detectLength]);
    % subplot(3,1,3);
    % spectrogram(analogData, 256, 64, 256, fs, 'yaxis');
    % hold on;
    % plot(ones(1,2)*courseStartIndex/fs, [0,15]);
    % ylim([0,15]);
    % xlim([0,2]);


    %% ======================================================
    %% PN: calculate corrcoeff
    %% ======================================================
    if DEBUG2, fprintf('PN: calculate corrcoeff\n'); end

    Np = numel(preamble);
    prev_idx = courseStartIndex;
    i = courseStartIndex;
    corr = zeros(1, T);
    norm_corr = corr;
    while i <= T
        show_progress(i-courseStartIndex+1, T-courseStartIndex+1, 1);

        if i+Np-1 > T
            break;
        end

        corr(i) = abs(analogData_pn(i:i+Np-1)*(preamble.'));

        if (i - prev_idx > (itvl_samp/4))

            %% used to find the offset between PN and FMCW
            % if prev_idx == courseStartIndex
            %     idx_offset = sync_pn_by_fmcw(analogData, analogData_pn, analogData_fmcw, corr, courseStartIndex, windowSize, tx_chirp, new_chirp_len);
            %     % idx_offset = 449.05;
            %     waitforbuttonpress
            % end


            idx = [prev_idx:i];
            norm_corr(idx) = corr(idx) / max(corr(idx));
            i = prev_idx + itvl_samp - 1;
            prev_idx = i;

            continue;
        end

        i = i + 1;
    end


    %% ============================================================
    %% FMCW: Roughly Sync Rx and Tx Chirt
    %% ============================================================
    if DEBUG2, fprintf('FMCW: Roughly Sync Rx and Tx Chirt\n'); end

    for idx = 1:chirp_len
        tmp = analogData_fmcw(idx:idx+chirp_len-1);
        r = corrcoef(tmp, tx_chirp);
        corr_fmcw(idx) = r(1,2);
    end
    [~,fmcw_sync_idx] = max(corr_fmcw);
    analogData_fmcw = analogData_fmcw(fmcw_sync_idx+chirp_len-200:end);


    %% ============================================================
    %% FMCW: Dechirp the Signal
    %% ============================================================
    if DEBUG2, fprintf('FMCW: Dechirp the Signal\n'); end

    gsize = 10; %% must be 1 ..
    num_rx_chirps = floor(length(analogData_fmcw) / new_chirp_len);
    num_groups = floor(num_rx_chirps/gsize);
    for gi = 1:num_groups
        std_idx = round((gi-1) * gsize * new_chirp_len + 1);
        end_idx = std_idx + gsize * chirp_len - 1;

        rx_chirp_time{gi} = std_idx:end_idx;
        rx_chirp{gi} = analogData_fmcw(std_idx:end_idx);
        n = 2^nextpow2(length(rx_chirp{gi}));
        Y = fft(rx_chirp{gi}.*(repmat(tx_chirp,gsize,1).'), n);

        P2 = abs(Y/n);
        P1 = P2(1:n/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        fmcw_fft{gi} = P1(1:n/2);

        [~, fmcw_peak_idx(gi)] = max(fmcw_fft{gi});
        freq = [0:(fs/n):(fs/2-fs/n)]';
        % peak_f = freq(idx);
        % range3(gi) = peak_f * 34600 * tm / bw;
        % fprintf('  %d: range=%f, idx=%d\n', gi, range3(gi), idx);

        % if gi == 1
        %     %% used to plot figure (fix the max freq to plot)
        %     sel_idx = 2*idx;
        % end


        % fh = figure(1); clf;
        % subplot(2,1,1)
        % times = [1:length(range3)] * tm;
        % plot(times, range3, '-b.');
        % hold on;
        % plot(times(gi), range3(gi), 'ro')
        % grid('on')

        % subplot(2,1,2)
        % spectrogram(rx_chirp,256,64,256,fs,'yaxis');
        % set(gca, 'ylim', [0 5])
        % waitforbuttonpress;
    end

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    plot(fmcw_peak_idx, '-b.');
    title('FMCW peak index');
    % waitforbuttonpress

    % if config <= 2
    %     % gt_fmcw_peak_idx = mean(fmcw_peak_idx);
    %     gt_fmcw_peak_idx = mean(fmcw_peak_idx(1:5));
    %     if abs(gt_fmcw_peak_idx-23) < abs(gt_fmcw_peak_idx-24)
    %         gt_fmcw_peak_idx = 23;
    %     else
    %         gt_fmcw_peak_idx = 24;
    %     end
    % elseif config == 3
    %     gt_fmcw_peak_idx = 14;
    % else
    %     error('  wrong config');
    % end
    gt_fmcw_peak_idx = median(fmcw_peak_idx(1:5));
    gt_fmcw_peak_dist = freq(gt_fmcw_peak_idx) * c * tm / bw;
    fmcw_dist = freq * c * tm / bw;


    %% plot distance
    dist_range = [-200 200];
    pn_dist_range = round(dist_range * fs / c);
    fmcw_dist_range = round(dist_range * bw / tm / c * n / fs);
    dist_unit = min(c/fs, fs/n*c*tm/bw);
    dist_ranges = [max(pn_dist_range(1)*fs/c, fmcw_dist_range(1)*fs/n*c*tm/bw):dist_unit:min(pn_dist_range(2)*fs/c, fmcw_dist_range(2)*fs/n*c*tm/bw)];
    fprintf('  range = %.2f~%.2fcm, unit=%.2fcm\n', dist_ranges(1), dist_ranges(end), dist_unit);
    fprintf('  pn range = %d~%d, fmcw range = %d~%d\n', pn_dist_range, fmcw_dist_range);



    %% ======================================================
    %% PN: Find peak
    %% ======================================================
    if DEBUG2, fprintf('PN: Find peak\n'); end

    idx = courseStartIndex;
    cnt = 1;
    %% corrcoeff
    peak_idx = [];
    %% normalized corrcoeff
    norm_peak_idx = [];
    norm_first_peak_idx = [];
    %% ground truth
    gt_peak_idx = [];
    %% thresholds
    norm_thresh = 0.95;
    status = [];

    while(idx + itvl_samp < length(corr))
        seg_idx = idx:idx+itvl_samp-1;
        seg = corr(seg_idx);

        %% max peak
        % [v, this_peak] = max(seg);
        [v, this_peak] = max(corr(idx:round(idx+windowSize+Np*1/2)));
        % [v, this_peak] = max(corr(idx:round(idx+windowSize+10)));
        peak_idx = [peak_idx idx+this_peak-1];

        %% max peak of normalized corr
        norm_seg = norm_corr(seg_idx);
        % [v, this_peak] = max(norm_seg);
        [v, this_peak] = max(norm_corr(idx:round(idx+windowSize+Np*1/2)));
        % [v, this_peak] = max(norm_corr(idx:round(idx+windowSize+10)));
        norm_peak_idx = [norm_peak_idx idx+this_peak-1];

        %% first peak of normalized corr > thresh
        [v, locs] = findpeaks(norm_seg);
        tmp = locs(find(v>norm_thresh));
        if length(tmp) > 0
            norm_first_peak_idx = [norm_first_peak_idx idx+tmp(1)-1];
            this_peak = tmp(1);
        end

        %% ground truth peak
        if cnt == 1
            tmp_first_idx = get_gt_first_peak_idx(filename, analogData_pn, corr, courseStartIndex, windowSize);
            if tmp_first_idx < 0
                tmp_first_idx = norm_peak_idx(1);
            end

            num_segs = floor(T / new_itvl_samp);
            gt_peak_idx = tmp_first_idx*ones(1, num_segs) + round([0:num_segs-1]*(new_itvl_samp));
            % gt_peak_idx = norm_first_peak_idx(1)*ones(1, num_segs) + round([0:num_segs-1]*(new_itvl_samp));
        end


        %% ======================================================
        %% FMCW: Find close chirps
        %% ======================================================
        chirp_cnt = floor(gt_peak_idx(cnt) / new_itvl_samp * (new_itvl_samp/new_chirp_len/gsize)) + 1;
        if chirp_cnt > length(fmcw_fft)
            break;
        end
        this_fmcw_fft = fmcw_fft{chirp_cnt};
        fprintf('  %d PN ~ %d FMCW\n', cnt, chirp_cnt);


        %% -------------------
        fh = figure(11); clf;
        subplot(6,2,[1:4]);
        spectrogram(analogData(seg_idx), 256, 64, 256, fs, 'yaxis');
        hold on;
        plot(ones(1,2)*(gt_peak_idx(cnt)-idx+1)/fs*1000, get(gca, 'YLim'), '-r');
        title(sprintf('RX signal %d', cnt))

        subplot(6,2,5);
        plot(seg_idx/fs, abs(analogData_pn(seg_idx)), '-b');
        hold on;
        plot(gt_peak_idx(cnt)/fs, abs(analogData_pn(gt_peak_idx(cnt))), 'ro');
        xlabel('time (s)');
        title('PN signal')

        subplot(6,2,7);
        plot(seg_idx/fs, corr(seg_idx), '-b');
        hold on;
        plot(gt_peak_idx(cnt)/fs, corr(gt_peak_idx(cnt)), 'ro');
        xlabel('time (s)'); ylabel('corr');
        title('PN: corr')

        subplot(6,2,[6,8]);
        spectrogram(analogData_pn(seg_idx), 256, 64, 256, fs, 'yaxis');
        hold on;
        plot(ones(1,2)*(gt_peak_idx(cnt)-idx+1)/fs*1000, get(gca, 'YLim'), '-r');
        ylim([0 5]);
        title('PN');

        subplot(6,2,9);
        plot(rx_chirp_time{chirp_cnt}/fs, abs(rx_chirp{chirp_cnt}), '-b');
        xlabel('time (s)');
        title('FMCW signal')

        subplot(6,2,11);
        plot(freq, this_fmcw_fft, '-b');
        hold on;
        plot(freq(gt_fmcw_peak_idx), this_fmcw_fft(gt_fmcw_peak_idx), 'ro');
        xlabel('frequency');
        title('FMCW: fft');

        subplot(6,2,[10,12]);
        spectrogram(rx_chirp{chirp_cnt}, 256, 64, 256, fs, 'yaxis');
        ylim([0 5]);
        title('FMCW')
        % print(fh, '-dpng', [fig_dir filename '.' combine_method '.spec.' num2str(cnt) '.png']);


        %%%%%%%%%%
        pn_cen_idx = [max(1,gt_peak_idx(cnt)+pn_dist_range(1)):min(T,gt_peak_idx(cnt)+pn_dist_range(2))];
        fmcw_cen_idx = [max(1,gt_fmcw_peak_idx+fmcw_dist_range(1)):min(T,gt_fmcw_peak_idx+fmcw_dist_range(2))];
        this_corr_dist = (pn_cen_idx-gt_peak_idx(cnt))/fs*c;
        this_fmcw_dist = fmcw_dist(fmcw_cen_idx)-gt_fmcw_peak_dist;

        fh = figure(12); clf;
        subplot(4,1,1);
        plot(pn_cen_idx/fs, corr(pn_cen_idx), '-b.');
        hold on;
        plot(gt_peak_idx(cnt)/fs, corr(gt_peak_idx(cnt)), 'ro');
        xlabel('time (ms)'); ylabel('corr');
        title('PN: corr');

        subplot(4,1,2);
        plot(this_corr_dist, corr(pn_cen_idx), '-b.');
        hold on;
        plot(0, corr(gt_peak_idx(cnt)), 'ro');
        xlabel('distance (cm)'); ylabel('corr');
        title('PN: corr');

        subplot(4,1,3);
        plot(freq(fmcw_cen_idx), this_fmcw_fft(fmcw_cen_idx), '-b.');
        hold on;
        plot(freq(gt_fmcw_peak_idx), this_fmcw_fft(gt_fmcw_peak_idx), 'ro');
        xlabel('frequency'); ylabel('corr');
        title('FMCW: fft');

        subplot(4,1,4);
        plot(this_fmcw_dist, this_fmcw_fft(fmcw_cen_idx), '-b.');
        hold on;
        plot(0, this_fmcw_fft(gt_fmcw_peak_idx), 'ro');
        xlabel('distance (cm)'); ylabel('corr');
        title('FMCW: fft');
        % print(fh, '-dpng', [fig_dir filename '.' combine_method '.dist.' num2str(cnt) '.png']);


        % status(cnt) = input('status of the signal (1: good, 2: bad):');
        if cnt <= 10
            status(cnt) = 1;
        else
            status(cnt) = 2;
        end




        %%%%%%%%%%
        % this_norm_corr = norm_corr(pn_cen_idx) / max(norm_corr(pn_cen_idx));
        % this_norm_fmcw_fft = this_fmcw_fft(fmcw_cen_idx) / max(this_fmcw_fft(fmcw_cen_idx));

        % this_range = [max(this_corr_dist(1), this_fmcw_dist(1)):dist_unit:min(this_corr_dist(end), this_fmcw_dist(end))];
        % this_range = unique(sort([this_range 0]));
        % corr_interp = interp1(this_corr_dist, this_norm_corr, this_range);
        % fmcw_interp = interp1(this_fmcw_dist, this_norm_fmcw_fft, this_range);
        % new_spikes = combine_pn_fmcw(corr_interp, fmcw_interp, combine_method);

        % [~,tmp1] = max(new_spikes);
        % new_error(cnt)  = abs(this_range(tmp1));
        % [~,tmp] = max(this_norm_corr);
        % pn_error(cnt)   = abs(this_corr_dist(tmp));
        % [~,tmp] = max(this_norm_fmcw_fft);
        % fmcw_error(cnt) = abs(this_fmcw_dist(tmp));
        % fprintf('  combined error = %.2fcm\n', new_error(cnt));
        % fprintf('  pn error = %.2fcm\n', pn_error(cnt));
        % fprintf('  fmcw error = %.2fcm\n', fmcw_error(cnt));

        % fh = figure(13); clf;
        % subplot(2,1,1);
        % lh(1) = plot(this_corr_dist, this_norm_corr, '-b.');
        % hold on;
        % plot(0, this_norm_corr(gt_peak_idx(cnt)-pn_cen_idx(1)+1), 'go');
        % lh(2) = plot(this_fmcw_dist, this_norm_fmcw_fft, '-r.');
        % plot(0, this_norm_fmcw_fft(gt_fmcw_peak_idx-fmcw_cen_idx(1)+1), 'yo');
        % legend(lh, {'PN', 'FMCW'});
        % xlabel('distance (cm)');
        % title(sprintf('signal %d', cnt));

        % subplot(2,1,2);
        % plot(this_range, new_spikes, '-b');
        % hold on;
        % plot(this_range(tmp1), new_spikes(tmp1), 'ro');
        % title('combine')
        % xlabel('distance (cm)');
        % print(fh, '-dpng', [fig_dir filename '.' combine_method '.combine.' num2str(cnt) '.png']);


        % waitforbuttonpress
        % return
        %% -------------------


        idx = idx + this_peak + floor(itvl_samp/2);
        cnt = cnt + 1;

    end


    dlmwrite([output_dir filename '.corr.txt'], corr', 'delimiter', '\t');
    dlmwrite([output_dir filename '.status.txt'], status', 'delimiter', '\t');
    dlmwrite([output_dir filename '.gt.txt'], gt_peak_idx', 'delimiter', '\t', 'precision', '%8d');



    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;

    % plot(new_error, '-ro'); hold on;
    % plot(pn_error, '-b+');
    % plot(fmcw_error, '-g*');
    % xlabel('PN index'); ylabel('error (cm)');
    % % ylim([0, 15]);
    % legend('combine', 'PN', 'FMCW');
    % % print(fh, '-dpng', [fig_dir filename '.' combine_method '.error.png']);

    % avg_errs = [mean(new_error), mean(pn_error), mean(fmcw_error)];
    % std_errs = [std(new_error), std(pn_error), std(fmcw_error)];

    % return
end




function [Nfft, Ncp, symbolRate, Ts, Fs, fc_pn, fc_fmcw, rollOff, nSamp, filterSpan, preambleSeq] = get_config_param(config)

    if config==1
        Nfft=128;
        Ncp=32;
        symbolRate=1/2205;
        Ts=1/44100;
        Fs=1/Ts;
        fc_pn = 12000;
        fc_fmcw = 5000;
        rollOff=0.25;
        nSamp=20;
        filterSpan=10;
        % preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
        preambleSeq=m_sequence([1 0 0 0 0 1 1 1 0 1 0 1 0 0 0 1 0]);
        % sampleInterval=0.16;
        % sampleInterval = 1/3*2 + 4263/Fs;
    elseif config==2
        Nfft=128;
        Ncp=32;
        symbolRate=1/2205;
        Ts=1/44100;
        Fs=1/Ts;
        fc_pn = 5000;
        fc_fmcw = 12000;
        rollOff=0.25;
        nSamp=20;
        filterSpan=10;
        % preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
        preambleSeq=m_sequence([1 0 0 0 0 1 1 1 0 1 0 1 0 0 0 1 0]);
        % sampleInterval=0.16;
        % sampleInterval = 1/3*2 + 4263/Fs;
    elseif config==3
        Nfft=128;
        Ncp=32;
        symbolRate=1/2205;
        Ts=1/44100;
        Fs=1/Ts;
        fc_pn = 5000;
        fc_fmcw = 12000;
        rollOff=0.25;
        nSamp=20;
        filterSpan=10;
        preambleSeq=m_sequence([1 0 0 0 0 1 1 1 0 1 0 1 0 0 0 1 0]);
    else
        error('wrong config');
    end
end


function [bw, tm, new_chirp_len] = get_fmcw_config_param(config)
    if config == 1
        bw = 2500;
        tm = 0.1;
        new_chirp_len = 4407.95;
    elseif config == 2
        bw = 2500;
        tm = 0.1;
        new_chirp_len = 4407.95;
    elseif config == 3
        bw = 2500;
        tm = 0.04;
        new_chirp_len = 4407.95 / 0.1 * 0.4;
    end
end


function [idx] = get_gt_first_peak_idx(filename, data, corr, courseStartIndex, windowSize)
    fs = 44100;

    % idx = find(corr > 0);
    % idx = idx(1);
    % range_idx = idx:courseStartIndex+windowSize+100;
    range_idx = courseStartIndex:courseStartIndex+4410;
    this_corr = corr(range_idx);
    this_data = abs(data(range_idx));

    fh = figure(101); clf;
    subplot(2,1,1)
    plot(this_data, '-b.'); hold on;
    plot(ones(1,2)*(courseStartIndex+windowSize-range_idx(1)+1), get(gca, 'ylim'), '-r')

    subplot(2,1,2)
    plot(this_corr, '-b.'); hold on;
    plot(ones(1,2)*(courseStartIndex+windowSize-range_idx(1)+1), get(gca, 'ylim'), '-r')
    % waitforbuttonpress

    if strcmp(filename, 'rx4')
        idx = 2968 + range_idx(1) - 1;
    % elseif strcmp(filename, 'rx5')
    %     idx = 3068 + range_idx(1) - 1;
    elseif strcmp(filename, 'rx6')
        idx = 2980 + range_idx(1) - 1;
    else
        idx = -1;
    end
end
