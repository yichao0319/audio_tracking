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
%%   ofdm_rx_itvl_longpn('ofdm.1.12000.128', 128, 1)
%%   ofdm_rx_itvl_longpn('rx1', 128, 1)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [norm_first_peak_idx] = ofdm_rx_itvl_longpn(filename, max_pn_len, config)
    % addpath('../utils');

    %% --------------------
    %% DEBUG
    %% --------------------
    DEBUG0 = 0;
    DEBUG1 = 1;
    DEBUG2 = 1;  %% progress
    DEBUG3 = 1;  %% verbose
    DEBUG4 = 1;  %% results

    OPT_PEAKS = 0;

    if nargin < 2, max_pn_len = 128; end
    if nargin < 3, config = 1; end


    %% --------------------
    %% Constant
    %% --------------------
    tx_dir  = './tx_sound/';
    input_dir  = './rx_sound/';
    fig_dir = './fig/';
    font_size = 18;


    %% --------------------
    %% Variable
    %% --------------------
    fig_idx = 0;


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, filename = 'ofdm.18000'; end
    if nargin < 2, config = 1; end

    [Nfft, Ncp, symbolRate, Ts, fs, fc, rollOff, nSamp, filterSpan, preambleSeq] = get_config_param(config);
    pn_len = Nfft*floor(max_pn_len/Nfft);
    fprintf('  PN length = %d\n', max_pn_len);

    preamble_filename = [tx_dir 'preamble.' num2str(config) '.' num2str(fc) '.' num2str(pn_len) '.mat'];
    tx_filename = [tx_dir 'txch_1signal.' num2str(config) '.' num2str(fc) '.' num2str(pn_len) '.mat'];


    %% --------------------
    %% Main starts
    %% --------------------

    %% ------------------
    %% DEBUG
    % all_peaks{1} = [1 5 10];
    % all_peaks{2} = [13 14 17];
    % all_peaks{3} = [21 25 27];
    % itvl = 10;
    % [min_err, min_idx] = select_peaks_minimize_err_rec(all_peaks, 1, 0, itvl);
    % min_err
    % min_idx
    % return;
    %% ------------------

    %% ====================
    %% Load data
    %% ====================
    if DEBUG2, fprintf('Load data\n'); end

    %% TX One Signal
    fprintf('  tx signal file: %s\n', tx_filename);
    tmp = load(tx_filename);
    tx_analogData = tmp.analog_one_symbol;
    itvl_samp = length(tx_analogData);
    fprintf('  tx signal size = %dx%d\n', size(tx_analogData));

    %% Preambles
    fprintf('  preamble file: %s\n', preamble_filename);
    load(preamble_filename);
    fprintf('  preamble size = %dx%d\n', size(preamble));

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
    else
        error('Cannot find the file');
    end
    fprintf('  audio file: %s\n', audio_filename);
    [analogData,~] = audioread(audio_filename);
    analogData = analogData.';



    %% ----------------
    %% XXX: need to modify
    % len = length(analogData) / fs
    if size(analogData,1) > 1
        analogData = analogData(2,:);
    end
    % analogData = analogData(:,(0*fs+1):min(end,100*fs));
    % analogData = analogData(:,(0.1*fs+1):min(end,itvl_samp*100));
    % analogData = analogData(:, 0*fs+1:min(end,4*fs));
    %% ----------------
    fprintf('  data size = %dx%d\n', size(analogData));



    %% ====================
    %% downconvert
    %% ====================
    if DEBUG2, fprintf('Downconvert\n'); end

    T = numel(analogData);
    analogData_orig = analogData;
    analogData = analogData .* exp(-1i*2*pi*fc*(1:T)*Ts);
    % analogDataFFT=fft(analogData);
    % plot(abs(analogDataFFT));

    if DEBUG2, fprintf('Low Pass Filter\n'); end
    analogData = lowPassFilterByFFT(analogData, fs, 1500, 0);


    %% ======================================================
    %% course timing synchornizing
    %% ======================================================
    if DEBUG2, fprintf('course timing synchornizing\n'); end

    windowSize = 3000;
    detectLength = ceil(itvl_samp * 1.5);
    courseStartIndex = findStartIndexByDoubleWin(analogData, windowSize, detectLength);
    % courseStartIndex = 15000;
    fprintf('  Course Start Index = %d (%fs)\n', courseStartIndex, courseStartIndex/fs);


    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(2,2,[1,2]);
    plot((1:T)/fs, analogData_orig, '-b');
    hold on;
    plot(ones(1,2)*courseStartIndex/fs, get(gca, 'YLim'), '-r');
    % set(gca, 'XLim', [0 3]);
    xlabel('Time (s)'); ylabel('Amplitude (v)');
    title('raw audio'); axis tight;

    subplot(2,2,3);
    spectrogram(analogData_orig(1:3*fs), 256, 64, 256, fs, 'yaxis');
    % set(gca, 'YLim', [0 5]);
    title('Before Downconvertion');

    subplot(2,2,4);
    spectrogram(analogData(1:3*fs), 256, 64, 256, fs, 'yaxis');
    hold on;
    set(gca, 'YLim', [0 22.05]);
    plot(ones(1,2)*courseStartIndex/fs, get(gca, 'YLim'), '-r');
    title('After Downconvertion');
    % return



    %% ======================================================
    %% calculate corrcoeff
    %% ======================================================
    if DEBUG2, fprintf('calculate corrcoeff\n'); end

    Np = numel(preamble);
    prev_idx = courseStartIndex;
    i = courseStartIndex;
    corr = zeros(1, length(analogData));
    norm_corr = corr;
    while i <= length(analogData)
        % fprintf('  %d/%d\n', i, length(analogData));
        % if mod(i-1, floor(length(analogData)/10)) == 0
        %     fprintf('%d\n', ceil((i-1)/length(analogData)*100));
        % end
        show_progress(i-courseStartIndex+1, length(analogData)-courseStartIndex+1, 1);

        if i+Np-1 > length(analogData)
            break;
        end

        corr(i) = abs(analogData(i:i+Np-1)*(preamble.'));

        if (i - prev_idx > (itvl_samp/4))
            idx = [prev_idx:i];
            norm_corr(idx) = corr(idx) / max(corr(idx));

            %% --------------------
            %% plot each signal and its correlation
            fh = figure(11); clf;

            subplot(4,1,1)
            idx = prev_idx:min(prev_idx+itvl_samp-1,length(corr));
            [~,max_idx] = max(corr(idx));
            max_idx = idx(1) + max_idx - 1;
            plot(idx/fs, abs(analogData(idx)), '-b.'); hold on;
            plot(max_idx/fs, abs(analogData(max_idx)), 'ro');
            title('analog data')

            subplot(4,1,2);
            plot(idx/fs, corr(idx), '-b.'); hold on;
            plot(max_idx/fs, corr(max_idx), 'ro');
            hold on;

            title('normalized corr')

            subplot(4,1,3)
            idx = max(1,max_idx-300):max_idx+300;
            plot(idx/fs, abs(analogData(idx)), '-b.'); hold on;
            plot(max_idx/fs, abs(analogData(max_idx)), 'ro');
            title('analog data')

            subplot(4,1,4);
            plot(idx/fs, corr(idx), '-b.'); hold on;
            plot(max_idx/fs, corr(max_idx), 'ro');
            title('corr')
            % waitforbuttonpress;
            %% --------------------

            i = prev_idx + itvl_samp - 1;
            prev_idx = i;

            continue;
        end

        i = i + 1;
    end


    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % subplot(2,1,1)
    % plot([1:length(analogData)]/fs, analogData);

    % subplot(2,1,2);
    % plot(corr);


    %% ======================================================
    %% Find peak
    %% ======================================================
    if DEBUG2, fprintf('Find peak\n'); end

    idx = courseStartIndex;
    %% corrcoeff
    peak_idx = [];
    first_peak_idx = [];
    %% normalized corrcoeff
    norm_peak_idx = [];
    norm_first_peak_idx = [];
    %% double window
    windowSize = 30;
    dwin_divisor = zeros(size(corr));
    dwin_peak_idx = [];
    %% Select Peaks which minimize error
    all_peaks = {};
    seg_num = 0;
    %% thresholds
    thresh = 0.01;
    % thresh = 0.14;
    norm_thresh = 0.95;
    while(idx + itvl_samp < length(corr))
        seg_idx = idx:idx+itvl_samp-1;
        seg = corr(seg_idx);

        %% max peak
        [v, this_peak] = max(seg);
        peak_idx = [peak_idx idx+this_peak-1];

        %% max peak of normalized corr
        norm_seg = norm_corr(seg_idx);
        [v, this_peak] = max(norm_seg);
        norm_peak_idx = [norm_peak_idx idx+this_peak-1];

        %% first peak > thresh
        [v, locs] = findpeaks(seg);
        tmp = locs(find(v>thresh));
        if length(tmp) > 0
            first_peak_idx = [first_peak_idx idx+tmp(1)-1];
            this_peak = tmp(1);
        end

        %% first peak of normalized corr > thresh
        [v, locs] = findpeaks(norm_seg);
        tmp = locs(find(v>norm_thresh));
        if length(tmp) > 0
            norm_first_peak_idx = [norm_first_peak_idx idx+tmp(1)-1];
            this_peak = tmp(1);
        end

        %% double window
        % std_idx = max(idx, 1+windowSize);
        std_idx = idx+windowSize+1;
        denom = sum(norm_corr(std_idx-windowSize:std_idx-1));
        nom = sum(norm_corr(std_idx+1:std_idx+windowSize));
        for ii = std_idx:idx+itvl_samp-1
            dwin_divisor(ii) = nom / denom;

            denom = denom - norm_corr(ii-windowSize) + norm_corr(ii);
            nom = nom - norm_corr(ii+1) + norm_corr(ii+windowSize+1);
        end
        dwin_divisor(dwin_divisor==Inf) = 0;
        % [v, this_dwin_peak] = max(dwin_divisor(std_idx:idx+itvl_samp-1));
        % dwin_peak_idx = [dwin_peak_idx std_idx+this_dwin_peak-1];
        [v, locs] = findpeaks(dwin_divisor(std_idx:idx+itvl_samp-1));
        [tmp, sort_v_idx] = sort(v, 'descend');
        if length(sort_v_idx) >= 2
            this_dwin_2nd_peak = locs(sort_v_idx(2));
            dwin_peak_idx = [dwin_peak_idx std_idx+this_dwin_2nd_peak-1];
        end

        %% select peaks which minimize error
        if OPT_PEAKS
            seg_num = seg_num + 1;
            [v, all_peaks{seg_num}] = findpeaks(norm_seg);
            all_peaks{seg_num} = idx+all_peaks{seg_num}-1;
            if length(all_peaks{seg_num}) == 0 && seg_num == 1
                seg_num = 0;
            end
        end


        %% --------------------
        %% plot
        fh = figure(12); clf;

        subplot(4,1,1)
        plot(seg_idx/fs, abs(analogData(seg_idx)), '-b.'); hold on;
        plot(peak_idx(end)/fs, abs(analogData(peak_idx(end))), 'mo');
        plot(first_peak_idx(end)/fs, abs(analogData(first_peak_idx(end))), 'm*');
        lh = plot(norm_peak_idx(end)/fs, abs(analogData(norm_peak_idx(end))), 'go');
        set(lh, 'MarkerSize', 1);
        lh = plot(norm_first_peak_idx(end)/fs, abs(analogData(norm_first_peak_idx(end))), 'g*');
        set(lh, 'MarkerSize', 1);
        title('analog data')

        subplot(4,1,2)
        plot(seg_idx/fs, corr(seg_idx), '-b.'); hold on;
        plot(peak_idx(end)/fs, corr(peak_idx(end)), 'mo');
        plot(first_peak_idx(end)/fs, corr(first_peak_idx(end)), 'm*');
        lh = plot(norm_peak_idx(end)/fs, corr(norm_peak_idx(end)), 'go');
        set(lh, 'MarkerSize', 1);
        lh = plot(norm_first_peak_idx(end)/fs, corr(norm_first_peak_idx(end)), 'g*');
        set(lh, 'MarkerSize', 1);
        title('corr')

        subplot(4,1,3)
        plot_seg_idx = [max(1,peak_idx(end)-2000):peak_idx(end)+2000];
        plot(plot_seg_idx/fs, corr(plot_seg_idx), '-b.'); hold on;
        plot(peak_idx(end)/fs, corr(peak_idx(end)), 'mo');
        plot(first_peak_idx(end)/fs, corr(first_peak_idx(end)), 'm*');
        lh = plot(norm_peak_idx(end)/fs, corr(norm_peak_idx(end)), 'go');
        set(lh, 'MarkerSize', 1);
        lh = plot(norm_first_peak_idx(end)/fs, corr(norm_first_peak_idx(end)), 'g*');
        set(lh, 'MarkerSize', 1);
        title('corr -- zoom in')

        subplot(4,1,4)
        plot(plot_seg_idx/fs, norm_corr(plot_seg_idx), '-b.'); hold on;
        plot(peak_idx(end)/fs, norm_corr(peak_idx(end)), 'mo');
        plot(first_peak_idx(end)/fs, norm_corr(first_peak_idx(end)), 'm*');
        lh = plot(norm_peak_idx(end)/fs, norm_corr(norm_peak_idx(end)), 'go');
        set(lh, 'MarkerSize', 1);
        lh = plot(norm_first_peak_idx(end)/fs, norm_corr(norm_first_peak_idx(end)), 'g*');
        set(lh, 'MarkerSize', 1);
        title('normalized corr -- zoom in')

        % waitforbuttonpress;
        %% --------------------


        idx = idx + this_peak + floor(itvl_samp/2);
        % idx = idx + this_peak + ceil(fs*4) - floor(itvl_samp/2);
    end


    %% ======================================================
    %% Deal wit the clock skew
    %% ======================================================
    if DEBUG2, fprintf('Deal with the clock skew'); end

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    plot(norm_peak_idx(2:end) - norm_peak_idx(1:end-1), '-b.'); hold on;
    plot((norm_peak_idx(2:end) - norm_peak_idx(1)) ./ (1:(length(norm_peak_idx)-1)), '-r*')
    grid('on');
    itvl_first = (norm_peak_idx(2:end) - norm_peak_idx(1)) ./ (1:(length(norm_peak_idx)-1));
    itvl_prev  = norm_peak_idx(2:end) - norm_peak_idx(1:end-1);
    fprintf('  new PN itvl: prev=%.2f, first = %.2f (last=%.2f)\n', mean(itvl_prev), mean(itvl_first), itvl_first(end));
    % new_itvl_samp = mean(itvl_prev);
    % new_itvl_samp = 44079.67;
    % new_itvl_samp = 44079.4;
    new_itvl_samp = 44079.5;


    num_segs = floor(length(corr) / new_itvl_samp);
    % gt_idx = norm_first_peak_idx(1)*ones(1, num_segs) + round([0:num_segs-1]*(itvl_samp));
    gt_idx = norm_peak_idx(1)*ones(1, num_segs) + round([0:num_segs-1]*(new_itvl_samp));
    while(gt_idx(end) > length(analogData))
        gt_idx = gt_idx(1:end-1);
    end


    %% ======================================================
    %% Select Peaks which minimize error
    %% ======================================================
    if OPT_PEAKS
        if DEBUG2, fprintf('Select Peaks which minimize error\n'); end

        % opt_peak_idx = select_peaks_minimize_err(all_peaks);
        [min_err, opt_peak_idx] = select_peaks_minimize_err_rec(all_peaks, 1, 0, itvl_samp, Inf);
    end


    %% ------------------------
    %% Plot figure
    %% ------------------------
    if DEBUG2, fprintf('Plot figure\n'); end

    idx = courseStartIndex;
    cnt = 1;
    while(idx + itvl_samp < length(corr))
        seg_idx = idx:idx+itvl_samp-1;

        fh = figure(13); clf;

        subplot(4,1,1)
        plot(seg_idx/fs, abs(analogData(seg_idx)), '-b.'); hold on;
        lh = plot(gt_idx(cnt)/fs, abs(analogData(gt_idx(cnt))), 'rx');
        set(lh, 'MarkerSize', 10);
        set(lh, 'LineWidth', 10);
        plot(peak_idx(cnt)/fs, abs(analogData(peak_idx(cnt))), 'mo');
        % plot(first_peak_idx(cnt)/fs, abs(analogData(first_peak_idx(cnt))), 'm*');
        lh = plot(norm_peak_idx(cnt)/fs, abs(analogData(norm_peak_idx(cnt))), 'go');
        set(lh, 'MarkerSize', 1);
        % lh = plot(norm_first_peak_idx(cnt)/fs, abs(analogData(norm_first_peak_idx(cnt))), 'g*');
        % set(lh, 'MarkerSize', 1);
        title('analog data')

        subplot(4,1,2)
        plot(seg_idx/fs, corr(seg_idx), '-b.'); hold on;
        lh = plot(gt_idx(cnt)/fs, corr(gt_idx(cnt)), 'rx');
        set(lh, 'MarkerSize', 10);
        set(lh, 'LineWidth', 10);
        plot(peak_idx(cnt)/fs, corr(peak_idx(cnt)), 'mo');
        % plot(first_peak_idx(cnt)/fs, corr(first_peak_idx(cnt)), 'm*');
        lh = plot(norm_peak_idx(cnt)/fs, corr(norm_peak_idx(cnt)), 'go');
        set(lh, 'MarkerSize', 1);
        % lh = plot(norm_first_peak_idx(cnt)/fs, corr(norm_first_peak_idx(cnt)), 'g*');
        % set(lh, 'MarkerSize', 1);
        title('corr')

        subplot(4,1,3)
        plot_seg_idx = [max(1,peak_idx(cnt)-100):peak_idx(cnt)+100];
        plot(plot_seg_idx/fs, corr(plot_seg_idx), '-b.'); hold on;
        lh = plot(gt_idx(cnt)/fs, corr(gt_idx(cnt)), 'rx');
        set(lh, 'MarkerSize', 10);
        set(lh, 'LineWidth', 10);
        plot(peak_idx(cnt)/fs, corr(peak_idx(cnt)), 'mo');
        % plot(first_peak_idx(cnt)/fs, corr(first_peak_idx(cnt)), 'm*');
        lh = plot(norm_peak_idx(cnt)/fs, corr(norm_peak_idx(cnt)), 'go');
        set(lh, 'MarkerSize', 1);
        % lh = plot(norm_first_peak_idx(cnt)/fs, corr(norm_first_peak_idx(cnt)), 'g*');
        % set(lh, 'MarkerSize', 1);
        title('corr -- zoom in')

        subplot(4,1,4)
        plot(plot_seg_idx/fs, norm_corr(plot_seg_idx), '-b.'); hold on;
        lh = plot(gt_idx(cnt)/fs, norm_corr(gt_idx(cnt)), 'rx');
        set(lh, 'MarkerSize', 10);
        set(lh, 'LineWidth', 10);
        plot(peak_idx(cnt)/fs, norm_corr(peak_idx(cnt)), 'mo');
        % plot(first_peak_idx(cnt)/fs, norm_corr(first_peak_idx(cnt)), 'm*');
        lh = plot(norm_peak_idx(cnt)/fs, norm_corr(norm_peak_idx(cnt)), 'go');
        set(lh, 'MarkerSize', 1);
        % lh = plot(norm_first_peak_idx(cnt)/fs, norm_corr(norm_first_peak_idx(cnt)), 'g*');
        % set(lh, 'MarkerSize', 1);
        title('normalized corr -- zoom in')

        waitforbuttonpress;
        %% --------------------


        idx = idx + this_peak + floor(itvl_samp/2);
        cnt = cnt + 1;
    end


    %% plot all PN seq
    if DEBUG2, fprintf('  plot all PN seq\n'); end
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    subplot(3,1,1);
    plot(abs(analogData));
    hold on;
    % lh = plot(peak_idx, abs(analogData(peak_idx)), 'ro');
    lh = plot(gt_idx, abs(analogData(gt_idx)), 'ro');
    set(lh, 'MarkerSize', 10);
    % lh = plot(first_peak_idx, abs(analogData(first_peak_idx)), 'g*');
    % set(lh, 'MarkerSize', 10);
    lh = plot(norm_first_peak_idx, abs(analogData(norm_first_peak_idx)), 'g^');
    set(lh, 'MarkerSize', 10);
    % lh = plot(dwin_peak_idx, abs(analogData(dwin_peak_idx)), 'ks');
    % set(lh, 'MarkerSize', 10);
    if OPT_PEAKS
        lh = plot(opt_peak_idx, abs(analogData(opt_peak_idx)), 'ks');
        set(lh, 'MarkerSize', 10);
    end
    set(gca, 'XLim', [1, length(analogData)]);
    set(gca, 'FontSize', font_size);
    title('Audio Signals', 'Color', 'r');

    subplot(3,1,2);
    plot(corr, '-b.');
    hold on;
    % lh = plot(peak_idx, corr(peak_idx), 'ro');
    lh = plot(gt_idx, corr(gt_idx), 'ro');
    set(lh, 'MarkerSize', 10);
    lh = plot(first_peak_idx, corr(first_peak_idx), 'g*');
    set(lh, 'MarkerSize', 10);
    set(gca, 'XLim', [1, length(analogData)]);
    set(gca, 'FontSize', font_size);
    title('CorrCoeff', 'Color', 'r');

    subplot(3,1,3);
    plot(norm_corr, '-b.');
    hold on;
    % lh = plot(norm_peak_idx, norm_corr(norm_peak_idx), 'ro');
    lh = plot(gt_idx, norm_corr(gt_idx), 'ro');
    set(lh, 'MarkerSize', 10);
    % lh = plot(first_peak_idx, norm_corr(first_peak_idx), 'g*');
    % set(lh, 'MarkerSize', 10);
    lh = plot(norm_first_peak_idx, norm_corr(norm_first_peak_idx), 'g^');
    set(lh, 'MarkerSize', 10);
    % lh = plot(dwin_peak_idx, norm_corr(dwin_peak_idx), 'ks');
    % set(lh, 'MarkerSize', 10);
    if OPT_PEAKS
        lh = plot(opt_peak_idx, norm_corr(opt_peak_idx), 'ks');
        set(lh, 'MarkerSize', 10);
    end
    set(gca, 'XLim', [1, length(analogData)]);
    % set(gca, 'XLim', [2, 2.2]*10^4);
    set(gca, 'FontSize', font_size);
    title('Normalized CorrCoeff', 'Color', 'r');

    % subplot(3,1,3);
    % plot(dwin_divisor, '-b.');
    % hold on;
    % % lh = plot(norm_peak_idx, dwin_divisor(norm_peak_idx), 'ro');
    % lh = plot(gt_idx, dwin_divisor(gt_idx), 'ro');
    % set(lh, 'MarkerSize', 10);
    % lh = plot(norm_first_peak_idx, dwin_divisor(norm_first_peak_idx), 'g^');
    % set(lh, 'MarkerSize', 10);
    % lh = plot(dwin_peak_idx, dwin_divisor(dwin_peak_idx), 'ks');
    % set(lh, 'MarkerSize', 10);
    % set(gca, 'XLim', [1, length(analogData)]);
    % % set(gca, 'XLim', [2, 2.2]*10^4);
    % set(gca, 'YLim', [0, 10]);
    % set(gca, 'FontSize', font_size);
    % title('Double Window Divisor', 'Color', 'r');


    %% Plot Distance
    if DEBUG2, fprintf('  plot distance\n'); end
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    dif = peak_idx(2:end) - peak_idx(1:end-1);
    dif_dist = cumsum(dif - new_itvl_samp)/fs * 34600;
    plot([0:length(peak_idx)-2]*new_itvl_samp/fs, dif_dist, '-r');
    hold on;

    first_dif = first_peak_idx(2:end) - first_peak_idx(1:end-1);
    first_dif_dist = cumsum(first_dif - new_itvl_samp)/fs * 34600;
    plot([0:length(peak_idx)-2]*new_itvl_samp/fs, first_dif_dist, '-b');

    norm_dif = norm_first_peak_idx(2:end) - norm_first_peak_idx(1:end-1);
    norm_dif_dist = cumsum(norm_dif - new_itvl_samp)/fs * 34600;
    plot([0:length(peak_idx)-2]*new_itvl_samp/fs, norm_dif_dist, '-m');

    dwin_dif = dwin_peak_idx(2:end) - dwin_peak_idx(1:end-1);
    % dwin_dif_dist = cumsum(dwin_dif - new_itvl_samp)/fs * 34600;
    % plot([0:length(peak_idx)-2]*new_itvl_samp/fs, dwin_dif_dist, '-g');

    title('displacement');
    xlabel('time (s)'); ylabel('distance (cm)');
    legend('largest', 'first', 'normalized first');

    if OPT_PEAKS
        opt_dif = opt_peak_idx(2:end) - opt_peak_idx(1:end-1);
    end

    return
    %% ======================================================


end


function [Nfft, Ncp, symbolRate, Ts, Fs, fc, rollOff, nSamp, filterSpan, preambleSeq] = get_config_param(config)

    if config==1
        Nfft=128;
        Ncp=32;
        symbolRate=1/2205;
        Ts=1/44100;
        Fs=1/Ts;
        fc=12000;
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
        fc=16000;
        rollOff=0.25;
        nSamp=20;
        filterSpan=10;
        % preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
        preambleSeq=m_sequence([1 0 0 0 0 1 1 1 0 1 0 1 0 0 0 1 0]);
    else
        error('wrong config');
    end
end


%% select_peaks_minimize_err: function description
function [sel_peaks] = select_peaks_minimize_err(peaks)

end


function [min_err, min_idx] = select_peaks_minimize_err_rec(peaks, tidx, idx, itvl, err)
    if tidx == length(peaks)
        for yi = 1:length(peaks{tidx})
            this_err(yi) = abs((peaks{tidx}(yi) - idx) - itvl);
        end

        [min_err, min_err_idx] = min(this_err);
        min_idx = peaks{tidx}(min_err_idx);
        return;
    end

    this_err = [];
    for yi = 1:length(peaks{tidx})
        if tidx > 1
            this_err(yi) = abs((peaks{tidx}(yi) - idx) - itvl);
        else
            % fprintf('  %d/%d -- %d\n', yi, length(peaks{1}), length(peaks));
            show_progress(yi, length(peaks{1}), 1);
            this_err(yi) = 0;
        end

        if yi > 1 && this_err(yi) > err
            this_err(yi) = Inf;
        else
            [this_min_err, this_min_idx{yi}] = select_peaks_minimize_err_rec(peaks, tidx+1, peaks{tidx}(yi), itvl, err);
            this_err(yi) = this_err(yi) + this_min_err;

            if this_err(yi) < err
                err = this_err(yi);
            end
        end
    end

    [min_err, min_err_idx] = min(this_err);
    min_idx = [peaks{tidx}(min_err_idx) this_min_idx{min_err_idx}];

end

