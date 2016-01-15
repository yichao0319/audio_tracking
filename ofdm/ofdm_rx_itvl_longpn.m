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
%%   ofdm_rx_itvl_longpn('1102.exp1', 256, 1)
%%   ofdm_rx_itvl_longpn('1102.exp2', 128, 1)
%%   ofdm_rx_itvl_longpn('1103.exp1', 3000, 1)
%%   ofdm_rx_itvl_longpn('1103.exp2', 3000, 1)
%%   ofdm_rx_itvl_longpn('1103.exp3', 3000, 1)
%%   ofdm_rx_itvl_longpn('1104.exp1', 128, 1)
%%   ofdm_rx_itvl_longpn('1111.exp1', 128, 1)
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ofdm_rx_itvl_longpn(filename, max_pn_len, config)
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
    tx_dir  = './gen_data/';
    input_dir  = './rx_sound/';
    % input_dir  = './gen_data/';
    % output_dir = '';
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

    [Nfft, Ncp, symbolRate, Ts, Fs, fc, rollOff, nSamp, filterSpan, preambleSeq, sampleInterval] = get_config_param(config);
    sampleDelay = filterSpan*nSamp/2+nSamp*Ncp;
    % max_pn_len = 300;
    % max_pn_len = 128;
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
    audio_filename = [input_dir filename '.wav'];
    if exist(audio_filename, 'file') ~= 2,
        audio_filename = [input_dir filename '.aac'];
    end
    fprintf('  audio file: %s\n', audio_filename);
    [analogData,~] = audioread(audio_filename);
    analogData = analogData.';

    

    %% ----------------
    %% XXX: need to modify
    % len = length(analogData) / Fs
    if size(analogData,1) > 1
        analogData = analogData(2,1:10*Fs);
    end
    % analogData = analogData(:,(0*Fs+1):min(end,100*Fs));
    % analogData = analogData(:,(0.1*Fs+1):min(end,itvl_samp*100));
    analogData = analogData(:, 0*Fs+1:min(end,100*Fs));
    %% ----------------
    fprintf('  data size = %dx%d\n', size(analogData));

    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % plot(analogData);
    % % hold on;
    % % plot(courseStartIndex, (analogData(courseStartIndex)), 'ro');
    % title('raw audio');
    % return


    %% ====================
    %% downconvert
    %% ====================
    if DEBUG2, fprintf('Downconvert\n'); end

    T = numel(analogData);
    analogData = analogData .* exp(-1i*2*pi*fc*(1:T)*Ts);
    % analogDataFFT=fft(analogData);
    % plot(abs(analogDataFFT));

    if DEBUG2, fprintf('Low Pass Filter\n'); end
    analogData=lowPassFilterByFFT(analogData,Fs,2000,500);

    % if DEBUG2, fprintf('FFT\n'); end
    % analogDataFFT=fft(analogData);

    % course timing synchornizing
    if DEBUG2, fprintf('course timing synchornizing\n'); end
    windowSize = 300;
    detectLength = ceil(itvl_samp * 1.5);
    courseStartIndex = findStartIndexByDoubleWin(analogData,windowSize,detectLength);

    %% ======================================================
    %% calculate corrcoeff
    Np=numel(preamble);
    % for i = 1:length(analogData)
    %     if mod(i-1, floor(length(analogData)/10)) == 0
    %         fprintf('%d\n', ceil((i-1)/length(analogData)*100));
    %     end
    %     if i+Np-1 > length(analogData)
    %         break;
    %     end
    %     corr(i) = abs(analogData(i:i+Np-1)*preamble');
    % end
    prev_idx = courseStartIndex;
    % for i = courseStartIndex:length(analogData)
    i = courseStartIndex;
    corr = zeros(1, length(analogData));
    norm_corr = corr;
    while i <= length(analogData)
        % fprintf('  %d/%d\n', i, length(analogData));
        if mod(i-1, floor(length(analogData)/10)) == 0
            fprintf('%d\n', ceil((i-1)/length(analogData)*100));
        end

        if i+Np-1 > length(analogData)
            break;
        end

        if (i - prev_idx > (itvl_samp/4))
            idx = [prev_idx:1:i];
            norm_corr(prev_idx:i) = corr(prev_idx:i) / max(corr(prev_idx:i));
            i = prev_idx + itvl_samp - 1;
            prev_idx = i;
            continue;
        end

        corr(i) = abs(analogData(i:i+Np-1)*preamble');
        i = i + 1;
    end

    %% ------------------------
    %% Find peak
    %% ------------------------
    if DEBUG2, fprintf('Find peak\n'); end

    idx = 1;
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
    thresh = 0.03;
    % thresh = 0.14;
    norm_thresh = 0.3;
    while(idx + itvl_samp < length(corr))
        seg = corr(idx:idx+itvl_samp-1);

        %% max peak
        [v, this_peak] = max(seg);
        peak_idx = [peak_idx idx+this_peak-1];

        %% max peak of normalized corr
        norm_seg = norm_corr(idx:idx+itvl_samp-1);
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


        idx = idx + this_peak + floor(itvl_samp/2);
    end

    %% get ground truth
    num_segs = floor(length(corr) / itvl_samp);
    gt_idx = norm_first_peak_idx(1)*ones(1, num_segs) + round([0:num_segs-1]*(itvl_samp-0.5));
    while(gt_idx(end) > length(analogData))
        gt_idx = gt_idx(1:end-1);
    end


    %% ------------------------
    %% Select Peaks which minimize error
    %% ------------------------
    if OPT_PEAKS
        if DEBUG2, fprintf('Select Peaks which minimize error\n'); end

        % opt_peak_idx = select_peaks_minimize_err(all_peaks);
        [min_err, opt_peak_idx] = select_peaks_minimize_err_rec(all_peaks, 1, 0, itvl_samp, Inf);
    end


    %% ------------------------
    %% Plot figure
    %% ------------------------
    if DEBUG2, fprintf('Plot figure\n'); end

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
    

    dif = peak_idx(2:end) - peak_idx(1:end-1)
    first_dif = first_peak_idx(2:end) - first_peak_idx(1:end-1)
    norm_dif = norm_first_peak_idx(2:end) - norm_first_peak_idx(1:end-1)
    dwin_dif = dwin_peak_idx(2:end) - dwin_peak_idx(1:end-1)
    disp('normalization');
    norm_dif - itvl_samp
    % disp('double window');
    % dwin_dif - itvl_samp
    if OPT_PEAKS
        opt_dif = opt_peak_idx(2:end) - opt_peak_idx(1:end-1)
        disp('opt window');
        opt_dif - itvl_samp
    end


    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;

    % ranges = [itvl_samp-50:itvl_samp+50];
    % bincnts = histc(norm_dif, ranges);
    % bincnts = bincnts / sum(bincnts);
    % lh = plot([itvl_samp itvl_samp], [min(bincnts) max(bincnts)], '-r');
    % set(lh, 'LineWidth', 2);
    % hold on;
    % plot(ranges, bincnts, '-b.');
    % set(gca, 'XLim', [ranges(1) ranges(end)]);
    % set(gca, 'FontSize', font_size);
    % xlabel('Peak Interval (samples)', 'FontSize', font_size);
    % ylabel('PDF', 'FontSize', font_size);
    
    return
    %% ======================================================


end


function [Nfft, Ncp, symbolRate, Ts, Fs, fc, rollOff, nSamp, filterSpan, preambleSeq, sampleInterval] = get_config_param(config)
    
    if config==1
        Nfft=128;
        Ncp=32;
        symbolRate=1/2205;
        Ts=1/44100;
        Fs=1/Ts;
        fc=18000;
        rollOff=0.25;
        nSamp=20;
        filterSpan=10;
        % preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
        preambleSeq=m_sequence([1 0 0 0 0 1 1 1 0 1 0 1 0 0 0 1 0]);
        % sampleInterval=0.16;
        sampleInterval = 1/3*2 + 4263/Fs;
    elseif config==2
        Nfft=256;
        Ncp=128;
        symbolRate=1/4410;
        Ts=1/44100;
        Fs=1/Ts;
        fc=18000;
        rollOff=0.25;
        nSamp=10;
        filterSpan=10;
        preambleSeq=m_sequence([0 0 0 0 1 0 0 0 1]);
    elseif config==3
        Nfft=64;
        Ncp=16;
        symbolRate=1/2205;
        Ts=1/44100;
        Fs=1/Ts;
        fc=18000;
        rollOff=0.25;
        nSamp=20;
        filterSpan=10;
        preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
        sampleInterval=0.08;
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
            fprintf('  %d/%d -- %d\n', yi, length(peaks{1}), length(peaks));
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

