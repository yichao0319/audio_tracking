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
%%   ofdm_rx_itvl_cum('1111.exp1', 128, 1)
%%   ofdm_rx_itvl_cum('1111.exp2', 128, 1)
%%   ofdm_rx_itvl_cum('1104.exp1', 128, 1)
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [orig_max, cum_max, err_orig_max, err_orig_median, err_orig_avg, err_cum] = ofdm_rx_itvl_cum(filename, max_pn_len, config, audio_len)
    % addpath('../utils');
    
    %% --------------------
    %% DEBUG
    %% --------------------
    DEBUG0 = 0;
    DEBUG1 = 1;
    DEBUG2 = 1;  %% progress
    DEBUG3 = 1;  %% verbose
    DEBUG4 = 1;  %% results

    DWIN_PEAK = 0;
    OPT_PEAKS = 0;

    if nargin < 2, max_pn_len = 128; end
    if nargin < 3, config = 1; end
    if nargin < 4, audio_len = 10000; end


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
    % analogData = analogData(:,(0.1*Fs+1):min(end,itvl_samp*12.4));
    % analogData = analogData(:, 0*Fs+1:min(end,100*Fs));
    analogData = analogData(:,(0.1*Fs+1):min(end,itvl_samp*audio_len));
    % analogData = analogData(:,500:min(end,itvl_samp*audio_len));
    %% ----------------
    %% start from a relative quite time
    if DEBUG2, fprintf('Start from a relative quite time\n'); end
    std_idx = find_1st_quiet_zone(analogData);
    analogData = analogData(std_idx:end);
    
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
    analogData = analogData .* exp(-1i*2*pi*(fc)*(1:T)*Ts);
    % analogDataFFT=fft(analogData);
    % plot(abs(analogDataFFT));

    if DEBUG2, fprintf('Low Pass Filter\n'); end
    analogData=lowPassFilterByFFT(analogData,Fs,2000,500);

    % if DEBUG2, fprintf('FFT\n'); end
    % analogDataFFT=fft(analogData);

    % course timing synchornizing
    if DEBUG2, fprintf('course timing synchornizing\n'); end
    windowSize = 300;
    detectLength = min(ceil(itvl_samp * 1.5), length(analogData));
    courseStartIndex = findStartIndexByDoubleWin(analogData,windowSize,detectLength);


    %% ====================
    %% calculate corrcoeff
    %% ====================
    Np = numel(preamble);
    prev_idx = courseStartIndex;
    i = courseStartIndex;
    corr = zeros(1, length(analogData));
    norm_corr = corr;
    max_peak_mag = [];
    % progress = ones(1,10);
    while i <= length(analogData)
        %% -------------
        show_progress(i, length(analogData), courseStartIndex);
        %% -------------

        if i+Np-1 > length(analogData)
            break;
        end

        if (i - prev_idx > (itvl_samp/4))
            %% normalize corr
            idx = [prev_idx:1:i];
            norm_corr(prev_idx:i) = corr(prev_idx:i) / max(corr(prev_idx:i));

            %% find the magnitude of max peak
            [v, idx] = max(corr(prev_idx:i));
            max_peak_mag = [max_peak_mag v];

            %% point to the next segment
            i = prev_idx + itvl_samp - 1;
            prev_idx = i;
            continue;
        end

        corr(i) = abs(analogData(i:i+Np-1)*preamble');
        i = i + 1;
    end


    %% ====================
    %% Find peak
    %% ====================
    if DEBUG2, fprintf('Find peak\n'); end

    %% normalized corrcoeff
    norm_first_peak_idx = [];
    norm_thresh = 0.3;
    
    idx = 1;
    while(idx + itvl_samp < length(corr))
        norm_seg = norm_corr(idx:idx+itvl_samp-1);
        
        %% first peak of normalized corr > thresh
        [v, locs] = findpeaks(norm_seg);
        tmp = locs(find(v>norm_thresh));
        if length(tmp) > 0
            norm_first_peak_idx = [norm_first_peak_idx idx+tmp(1)-1];
            this_peak = tmp(1);
        else
            this_peak = idx + floor(itvl_samp/2);
        end

        idx = idx + this_peak + floor(itvl_samp/2);
    end
    peak_idx = norm_first_peak_idx;

    
    %% ====================
    %% Correct the signal interval
    %% ====================
    if DEBUG2, fprintf('Correct the signal interval\n'); end
    
    [itvl_samp, base_sig] = correct_itvl_sample(peak_idx);

    
    %% ====================
    %% Get ground truth peak index
    %% ====================
    if DEBUG2, fprintf('Get ground truth peak index\n'); end

    gt_idx = round(peak_idx(base_sig) + ([1:length(peak_idx)]-base_sig)*itvl_samp);
    peak_err = abs(gt_idx - peak_idx);
    fprintf('  avg error = %2.2f\n', mean(peak_err));
    fprintf('  max = %d, min = %d\n', max(peak_err), min(peak_err));

    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % plot(peak_err);


    %% ====================
    %% Correct frequency offset
    %% ====================
    if DEBUG2, fprintf('Correct frequency offset\n'); end
    
    % [analogData] = correct_freq_offset(analogData, gt_idx, preamble);
    % return


    %% ====================
    %% Accumulate Signals
    %% ====================
    if DEBUG2, fprintf('Accumulate Signals\n'); end

    num_segs = floor(length(analogData) / itvl_samp);
    cum_analogData = zeros(1, round(itvl_samp));
    for si = 1:num_segs
        cum_analogData = cum_analogData + analogData(round((si-1)*itvl_samp+1):round((si-1)*itvl_samp)+round(itvl_samp));
        % cum_analogData = cum_analogData + abs(analogData(round((si-1)*itvl_samp+1):round((si-1)*itvl_samp)+round(itvl_samp)));
    end


    %% ====================
    %% calculate corrcoeff of Accumulated Signals
    %% ====================
    Np = numel(preamble);
    cum_corr = zeros(1, length(analogData));
    for i = 1:length(cum_analogData)
        if i+Np-1 > length(cum_analogData)
            break;
        end
        cum_corr(i) = abs(cum_analogData(i:i+Np-1)*preamble');
        % r = corrcoef(cum_analogData(i:i+Np-1), abs(preamble));
        % cum_corr(i) = r(1,2);
    end

    [cum_max_peak_mag, cum_max_peak_idx] = max(cum_corr);
    cum_gt_idx = mod(gt_idx(1), round(itvl_samp));
    fprintf('  peak mag = %.4f\n', cum_max_peak_mag);
    fprintf('  peak idx = %d, gt idx = %d -- err = %d\n', cum_max_peak_idx, cum_gt_idx, abs(cum_max_peak_idx-cum_gt_idx));

    norm_cum_corr = [cum_corr / max(cum_corr)];
    [v, locs] = findpeaks(norm_cum_corr);
    tmp = locs(find(v>norm_thresh));
    if length(tmp) > 0
        cum_first_peak_idx = tmp(1);
        cum_first_peak_mag = cum_corr(cum_first_peak_idx)
    end
    fprintf('  first peak mag = %.4f\n', cum_first_peak_mag);
    fprintf('  first peak idx = %d, gt idx = %d -- err = %d\n', cum_first_peak_idx, cum_gt_idx, abs(cum_first_peak_idx-cum_gt_idx));


    %% ====================
    %% Plot figure
    %% ====================
    if DEBUG2, fprintf('Plot figure\n'); end

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    subplot(4,1,1);
    plot(abs(analogData));
    hold on;
    plot(gt_idx, abs(analogData(gt_idx)), 'ro');
    set(gca, 'XLim', [1, length(analogData)]);
    set(gca, 'FontSize', font_size);
    title('Original Audio Signals', 'Color', 'r');
    
    subplot(4,1,2);
    plot(corr, '-b.');
    hold on;
    plot(gt_idx, corr(gt_idx), 'ro');
    plot(peak_idx, corr(peak_idx), 'g^');
    set(gca, 'XLim', [1, length(analogData)]);
    set(gca, 'YLim', [0, max(corr)]);
    set(gca, 'FontSize', font_size);
    title('CorrCoeff', 'Color', 'r');
    
    subplot(4,1,3);
    plot(abs(cum_analogData));
    hold on;
    plot(cum_gt_idx, abs(cum_analogData(cum_gt_idx)), 'ro');
    plot(cum_first_peak_idx, abs(cum_analogData(cum_first_peak_idx)), 'g^');
    set(gca, 'XLim', [1, length(cum_analogData)]);
    set(gca, 'FontSize', font_size);
    title('Accumulated Audio Signals', 'Color', 'r');
    
    subplot(4,1,4);
    plot(cum_corr, '-b.'); 
    hold on;
    plot(cum_gt_idx, cum_corr(cum_gt_idx), 'ro');
    plot(cum_first_peak_idx, cum_corr(cum_first_peak_idx), 'g^');
    set(gca, 'XLim', [1, length(cum_analogData)]);
    set(gca, 'YLim', [0, max(cum_corr)]);
    set(gca, 'FontSize', font_size);
    title('CorrCoeff', 'Color', 'r');
    

    orig_max = mean(max_peak_mag);
    cum_max = cum_max_peak_mag;
    fprintf('  avg orig peak mag: %f\n', orig_max);
    fprintf('  cum peak mag: %f\n', cum_max);

    err_orig_max = max(peak_err);
    err_orig_median = median(peak_err);
    err_orig_avg = mean(peak_err);
    err_cum = abs(cum_first_peak_idx-cum_gt_idx)

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

