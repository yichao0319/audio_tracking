%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%% MAIN
%% - Input:
%%
%%
%% - Output:
%%
%%
%% example:
%%   visualize_sound_pn('0908.exp1.pn2047.fix.1m.noblock', 8000, 2047, 1)
%%   visualize_sound_pn('0908.exp2.pn2047.fix.1m.block', 8000, 2047, 1)
%%   visualize_sound_pn('0913.exp1.pn2047.fix.self', 7000, 2047, 1, 'wav')
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visualize_sound_pn(filename, freq, pn_len, init_dist, ext)
    % addpath('../utils');
    
    %% --------------------
    %% DEBUG
    %% --------------------
    DEBUG0 = 0;
    DEBUG1 = 1;
    DEBUG2 = 1;  %% progress
    DEBUG3 = 1;  %% verbose
    DEBUG4 = 1;  %% results

    warning('off','comm:commsrc:pn:GenPolyNotPrimitive');


    %% --------------------
    %% Constant
    %% --------------------
    input_dir = './raw/';
    fig_dir = './fig/';

    sound_speed = 331;

    frame_len = 1;
    thresh = 0.8;

    fig_idx = 0;
    font_size = 28;




    %% --------------------
    %% Variable
    %% --------------------
    


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, filename = 'pn.8000.1'; end
    if nargin < 2, freq = 8000; end
    if nargin < 3, pn_len = 2047; end
    if nargin < 4, init_dist = 2.7; end
    if nargin < 5, ext = 'aac'; end


    %% --------------------
    %% Main starts
    %% --------------------
    %% ====================================
    %% Get sound info
    %% ====================================
    % [pn_len, freq, frame_len, dist] = get_sound_info(filename);

    
    %% ====================================
    %% Read the audio file
    %% ====================================
    if DEBUG2, fprintf('Read file\n'); end

    % ext = 'aac';
    % if strcmp(filename, 'tx'), ext = 'wav'; end
    file_path_name = [input_dir filename '.' ext];
    fprintf('  file: %s\n', file_path_name)
    [wav_data, Fs] = audioread(file_path_name);
    Ts = 1/Fs;
    nbits = 16;


    %% ----------------------------
    %% only use part of the trace
    % wav_data = wav_data(1:min(end,20*Fs), :);
    %% ----------------------------

    wav_len = length(wav_data);
    wav_time = 0:1/Fs:(wav_len-1)/Fs;
    if DEBUG4, 
        fprintf('- wav_data %d x %d (Fs=%d, nbits=%d)\n', size(wav_data), Fs, nbits); 
        fprintf('  duration = %f\n', wav_time(end));
    end
    
    if strcmp(filename, 'tx')
        wav_data = wav_data(:,2);
    else
        wav_data = wav_data(:,1);
    end

    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % plot(wav_time, wav_data);
    % title('sound')
    % print(fh, '-dpsc', [fig_dir filename '.sound.eps']);



    %% ====================================
    %% Down-convert
    %% ====================================
    % if DEBUG2, fprintf('Down-convert\n'); end

    % wav_data_base = wav_data .* sin(2*pi*freq*Ts * [1:length(wav_data)]');
    % wav_data_base = lowPassFilterByFFT(wav_data_base', Fs, 3000, 10)';


    %% ====================================
    %% generate base sound
    %%   PN sequence
    %% ====================================
    if DEBUG2, fprintf('generate PN sequence\n'); end
    
    %% --------------
    %% Method 1
    if pn_len == 511
        fbconnection = [0 0 0 0 0 1 0 0 1];
    elseif pn_len == 1023
        fbconnection = [0 0 0 0 0 0 1 0 0 1];
    elseif pn_len == 2047
        fbconnection = [0 0 0 0 0 0 0 1 0 0 1];
    elseif pn_len == 4095
        fbconnection = [0 0 0 0 0 0 0 0 1 0 0 1];
    elseif pn_len == 8191
        fbconnection = [0 0 0 0 0 0 0 0 0 1 0 0 1];
    else
        error('wrong number of bits');
    end
    pnseq = m_sequence(fbconnection)';

    %% --------------
    %% Method 2
    % h = commsrc.pn('GenPoly', [[8 2 0]], 'Mask', [1 0 0 0 0 0 1 0]);
    % set(h, 'NumBitsOut', pn_len);
    % pnseq = generate(h);
    % for k=1:length(pnseq)
    %     if(pnseq(k) == 0)
    %         pnseq(k) = -1;
    %     end
    % end
    
    sound_sample = pnseq;
    fprintf('  sound_sample: %d x %d\n', size(sound_sample));

    % up-convert
    sound_sample = sound_sample .* sin(2*pi*freq*Ts * [1:length(sound_sample)])';
    n_sound = length(sound_sample);

    
    %% ====================================
    %% Find xcorr
    %% ====================================
    if DEBUG2, fprintf('Find xcorr\n'); end

    maxpeak = xcorr(sound_sample, wav_data);
    maxpeak = maxpeak(1:floor(end/2)+1);
    maxpeak = -maxpeak;
    % size(0:1/Fs:length(maxpeak))
    % plot(wav_time, maxpeak)
    % set(gca, 'XLim', [0 10]);
    

    [v, idx] = max(maxpeak(1:1*frame_len*Fs));

    % fh = figure(1); clf;
    % plot(wav_time, maxpeak);
    % hold on;
    % plot(wav_time(idx), maxpeak(idx), 'ro');
    % set(gca, 'XLim', [0 10]);
    % return


    %% ====================================
    %% Find intervals
    %% ====================================
    if DEBUG2, fprintf('Find intervals\n'); end

    win_std = max(1,idx - floor(frame_len*Fs/2));
    win_end = win_std + frame_len*Fs - 1;
    peaks_idx = [];
    while(1)
        [v, idx] = max(maxpeak(win_std:win_end));
        peaks_idx = [peaks_idx, win_std + idx - 1];

        win_std = win_end + 1;
        win_end = win_std + frame_len*Fs - 1;
        if win_end > length(maxpeak)
            break;
        end
    end
    
    itvl_pre = peaks_idx(2:end) - peaks_idx(1:end-1);

    fprintf('  itvl to the previous signal: \n');
    fprintf('    %.2f\n', itvl_pre);
    fprintf('\n\n');

    itvl_first = (peaks_idx(2:end) - peaks_idx(1)) ./ [1:length(peaks_idx)-1];
    fprintf('  itvl to the first signal: \n');
    fprintf('    %.2f\n', itvl_first);
    fprintf('\n');


    range = [min(itvl_pre):0.2:max(itvl_pre)];
    hist = histc(itvl_pre, range);
    hist = hist / sum(hist);
    
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    bh1 = bar(range/1000, hist);
    xlabel('Interval (K samples)', 'FontSize', font_size)
    ylabel('Ratio', 'FontSize', font_size);
    tmpx = get(gca,'XTick');
    % set(gca, 'XLim', [range(1)-0.1 range(end)+0.1]);
    % set(gca, 'XLim', [44098 44102]);
    % set(gca,'XTick', [range(1):(range(end)-range(1))/10:range(end)]);
    % set(gca,'XTickLabel',sprintf('%.1f|', range));
    set(gca, 'FontSize', font_size);
    print(fh, '-dpsc', [fig_dir filename '.itvl2pre.eps']);

    % [f,x] = ecdf(itvl_pre);
    
    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % lh = plot(x/1000, f, '-b.');
    % set(lh, 'LineWidth', 2);
    % set(lh, 'MarkerSize', 20);
    % xlabel('Interval (*1000 samples)', 'FontSize', font_size)
    % ylabel('CDF', 'FontSize', font_size);
    % tmpx = get(gca,'XTick');
    % set(gca, 'FontSize', font_size);
    % print(fh, '-dpsc', [fig_dir filename '.itvl2pre.eps']);


    %% ====================================
    %% learn interval length
    %% ====================================
    if DEBUG2, fprintf('learn interval length\n'); end

    num_train = 10;
    train_itvls = (peaks_idx(2:num_train) - peaks_idx(1)) ./ [1:num_train-1] / Fs;
    rcv_itvl = median(train_itvls);
    % std(train_itvls)
    % mean(train_itvls)
    % rcv_itvl
    fprintf('  learned interval: %f (%.2f)\n', rcv_itvl, rcv_itvl*Fs);
    % return


    %% ====================================
    %% learn initial time
    %% ====================================
    if DEBUG2, fprintf('learn initial time\n'); end

    num_train = 10;
    dist = init_dist;
    init_times = peaks_idx(1:num_train) / Fs - [0:num_train-1] * rcv_itvl - dist/sound_speed;
    init_time = median(init_times);
    % std(init_times)
    % mean(init_times)
    % init_time


    %% ====================================
    %% Plot xcorr
    %% ====================================
    if DEBUG2, fprintf('Plot xcorr\n'); end

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    exp_peaks_idx = [peaks_idx(1) peaks_idx(2:end-1)];
    exp_peaks_idx(2:end) = exp_peaks_idx(2:end) + 44100;
    
    % xlim = [0.95 1.05];
    % xlim = [2.25 2.4];

    subplot(2, 1, 1)
    lh = plot(wav_time, wav_data);
    set(lh, 'LineWidth', 1);
    xlabel('Time (s)', 'FontSize', font_size);
    ylabel('RX Sound', 'FontSize', font_size);
    set(gca, 'FontSize', font_size);
    % set(gca, 'XLim', xlim);

    subplot(2, 1, 2)
    plot(wav_time, maxpeak, '-b.');
    hold on;
    lh = plot(wav_time(peaks_idx), maxpeak(peaks_idx), 'ro');
    set(lh, 'LineWidth', 4);
    set(lh, 'MarkerSize', 10);
    lh = plot(wav_time(exp_peaks_idx), maxpeak(exp_peaks_idx), 'gx');
    set(lh, 'LineWidth', 4);
    set(lh, 'MarkerSize', 10);
    xlabel('Time (s)', 'FontSize', font_size);
    ylabel('xcorr', 'FontSize', font_size);
    set(gca, 'FontSize', font_size);
    % set(gca, 'XLim', xlim);


    return


    %% ====================================
    %% calculate distance
    %% ====================================
    if DEBUG2, fprintf('calculate distance\n'); end

    dists = (peaks_idx / Fs - init_time - [0:length(peaks_idx)-1]*rcv_itvl) * sound_speed;
    dist_err = abs(dists - dist);
    mean(dist_err) * 100
    
    time_len = wav_time(end);
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    
    subplot(3, 1, 1)
    lh = plot([0:length(dists)-1]*rcv_itvl, dists);
    set(lh, 'LineWidth', 4);
    set(gca, 'XLim', [0 time_len]);
    % set(gca, 'YLim', [min(dists)-5 max(dists)+5]);
    xlabel('Time (s)', 'FontSize', font_size);
    ylabel('Distance (m)', 'FontSize', font_size);
    set(gca, 'FontSize', font_size);

    subplot(3, 1, 2)
    lh = plot([0:length(dists)-1]*rcv_itvl, dist_err);
    set(lh, 'LineWidth', 4);
    set(gca, 'XLim', [0 time_len]);
    xlabel('Time (s)', 'FontSize', font_size);
    ylabel('Distance (m)', 'FontSize', font_size);
    set(gca, 'FontSize', font_size);

    subplot(3, 1, 3)
    lh = plot(wav_time, maxpeak);
    set(lh, 'LineWidth', 4);
    hold on;
    plot(wav_time(peaks_idx), maxpeak(peaks_idx), 'ro');
    set(gca, 'XLim', [0 time_len]);
    xlabel('Time (s)', 'FontSize', font_size);
    ylabel('xcorr', 'FontSize', font_size);
    set(gca, 'FontSize', font_size);

    print(fh, '-dpsc', [fig_dir filename '.sound_dist.eps']);
    
    
    
    % mean(dist_err) * 100
    % [p, v] = ecdf(dist_err);

    % fh = figure; clf;
    % plot(v, p)
end


%% 08.10.pn.511.8000.1.dist30
function [pn_len, freq, frame_len, dist] = get_sound_info(filename)
    expression = ['.*\.pn\.(?<pn_len>\d+)\.(?<freq>\d+)\.(?<frame_len>\d+)\.dist(?<dist>\d+\.*\d*)'];
    tokenNames = regexp(filename, expression, 'names');
    if(length(tokenNames) > 0)
        % freq = str2num(tokenNames(1).freq) * 1000;
        pn_len = str2num(tokenNames(1).pn_len);
        freq = str2num(tokenNames(1).freq);
        dist = str2num(tokenNames(1).dist);
        frame_len = str2num(tokenNames(1).frame_len);
        
        fprintf('  pn len = %f\n', pn_len);
        fprintf('  freq = %d\n', freq);
        fprintf('  dist = %f\n', dist);
        fprintf('  frame_len = %f\n', frame_len);
    else
        error(['wong format: ' filename]);
    end
end