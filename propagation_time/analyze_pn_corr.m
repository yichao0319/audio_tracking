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
%%  analyze_pn_corr('./rx_sound/', '08.10.pn.511.8000.1.dist30', 'rx')
%%  analyze_pn_corr('./rx_sound/', '08.10.pn.4095.8000.1.dist30', 'rx')
%%  analyze_pn_corr('./rx_sound/', '08.10.pn.511.8000.1.dist50', 'rx')
%%  analyze_pn_corr('./rx_sound/', '08.10.pn.1023.8000.1.dist50', 'rx')
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function analyze_pn_corr(input_dir, filename, txrx)
    % addpath('../utils');
    
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
    sound_speed = 331;

    freq = 8000;
    amp = 2;
    frame_len = 1;
    thresh = 0.8;
    font_size = 16;


    %% --------------------
    %% Variable
    %% --------------------
    output_dir = './fig/';


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, input_dir = './tx_sound/'; end
    if nargin < 2, filename = 'pn.8000.1'; end
    if nargin < 3, txrx = 'tx'; end


    %% --------------------
    %% Main starts
    %% --------------------
    %% ====================================
    %% Get sound info
    %% ====================================
    [pn_len, freq, frame_len, dist] = get_sound_info(filename);

    
    %% ====================================
    %% Read the audio file
    %% ====================================
    if DEBUG2, fprintf('Read file\n'); end

    if strcmp(txrx, 'tx')
        ext = 'wav';
    else
        ext = 'aac';
    end
    file_path_name = [input_dir filename '.' ext];
    file_path_name
    [wav_data, Fs] = audioread(file_path_name);
    Ts = 1/Fs;
    nbits = 16;


    %% ----------------------------
    %% only use part of the trace
    wav_data = wav_data(1.4*Fs:end, :);
    %% ----------------------------

    wav_len = length(wav_data);
    wav_time = 0:1/Fs:(wav_len-1)/Fs;
    if DEBUG4, 
        fprintf('- wav_data %d x %d (Fs=%d, nbits=%d)\n', size(wav_data), Fs, nbits); 
        fprintf('  duration = %f\n', wav_time(end));
    end
    
    if strcmp(txrx, 'tx')
        wav_data = wav_data(:,2);
    else
        wav_data = wav_data(:,1);
    end


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
    
    % sound_sample = m_sequence([0 0 0 0 0 0 1 0 0 1])';
    % % up-convert
    % sound_sample_pass = sound_sample .* sin(2*pi*freq*Ts * [1:length(sound_sample)])';
    % n_sound = length(sound_sample);
    h = commsrc.pn('GenPoly', [[8 2 0]], 'Mask', [1 0 0 0 0 0 1 0]);
    set(h, 'NumBitsOut', pn_len);
    pnseq = generate(h);
    
    for k=1:length(pnseq)
        if(pnseq(k) == 0)
            pnseq(k) = -1;
        end
    end
    
    sound_sample = pnseq;
    fprintf('  sound_sample: %d x %d\n', size(sound_sample));

    % up-convert
    sound_sample = sound_sample .* sin(2*pi*freq*Ts * [1:length(sound_sample)])';
    n_sound = length(sound_sample);

    
    maxpeak = xcorr(sound_sample, wav_data);
    maxpeak = maxpeak(1:floor(end/2)+1);
    % size(0:1/Fs:length(maxpeak))
    % plot(wav_time, maxpeak)
    % set(gca, 'XLim', [0 10]);
    

    [v, idx] = max(maxpeak(1:5*frame_len*Fs));

    % fh = figure(1); clf;
    % plot(wav_time, maxpeak);
    % hold on;
    % plot(wav_time(idx), maxpeak(idx), 'ro');
    % set(gca, 'XLim', [0 10]);
    % return


    %% ====================================
    %% Find intervals
    %% ====================================
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

    % plot(wav_time, maxpeak);
    % hold on;
    % plot(wav_time(peaks_idx), maxpeak(peaks_idx), 'ro');
    % set(gca, 'XLim', [0 30]);


    %% ====================================
    %% learn interval length
    %% ====================================
    num_train = 10;
    train_itvls = (peaks_idx(2:num_train) - peaks_idx(1)) ./ [1:num_train-1] / Fs;
    rcv_itvl = median(train_itvls);
    % std(train_itvls)
    % mean(train_itvls)
    % rcv_itvl


    %% ====================================
    %% learn initial time
    %% ====================================
    num_train = 10;
    init_dist = dist;
    init_times = peaks_idx(1:num_train) / Fs - [0:num_train-1] * rcv_itvl - dist/sound_speed;
    init_time = median(init_times);
    % std(init_times)
    % mean(init_times)
    % init_time


    %% ====================================
    %% calculate distance
    %% ====================================
    dists = (peaks_idx / Fs - init_time - [0:length(peaks_idx)-1]*rcv_itvl) * sound_speed;
    dist_err = abs(dists - dist);
    mean(dist_err) * 100
    
    time_len = 120;
    fh = figure; clf;
    subplot(3, 1, 1)
    lh = plot([0:length(dists)-1]*rcv_itvl, dists);
    set(lh, 'LineWidth', 4);
    set(gca, 'XLim', [0 time_len]);
    set(gca, 'YLim', [min(dists)-5 max(dists)+5]);
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

    print(fh, '-dpng', [output_dir filename '.png']);
    
    
    
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