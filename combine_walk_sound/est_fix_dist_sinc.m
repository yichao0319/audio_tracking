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
%%   est_fix_dist_sinc('0903.exp7.fix.sinc.1m', 8000, 1, 1)
%%   est_fix_dist_sinc('0903.exp6.fix.sinc.0.5m', 8000, 0.5, 1)
%%   est_fix_dist_sinc('0904.exp2.fix.sinc.loc1', 8000, sqrt(2.8^2+1^2), 1)
%%   est_fix_dist_sinc('0904.exp2.fix.sinc.loc2', 8000, 1, 1)
%%   est_fix_dist_sinc('0904.exp2.fix.sinc.loc3', 8000, sqrt(3.1^2+1^2), 1)
%%   est_fix_dist_sinc('0904.exp2.fix.sinc.loc4', 8000, sqrt(3.1^2+2^2), 1)
%%   est_fix_dist_sinc('0904.exp2.fix.sinc.loc5', 8000, 2, 1)
%%   est_fix_dist_sinc('0904.exp2.fix.sinc.loc6', 8000, sqrt(2.8^2+2^2), 1)
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [err] = est_fix_dist_sinc(filename, freq, init_dist, frame_len)
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
    input_dir = './sound/';
    fig_dir = './fig/';

    sound_speed = 331;

    thresh = 0.8;

    fig_idx = 0;
    font_size = 16;




    %% --------------------
    %% Variable
    %% --------------------
    


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, filename = 'pn.8000.1'; end
    if nargin < 2, freq = 8000; end
    if nargin < 3, init_dist = 1; end
    if nargin < 4, frame_len = 1; end


    %% --------------------
    %% Main starts
    %% --------------------
    
    %% ====================================
    %% Read the audio file
    %% ====================================
    if DEBUG2, fprintf('Read file\n'); end

    ext = 'aac';
    file_path_name = [input_dir filename '.' ext];
    [wav_data, Fs] = audioread(file_path_name);
    Ts = 1/Fs;
    nbits = 16;


    %% ----------------------------
    %% only use part of the trace
    % wav_data = wav_data(1.4*Fs:end, :);
    %% ----------------------------

    wav_len = length(wav_data);
    wav_time = 0:1/Fs:(wav_len-1)/Fs;
    if DEBUG4, 
        fprintf('- wav_data %d x %d (Fs=%d, nbits=%d)\n', size(wav_data), Fs, nbits); 
        fprintf('  duration = %f\n', wav_time(end));
    end
    
    wav_data = wav_data(:,1);

    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % plot(wav_time, wav_data);
    % print(fh, '-dpsc', [fig_dir filename '.sound.eps']);

    % return;


    %% ====================================
    %% Find intervals
    %% ====================================
    [v, idx] = max(wav_data(1:1*frame_len*Fs));

    win_std = max(1,idx - floor(frame_len*Fs/2));
    win_end = win_std + frame_len*Fs - 1;
    peaks_idx = [];
    while(1)
        [v, idx] = max(wav_data(win_std:win_end));
        peaks_idx = [peaks_idx, win_std + idx - 1];

        win_std = win_end + 1;
        win_end = win_std + frame_len*Fs - 1;
        if win_end > length(wav_data)
            break;
        end
    end


    %% ====================================
    %% Calculate base, interval, and t1
    %% ====================================
    peaks_idx(2:end) - peaks_idx(1:end-1)
    [base, t1, itvl] = cal_t1_itvl(peaks_idx, init_dist, Fs, frame_len);
    

    %% ====================================
    %% Calculate distance
    %% d_i = (t_i' - t_b - (i-b)*itvl) * v
    %% ====================================
    dists = (peaks_idx/Fs - t1 - ([1:length(peaks_idx)]-base)*itvl) * sound_speed
    median(dists)
    mean(dists)
    err = abs(dists - init_dist)
end

