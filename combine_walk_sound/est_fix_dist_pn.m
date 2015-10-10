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
%%   est_fix_dist_pn('0903.exp4.fix.pn2047.1m', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0903.exp5.fix.pn2047.0.5m', 8000, 2047, 0.5, 1)
%%
%%   est_fix_dist_pn('0904.exp1.fix.pn2047.loc1', 8000, 2047, sqrt(2.8^2+1^2), 1)
%%   est_fix_dist_pn('0904.exp1.fix.pn2047.loc2', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0904.exp1.fix.pn2047.loc3', 8000, 2047, sqrt(3.1^2+1^2), 1)
%%   est_fix_dist_pn('0904.exp1.fix.pn2047.loc4', 8000, 2047, sqrt(3.1^2+2^2), 1)
%%   est_fix_dist_pn('0904.exp1.fix.pn2047.loc5', 8000, 2047, 2, 1)
%%   est_fix_dist_pn('0904.exp1.fix.pn2047.loc6', 8000, 2047, sqrt(2.8^2+2^2), 1)
%%
%%   est_fix_dist_pn('0906.exp4.fix.pn2047.loc1', 8000, 2047, sqrt(2.8^2+1^2), 1)
%%   est_fix_dist_pn('0906.exp4.fix.pn2047.loc2', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp4.fix.pn2047.loc3', 8000, 2047, sqrt(3.1^2+1^2), 1)
%%   est_fix_dist_pn('0906.exp4.fix.pn2047.loc4', 8000, 2047, sqrt(3.1^2+2^2), 1)
%%   est_fix_dist_pn('0906.exp4.fix.pn2047.loc5', 8000, 2047, 2, 1)
%%   est_fix_dist_pn('0906.exp4.fix.pn2047.loc6', 8000, 2047, sqrt(2.8^2+2^2), 1)
%%     
%%   est_fix_dist_pn('0906.exp1.1m.degree0', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp1.1m.degree45', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp1.1m.degree90', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp1.1m.degree135', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp1.1m.degree180', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp1.1m.degree225', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp1.1m.degree270', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp1.1m.degree315', 8000, 2047, 1, 1)
%%
%%   est_fix_dist_pn('0906.exp2.3m.wall.degree0', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp2.3m.wall.degree45', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp2.3m.wall.degree90', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp2.3m.wall.degree135', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp2.3m.wall.degree180', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp2.3m.wall.degree225', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp2.3m.wall.degree270', 8000, 2047, 1, 1)
%%   est_fix_dist_pn('0906.exp2.3m.wall.degree315', 8000, 2047, 1, 1)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [err] = est_fix_dist_pn(filename, freq, pn_len, init_dist, frame_len)
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
    if nargin < 3, pn_len = 2047; end
    if nargin < 4, init_dist = 1; end
    if nargin < 5, frame_len = 1; end


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
    % waitforbuttonpress
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

