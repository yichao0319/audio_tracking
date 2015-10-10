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
%%  visualize_beep('0828.exp1')
%%  visualize_beep('0828.exp2')
%%  visualize_beep('0901.exp1')
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visualize_beep(filename)
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
    input_dir  = './raw/';
    output_dir = '';

    sound_speed = 331;

    % freq = 17000;
    freqs = [8000, 6000];
    amp = 2;
    frame_len = 1;
    thresh = 0.8;

    pn_len = 2047;

    font_size = 16;


    %% --------------------
    %% Variable
    %% --------------------
    fig_idx = 0;


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, filename = 'tmp'; end
    % if nargin < 1, arg = 1; end


    %% --------------------
    %% Main starts
    %% --------------------

    %% --------------------
    %% Read Audio
    %% --------------------
    if DEBUG2, fprintf('Read Audio\n'); end

    file1 = [input_dir filename '.pc1.wav'];
    [data1, Fs] = audioread(file1);
    Ts = 1/Fs;
    nbits = 16;
    fprintf('  file = %s\n', file1);
    fprintf('  size = %dx%d\n', size(data1));

    file2 = [input_dir filename '.pc2.wav'];
    [data2, Fs] = audioread([input_dir filename '.pc2.wav']);
    Ts = 1/Fs;
    nbits = 16;
    fprintf('  file = %s\n', file2);
    fprintf('  size = %dx%d\n', size(data2));


    
    %% --------------------
    %% generate base sound
    %%   PN sequence
    %% --------------------
    if DEBUG2, fprintf('Generate PN sequence\n'); end
    
    h = commsrc.pn('GenPoly', [[8 2 0]], 'Mask', [1 0 0 0 0 0 1 0]);
    set(h, 'NumBitsOut', pn_len);
    pnseq = generate(h);
    
    for k=1:length(pnseq)
        if(pnseq(k) == 0)
            pnseq(k) = -1;
        end
    end
    
    % up-convert
    sound_sample1 = pnseq .* sin(2*pi*freqs(1)*Ts * [1:length(pnseq)])';
    sound_sample2 = pnseq .* sin(2*pi*freqs(2)*Ts * [1:length(pnseq)])';

    
    %% --------------------
    %% Cross Correlation
    %% --------------------
    if DEBUG2, fprintf('Cross Correlation\n'); end

    % ts_corr{1,1} = xcorr(sound_sample1, data1);
    % ts_corr{1,1} = ts_corr{1,1}(1:floor(end/2)+1);

    % ts_corr{1,2} = xcorr(sound_sample2, data1);
    % ts_corr{1,2} = ts_corr{1,2}(1:floor(end/2)+1);

    % ts_corr{2,1} = xcorr(sound_sample1, data2);
    % ts_corr{2,1} = ts_corr{2,1}(1:floor(end/2)+1);

    % ts_corr{2,2} = xcorr(sound_sample2, data2);
    % ts_corr{2,2} = ts_corr{2,2}(1:floor(end/2)+1);

    ts_corr{1,1} = my_xcorr(data1, sound_sample1);
    ts_corr{1,2} = my_xcorr(data1, sound_sample2);
    ts_corr{2,1} = my_xcorr(data2, sound_sample1);
    ts_corr{2,2} = my_xcorr(data2, sound_sample2);


    %% --------------------
    %% Find Peaks
    %% --------------------
    if DEBUG2, fprintf('Find Peaks\n'); end

    for i = 1:size(ts_corr,1)
        for j = 1:size(ts_corr,2)
            fprintf('  mic%d, spk%d\n', i, j);

            ti = 1;
            peak_idx{i,j} = [];
            while(ti < length(ts_corr{i,j}))
                ts_seg = ts_corr{i,j}(ti:min(ti+frame_len*Fs-1,end));
                [v,idx] = max(ts_seg);
                idx = idx + ti - 1;
                peak_idx{i,j} = [peak_idx{i,j} idx];

                ti = ti + frame_len*Fs;
            end

            % size(peak_idx{i,j})
            peak_idx{i,j}(2:end) - peak_idx{i,j}(1:end-1)
        end
    end

    

    %% --------------------
    %% Plot
    %% --------------------
    if DEBUG2, fprintf('Plot\n'); end

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(3,1,1);
    plot(data1);
    hold on;
    plot(peak_idx{1,1}, data1(peak_idx{1,1}), 'ro');
    plot(peak_idx{1,2}, data1(peak_idx{1,2}), 'gx');
    subplot(3,1,2);
    plot(ts_corr{1,1});
    subplot(3,1,3);
    plot(ts_corr{1,2});
    

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(3,1,1);
    plot(data2);
    hold on;
    plot(peak_idx{2,1}, data2(peak_idx{2,1}), 'ro');
    plot(peak_idx{2,2}, data2(peak_idx{2,2}), 'gx');
    subplot(3,1,2);
    plot(ts_corr{2,1});
    subplot(3,1,3);
    plot(ts_corr{2,2});


    %% --------------------
    %% Calculate distance
    %% --------------------
    if DEBUG2, fprintf('Calculate distance\n'); end
    
    for i = 1:size(peak_idx, 1)
        itvl{i} = abs(peak_idx{i,1} - peak_idx{i,2});
        % itvl{i}
    end

    (itvl{1} - itvl{2}) / Fs * sound_speed / 2
end



