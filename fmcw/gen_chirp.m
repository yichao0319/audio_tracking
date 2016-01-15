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
%%   gen_chirp(15000, 10)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function gen_chirp(fc, time_len)

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
    input_dir  = '';
    output_dir = './tx_sound/';


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, fc = 17000; end
    if nargin < 2, time_len = 60; end


    %% --------------------
    %% Variable
    %% --------------------
    fig_idx = 0;

    % c = 300;
    % lambda = c/fc;
    % range_max = 10;
    % % tm = 5.5*range2time(range_max,c);
    % tm = 0.1;
    % range_res = 0.01;
    % % bw = range2bw(range_res, c);
    % bw = 2500;
    % sweep_slope = bw/tm;
    % fr_max = range2beat(range_max,sweep_slope,c);

    % fd_max = speed2dop(5,lambda);
    % fb_max = fr_max + fd_max;

    c = 300;
    lambda = c/fc;
    range_max = 10;
    % tm = 5.5*range2time(range_max,c);
    tm = 0.1;
    range_res = 0.01;
    % bw = range2bw(range_res, c);
    bw = 3000;
    sweep_slope = bw/tm;
    fr_max = range2beat(range_max,sweep_slope,c);

    fd_max = speed2dop(5,lambda);
    fb_max = fr_max + fd_max;

    % fs = max(2*fb_max,bw);
    fs = 44100;




    %% --------------------
    %% Main starts
    %% --------------------

    times = [0:1/fs:tm-1/fs]';
    %% --------------------
    %% Generate Chirp -- Method 1
    %% --------------------
    hw = phased.FMCWWaveform('SweepBandwidth', bw,...
       'SampleRate',fs,'SweepDirection','Up','SweepTime',tm, ...
       'NumSweeps',1);
    x_base = step(hw);
    x = x_base .* exp(1i*2*pi*fc*(1:(tm*fs))/fs)';


    %% --------------------
    %% Generate Chirp -- Method 1
    %% --------------------
    % x = chirp(times,fc,tm,fc+bw);

    fh = figure(1); clf;
    subplot(211); plot(times,real(x));
    xlabel('Time (s)'); ylabel('Amplitude (v)');
    title('FMCW signal'); axis tight;

    subplot(212); spectrogram(x,256,64,256,fs,'yaxis');
    % set(gca, 'ylim', [0, 5000])
    % set(gca, 'ylim', [0, 5000]+fc-1000)
    title('FMCW signal spectrogram');


    num_chirp = floor(time_len / tm);
    time_len = num_chirp * tm;
    fprintf('  #chirps=%d\n', num_chirp);

    %% scaling
    chirp_signal = x/max(abs(x))*0.9;
    chirp_signal_base = x_base/max(abs(x_base))*0.9;

    %% repeating
    analogData = repmat(chirp_signal, num_chirp, 1);


    %% save file
    filename = sprintf('tx_chirp.%d.B%d.T%.2f', fc, bw, tm);
    save([output_dir filename '.mat'], 'chirp_signal_base');

    signal_d = zeros(time_len*fs, 2);
    signal_d(:,2) = analogData;
    audiowrite([output_dir filename '.wav'], signal_d, fs, 'BitsPerSample', 16);
end

