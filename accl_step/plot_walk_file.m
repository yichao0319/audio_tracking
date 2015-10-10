%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%%
%% - Input:
%%     - marix format:
%%       1. locationHeadingTimestamp_since1970
%%       2. locationHeadingX
%%       3. locationHeadingY
%%       4. locationHeadingZ
%%       5. locationTrueHeading
%%       6. locationMagneticHeading
%%       7. locationHeadingAccuracy
%%       8. accelerometerTimestamp_sinceReboot
%%       9. accelerometerAccelerationX
%%       10. accelerometerAccelerationY
%%       11. accelerometerAccelerationZ
%%       12. gyroTimestamp_sinceReboot
%%       13. gyroRotationX
%%       14. gyroRotationY
%%       15. gyroRotationZ
%%
%% - Output:
%%
%%
%% example:
%%
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_walk_file()
    % addpath('../utils');
    
    %% --------------------
    %% DEBUG
    %% --------------------
    DEBUG0 = 0;
    DEBUG1 = 1;
    DEBUG2 = 1;  %% progress
    DEBUG3 = 1;  %% verbose
    DEBUG4 = 1;  %% results

    STFT = 0;


    %% --------------------
    %% Constant
    %% --------------------


    %% --------------------
    %% Variable
    %% --------------------
    input_dir  = './accl_mat/';
    fig_dir = './fig/';
    filename = '0820.exp1.accl.walk.50m';
    % filename = '0820.exp2.accl.walk.55m';
    % filename = '0820.exp3.accl.walk.45m';
    % filename = '0820.exp4.accl.walk.50m';
    % filename = '0820.exp5.accl.walk.square';
    fig_idx = 0;

    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, arg = 1; end
    if nargin < 1, arg = 1; end


    %% --------------------
    %% Main starts
    %% --------------------
    
    
    %% --------------------
    %% load data
    %% --------------------
    if DEBUG2, fprintf('load data\n'); end

    data = load([input_dir filename '.txt'])';
    % size(data)
    data(1,:) = data(1,:) - data(1,1);
    data(8,:) = data(8,:) - data(8,1);
    data(12,:) = data(12,:) - data(12,1);

    Fs = size(data,2) / (data(8,end)-data(8,1));
    fprintf('  Fs = %f\n', Fs);
    Fs = ceil(Fs);
    
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    subplot(3, 1, 1);
    plot(data(8,:), data(9,:));
    subplot(3, 1, 2);
    plot(data(8,:), data(10,:));
    subplot(3, 1, 3);
    plot(data(8,:), data(11,:));

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    subplot(3, 1, 1);
    plot(data(12,:), data(13,:));
    subplot(3, 1, 2);
    plot(data(12,:), data(14,:));
    subplot(3, 1, 3);
    plot(data(12,:), data(15,:));

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    subplot(5, 1, 1);
    plot(data(1,:), data(2,:));
    subplot(5, 1, 2);
    plot(data(1,:), data(3,:));
    subplot(5, 1, 3);
    plot(data(1,:), data(4,:));
    subplot(5, 1, 4);
    plot(data(1,:), data(5,:));
    subplot(5, 1, 5);
    plot(data(1,:), data(6,:));


    %% --------------------
    %% Remove DC (zero frequency component)
    %% --------------------
    if DEBUG2, fprintf('Remove DC\n'); end

    [B,A] = butter(10, 0.15/(Fs/2), 'high');
    for fi = 2:size(data,1)
        tmp = filter(B,A,data(fi,:));
        data(fi,:) = tmp;
        % data(fi,:) = data(fi,:) - mean(data(fi,:));
    end


    %% --------------------
    %% STFT
    %% --------------------
    if STFT
        if DEBUG2, fprintf('STFT\n'); end
        
        window = Fs*5;
        noverlap = floor(window/2);
        Nfft = Fs*2000;
        % Nfft = window;
        
        % Spectrogram takes the STFT of the signal
        % P matrix contains the power spectral density of each segment of the STFT
        [S,F,T,P] = spectrogram(data(9,:), window, noverlap, Nfft, Fs);

        fig_idx = fig_idx + 1;
        fh = figure(fig_idx); clf;
        imagesc(T, F, 10*log10(P)); % frequency-time Plot of the signal
        colorbar;
        % ylim([f_min f_max]);
        xlabel('Time (s)');
        ylabel('Power/Frequency (dB/Hz)');
        print(fh, '-dpsc', [fig_dir filename '.psd.eps']);
    end


    ts = data(9,:);
    ts_len = length(ts);
    [step_idx, periodic_idx, autocorr] = get_step_idx(ts, [], fig_idx);

    
    %% plot ts and corrcoef
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    subplot(2,1,1);
    plot(ts);
    hold on;
    plot(periodic_idx, ts(periodic_idx), 'ro');
    hold on;
    plot(step_idx, ts(step_idx), 'g.');
    set(gca, 'XLim', [1 ts_len]);
    subplot(2,1,2);
    plot(autocorr);
    hold on;
    plot(periodic_idx, autocorr(periodic_idx), 'ro');
    % hold on;
    % plot(step_idx, autocorr(step_idx), 'g.');
    set(gca, 'XLim', [1 ts_len]);
    print(fh, '-dpsc', [fig_dir filename '.xcorr.eps']);

end






