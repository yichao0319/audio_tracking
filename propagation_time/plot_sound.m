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
%%  plot_sound('./tmp/', 'sinc.8000.1')
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_sound(input_dir, filename, txrx)
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


    %% --------------------
    %% Variable
    %% --------------------
    % input_dir  = '';
    output_dir = './fig/';
    freq = 8000;
    font_size = 12;
    frame_len = 0.1;


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, input_dir = './tx_sound/'; end
    if nargin < 2, filename = 'sinc.8000.1'; end
    if nargin < 3, txrx = 'tx'; end


    %% --------------------
    %% Main starts
    %% --------------------

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
    [wav_data, Fs] = audioread(file_path_name);
    Ts = 1/Fs;
    nbits = 16;


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

    wav_data_down = wav_data .* sin(2*pi*freq*Ts * [1:length(wav_data)]');
    wav_data_down = lowPassFilterByFFT(wav_data_down', Fs, 3000, 10)';
    
    % wav_data_down = [];
    % f_samples = ceil(frame_len*Fs);
    % for fi = 1:ceil(length(wav_data)/frame_len/Fs)
    %     wav_tmp = wav_data((fi-1)*f_samples+1:min(fi*f_samples,end));
    %     wav_tmp_down = wav_tmp .* sin(2*pi*freq*Ts * [1:length(wav_tmp)]');
    %     wav_tmp_down = lowPassFilterByFFT(wav_tmp_down', Fs, 3000, 10)';
    %     wav_data_down = [wav_data_down; wav_tmp_down];
    % end


    %% ====================================
    %% plot wav in time domain
    %% ====================================
    fh = figure(1);
    clf;

    if strcmp(txrx, 'tx')
        subplot(3,1,1);
        plot(wav_time, wav_data);
        set(gca, 'XLim', [0 10]);
        ylabel('Amplitude', 'FontSize', font_size);
        xlabel('Time (s)', 'FontSize', font_size);
        set(gca, 'FontSize', font_size);
        title('tx signals');

        subplot(3,1,2);
        plot(wav_time, wav_data_down, 'r');
        % set(gca, 'XLim', [0.11 0.13]);
        set(gca, 'XLim', [0.99 1.03]);
        set(gca, 'YLim', [min(wav_data_down) Inf]);
        ylabel('Amplitude', 'FontSize', font_size);
        xlabel('Time (s)', 'FontSize', font_size);
        set(gca, 'FontSize', font_size);
        title('1 pulse (baseband)');

        subplot(3,1,3);
        plot(wav_time, wav_data, 'r');
        set(gca, 'XLim', [0.99 1.03]);
        ylabel('Amplitude', 'FontSize', font_size);
        xlabel('Time (s)', 'FontSize', font_size);
        set(gca, 'FontSize', font_size);
        title('1 pulse (passband)');
    else
        subplot(3,1,1);
        plot(wav_time, wav_data);
        % set(gca, 'XLim', [5 6]);
        set(gca, 'XLim', [5 10]);
        ylabel('Amplitude', 'FontSize', font_size);
        xlabel('Time (s)', 'FontSize', font_size);
        set(gca, 'FontSize', font_size);
        title('rx signals');

        subplot(3,1,2);
        plot(wav_time, wav_data, 'r');
        % set(gca, 'XLim', [5.56 5.58]);
        set(gca, 'XLim', [5 10]);
        ylabel('Amplitude', 'FontSize', font_size);
        xlabel('Time (s)', 'FontSize', font_size);
        set(gca, 'FontSize', font_size);
        title('1 pulse (passband)');

        subplot(3,1,3);
        plot(wav_time, wav_data_down, 'r');
        % set(gca, 'XLim', [5.56 5.58]);
        set(gca, 'XLim', [5 10]);
        set(gca, 'YLim', [min(wav_data_down) Inf]);
        ylabel('Amplitude', 'FontSize', font_size);
        xlabel('Time (s)', 'FontSize', font_size);
        set(gca, 'FontSize', font_size);
        title('1 pulse (baseband)');
    end

    
    print(fh, '-dpsc', [output_dir filename '.time.ps']);
    

    return


    %% ====================================
    %% Short-Time Fourier Transform
    %% ====================================
    if DEBUG2, fprintf('Short-time Fourier transform\n'); end

    % window = Fs/2; % Should be minimum twice the maximum frequency we want to analyze
    window = floor(Fs/32);
    % window = floor(Fs/128);
    noverlap = floor(window/4); % 75% overlap
    Nfft = Fs;
    % Nfft = window;
    
    % Spectrogram takes the STFT of the signal
    % P matrix contains the power spectral density of each segment of the STFT
    [S,F,T,P] = spectrogram(wav_data, window, noverlap, Nfft, Fs);


    %% ====================================
    %% Plotting Power Spectral Density
    %% ====================================
    if DEBUG4, fprintf('- plot P: %d x %d\n', size(P)); end;

    fh = figure(3);
    clf;
    imagesc(T, F, 10*log10(P)); % frequency-time Plot of the signal
    colorbar;
    % f_min = 17500; %16000;  %
    % f_max = 18500; %22000;  %
    f_min = 0;
    f_max = 10000;
    ylim([f_min f_max]);
    xlabel('Time (s)');
    ylabel('Power/Frequency (dB/Hz)');
    print(fh, '-dpsc', [output_dir filename '.psd.ps']);
end