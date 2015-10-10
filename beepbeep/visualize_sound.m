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
%%
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visualize_sound(filename)
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
    input_dir  = './';
    output_dir = '';

    freq = 17000;


    %% --------------------
    %% Variable
    %% --------------------
    fig_idx = 0;


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, filename = 'tmp'; end
    if nargin < 1, arg = 1; end


    %% --------------------
    %% Main starts
    %% --------------------

    %% --------------------
    %% Read Audio
    %% --------------------
    if DEBUG2, fprintf('Read Audio\n'); end

    find_file = 0;
    while(find_file < 2)
        if find_file == 0
            full_filename = [input_dir filename '.wav'];
        elseif find_file == 1
            full_filename = [input_dir filename '.aac'];
        end
        find_file = find_file + 1;
        
        if exist(full_filename, 'file') == 2
            break;
        end
    end
    if find_file > 2
        error(['  cannot find file: ' full_filename]);
    end
    
    [data, Fs] = audioread(full_filename);
    Ts = 1/Fs;
    nbits = 16;
    fprintf('  file = %s\n', [input_dir filename '.wav']);
    fprintf('  size = %dx%d\n', size(data));


    %% --------------------
    %% High Pass Filter
    %% --------------------
    if DEBUG2, fprintf('High Pass Filter\n'); end

    [B,A] = butter(10, (freq-2000)/(Fs/2), 'high');
    data_filter = filter(B, A, data);
    

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(2,1,1);
    plot(data);

    subplot(2,1,2);
    plot(data_filter);



    %% --------------------
    %% Short-Time Fourier Transform
    %% --------------------
    if DEBUG2, fprintf('Short-time Fourier transform\n'); end

    % window = floor(Fs/2);
    window = floor(Fs/10);
    noverlap = floor(window/4); % 75% overlap
    Nfft = Fs;
    
    [S,F,T,P] = spectrogram(data, window, noverlap, Nfft, Fs);
    [Sf,Ff,Tf,Pf] = spectrogram(data_filter, window, noverlap, Nfft, Fs);


    % ====================================
    % Plotting frequency-time Plot
    % ====================================
    if DEBUG2, fprintf('Plotting frequency-time Plot\n'); end;
    
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    
    imagesc(T, F, real(S));
    colorbar;
    % ylim([15000 18000]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    % title('Time-Frequency plot of a Audio signal');
    % print(fh, '-dpsc', [output_dir filename '.freq-time.ps']);


    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    
    imagesc(Tf, Ff, real(Sf));
    colorbar;
    % ylim([15000 18000]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    % title('Time-Frequency plot of a Audio signal');
    % print(fh, '-dpsc', [output_dir filename '.freq-time.ps']);
end