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
%%   ofdm_pn_sim_noise('ofdm.18000', 1)
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ofdm_pn_sim_noise(filename, config)
    addpath('../ofdm');
    
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
    tx_dir  = '../../data/ofdm/gen_data/';
    % input_dir  = './gen_data/';
    % output_dir = '';
    fig_dir = './fig/';


    %% --------------------
    %% Variable
    %% --------------------
    fig_idx = 0;
    data_std_idx = 501;
    

    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, filename = 'ofdm.18000'; end
    if nargin < 2, config = 1; end

    [Nfft, Ncp, symbolRate, Ts, Fs, fc, rollOff, nSamp, filterSpan, preambleSeq, sampleInterval] = get_config_param(config);
    sampleDelay = filterSpan*nSamp/2+nSamp*Ncp;
    preamble_filename = [tx_dir 'preambleC' num2str(config) '.mat'];



    %% --------------------
    %% Main starts
    %% --------------------

    %% ====================
    %% Load data
    %% ====================
    if DEBUG2, fprintf('Load data\n'); end
    
    fprintf('  preamble file: %s\n', preamble_filename);
    load(preamble_filename);
    fprintf('  preamble size = %dx%d\n', size(preamble));


    %% ====================
    %% Generate artificial data
    %% ====================
    data = zeros(1, length(preamble) + data_std_idx*2-2);
    data(data_std_idx:data_std_idx+length(preamble)-1) = preamble;


    fig_idx = fig_idx + 1;
    snrs = [30:-5:0];
    for ni = 1:length(snrs)
        snr = snrs(ni);
        noisy_data = awgn(data, snr, 'measured');


        %% ====================
        %% Calculate correlation coefficient
        %% ====================
        for ti = 1:length(noisy_data)-length(preamble)+1
            corr(ti) = abs(noisy_data(ti:ti+length(preamble)-1) * preamble');
        end


        %% ====================
        %% Visualize
        %% ====================
        fh = figure(fig_idx); clf;

        subplot(2,1,1);
        plot(abs(noisy_data), '-r');
        hold on;
        plot(abs(data), '-b');
        plot(abs(noisy_data-data), '-g');
        set(gca, 'XLim', [1 length(noisy_data)]);
        legend('mix', 'orig', 'noise');

        subplot(2,1,2);
        plot(corr)
        set(gca, 'XLim', [1 length(noisy_data)]);

        title(['noise=' num2str(snr)]);
        print(fh, '-dpsc', [fig_dir filename '.noise' num2str(snr) '.eps']);

        % waitforbuttonpress;
    end


    fig_idx = fig_idx + 1;
    % mp_power_ratios = [0:0.05:0.5];
    % mp_power_ratios = [0:0.1:0.5];
    mp_power_ratios = [0.3];
    mp_shift = -50;
    for si = 1:length(mp_power_ratios)
        mp_power_ratio = mp_power_ratios(si);

        mp = zeros(1,length(data));
        mp(data_std_idx+mp_shift:data_std_idx+mp_shift+length(preamble)-1) = preamble * mp_power_ratio;
        mp_data = data + mp;
        mp_data = mp_data * 0.8;
        data = data * 0.8;
        mp = mp * 0.8;


        %% ====================
        %% Calculate correlation coefficient
        %% ====================
        for ti = 1:length(mp_data)-length(preamble)+1
            corr(ti) = abs(mp_data(ti:ti+length(preamble)-1) * preamble');
        end


        %% ====================
        %% Visualize
        %% ====================
        fh = figure(fig_idx); clf;

        subplot(2,1,1);
        plot(abs(mp_data), '-r.');
        hold on;
        plot(abs(data), '-b.');
        plot(abs(mp), '-g.');
        % set(gca, 'XLim', [1 length(mp_data)]);
        set(gca, 'XLim', [200 1200]);
        legend('mix', 'orig', 'multi-path');

        subplot(2,1,2);
        plot(corr, '-b.');
        % set(gca, 'XLim', [1 length(mp_data)]);
        set(gca, 'XLim', [200 1200]);

        title(['mp power ratio=' num2str(mp_power_ratio)]);

        print(fh, '-dpsc', [fig_dir filename '.mp' num2str(mp_power_ratio) '.eps']);

        % waitforbuttonpress;
    end
    
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
        preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
        % sampleInterval=0.16;
        sampleInterval = 1/1*2 + 4263/Fs;
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
