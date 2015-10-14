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
%%   ofdm_rx_itvl('ofdm.18000', 1)
%%   ofdm_rx_itvl('0914.exp1', 1)
%%   ofdm_rx_itvl('0914.exp2', 1)
%%   ofdm_rx_itvl('1014.exp1', 1)
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ofdm_rx_itvl(filename, config)
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
    tx_dir  = './gen_data/';
    input_dir  = './rx_sound/';
    % input_dir  = './gen_data/';
    % output_dir = '';
    fig_dir = './fig/';


    %% --------------------
    %% Variable
    %% --------------------
    fig_idx = 0;
    

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
    audio_filename = [input_dir filename '.wav'];
    if exist(audio_filename, 'file') ~= 2,
        audio_filename = [input_dir filename '.aac'];
    end
    fprintf('  audio file: %s\n', audio_filename);
    [analogData,~] = audioread(audio_filename);
    analogData = analogData.';

    fprintf('  preamble file: %s\n', preamble_filename);
    load(preamble_filename);
    fprintf('  preamble size = %dx%d\n', size(preamble));


    %% ----------------
    %% XXX: need to modify
    if size(analogData,1) > 1
        analogData = analogData(2,1:10*Fs);
    end
    analogData = analogData(:,(0*Fs+1):min(end,100*Fs));
    %% ----------------
    fprintf('  data size = %dx%d\n', size(analogData));

    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % plot(analogData);
    % hold on;
    % plot(courseStartIndex, (analogData(courseStartIndex)), 'ro');
    % title('raw audio');


    %% ====================
    %% downconvert
    %% ====================
    if DEBUG2, fprintf('Downconvert\n'); end

    T = numel(analogData);
    analogData = analogData .* exp(-1i*2*pi*fc*(1:T)*Ts);
    % analogDataFFT=fft(analogData);
    % plot(abs(analogDataFFT));

    if DEBUG2, fprintf('Low Pass Filter\n'); end
    analogData=lowPassFilterByFFT(analogData,Fs,2000,500);

    % if DEBUG2, fprintf('FFT\n'); end
    % analogDataFFT=fft(analogData);

    % Np=numel(preamble);
    % for i = 1:length(analogData)
    %     if mod(i-1,length(analogData)/10) == 0
    %         fprintf('%d\n', (i-1)/length(analogData)*100);
    %     end
    %     if i+Np-1 > length(analogData)
    %         break;
    %     end
    %     corr(i) = abs(analogData(i:i+Np-1)*preamble');
    % end

    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % subplot(2,1,1);
    % plot(abs(analogData));
    % subplot(2,1,2);
    % plot(corr);
    % return


    % course timing synchornizing
    if DEBUG2, fprintf('course timing synchornizing\n'); end
    windowSize=300;
    detectLength=sampleInterval/Ts*2;
    courseStartIndex=findStartIndexByDoubleWin(analogData,windowSize,detectLength);
    
    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % plot(abs(analogData));
    % hold on;
    % plot(courseStartIndex, abs(analogData(courseStartIndex)), 'ro');
    % print(fh, '-dpsc', [fig_dir filename '.eps']);


    %% ====================
    %% split the received signals
    %% ====================
    if DEBUG2, fprintf('split the received signals\n'); end

    dataindex=courseStartIndex;
    dataSpan=(Nfft+Ncp-1)*nSamp+nSamp*filterSpan+1+Fs/100*2;
    fprintf('  dataSpan=%d\n', dataSpan);
    dataInterval = sampleInterval/Ts-dataSpan;
    fprintf('  dataInterval=%d\n', dataInterval);
    Ns=400;
    analogData_orig = analogData;
    [analogData, split_idx]=split(analogData,dataindex,dataSpan,dataInterval,Ns);
    Ns = size(analogData, 1);
    fprintf('  data after split = %dx%d\n', size(analogData));
    fprintf('  #split = %dx%d\n', size(split_idx));

    Nd=numel(analogData(1,:));
    Np=numel(preamble);
    fprintf('  #data bits=%d, #preamble bits=%d\n', Nd, Np);
    h=zeros(Ns,Nfft);
    startIndex=zeros(1,Ns);
    startIndexCorr=zeros(1,Ns);

    fig_idx = fig_idx + 1;
    for j=1:Ns
        %% ====================
        %% refined timing sychronization  
        %% ====================
        if DEBUG2, fprintf('refined timing sychronization: %d\n', j); end
        
        corr=zeros(1,Nd-Np+1);
        maxCorr=0;
        maxCorrIndex=0;
        for i=0:Nd-Np
            corr(i+1)=abs(analogData(j,i+1:i+Np)*preamble');
            if corr(i+1)>maxCorr
                maxCorr=corr(i+1);
                maxCorrIndex=i+1;
            end
        end
        % disp(maxCorrIndex);

        % fh = figure(fig_idx); clf;
        % plot(corr);
        % waitforbuttonpress;

        fh = figure(fig_idx); clf;
        subplot(2,1,1);
        plot(abs(analogData(j,:)));
        hold on;
        plot(maxCorrIndex, abs(analogData(j,maxCorrIndex)), 'ro');
        set(gca, 'XLim', [1 size(analogData,2)]);
        subplot(2,1,2);
        plot(corr);
        hold on;
        plot(maxCorrIndex, corr(maxCorrIndex), 'ro');
        set(gca, 'XLim', [1 size(analogData,2)]);
        waitforbuttonpress;
        print(fh, '-dpsc', [fig_dir filename '.' num2str(j) '.eps']);
        
        startIndex(j)=maxCorrIndex+sampleDelay;
        startIndexCorr(j)=0;
        minh=100;
        data_orgn=2*preambleSeq-1;
        data_orgn=data_orgn(1:Nfft);
        for i=startIndex(j)-10:startIndex(j)+10
            % downsampling
            data_fft=analogData(j,i:nSamp:i+nSamp*(Nfft-1));
            
            % convert to frequency domain
            data=fft(data_fft);
            
            % calculate channel gain        
            Hf=data./data_orgn;
            htemp=ifft(Hf);
            if abs(htemp(Nfft))<minh
                minh=abs(htemp(Nfft));
                ht=htemp;
                startIndexCorr(j)=i;
            end
        end
        % disp(startIndexCorr(j)-sampleDelay);
        % stem(abs(ht));
        h(j,:)=ht;

        split_idx(j) = split_idx(j) + startIndexCorr(j) - sampleDelay;
    end


    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;
    plot(abs(analogData_orig));
    hold on;
    plot(split_idx, abs(analogData_orig(split_idx)), 'ro');
    print(fh, '-dpsc', [fig_dir filename '.eps']);

    split_idx(2:end) - split_idx(1:end-1)
    
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
