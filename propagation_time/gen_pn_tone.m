%% gen_pure_tone
function gen_pn_tone()
    
    output_dir = './tx_sound/';
    freq = 8000;
    total_len = 100;
    frame_len = 1;
    amp = 1;
    Fs = 44100;
    Ts = 1/Fs;

    % pn_len = 511;
    % pn_len = 1023;
    pn_len = 2047;
    % pn_len = 4095;
    % pn_len = 8191;

    %% ===================================
    %% generate base sound
    %%   PN sequence
    fprintf('generate base sound\n');
    
    %% --------------
    %% Method 1
    if pn_len == 511
        fbconnection = [0 0 0 0 0 1 0 0 1];
    elseif pn_len == 1023
        fbconnection = [0 0 0 0 0 0 1 0 0 1];
    elseif pn_len == 2047
        fbconnection = [0 0 0 0 0 0 0 1 0 0 1];
    elseif pn_len == 4095
        fbconnection = [0 0 0 0 0 0 0 0 1 0 0 1];
    elseif pn_len == 8191
        fbconnection = [0 0 0 0 0 0 0 0 0 1 0 0 1];
    else
        error('wrong number of bits');
    end
    pnseq = m_sequence(fbconnection)';
    
    %% --------------
    %% Method 2
    % h = commsrc.pn('GenPoly', [[8 2 0]], 'Mask', [1 0 0 0 0 0 1 0]);
    % set(h, 'NumBitsOut', pn_len);
    % pnseq = generate(h);
    % for k=1:length(pnseq)
    %     if(pnseq(k) == 0)
    %         pnseq(k) = -1;
    %     end
    % end
    
    fprintf('  size of PN seq = %dx%d\n', size(pnseq));
    sound_sample = pnseq;
    
    %% cross correlation
    maxpeak = xcorr(pnseq,[pnseq; zeros(1000,1); pnseq]);
    maxpeak = maxpeak(1:floor(end/2)+1);
    fh = figure(1); clf;
    plot(maxpeak);
    
    % return
    %% ===================================
    

    % up-convert
    sound_sample = sound_sample .* sin(2*pi*freq*Ts * [1:length(sound_sample)])';
    n_sound = length(sound_sample);
    fprintf('  sound_sample: %d x %d\n', size(sound_sample));


    %% concatenate base sound
    fprintf('embed base sound\n');
    n_frame = floor(total_len / frame_len);
    for fi = 1:n_frame
        fprintf('  frame %d\n', fi);
        std_idx = int32((fi-1) * frame_len * Fs) + 1;
        end_idx = int32((fi-1) * frame_len * Fs) + n_sound;
        signal(std_idx:end_idx, 2) = sound_sample;
    end

    fprintf('write to the file\n');
    audiowrite([output_dir 'pn.' num2str(pn_len) '.' num2str(freq) '.' num2str(frame_len) '.wav'], amp*signal, Fs);

    length(sound_sample)


    %% cross correlation
    maxpeak = xcorr(sound_sample, signal(1:5*Fs,2));
    maxpeak = maxpeak(1:floor(end/2)+1);
    fh = figure(2); clf;
    plot(maxpeak);
    
