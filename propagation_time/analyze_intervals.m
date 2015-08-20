%% ====================================
%% Yi-Chao@UT Austin
%%
%% e.g.
%%   analyze_intervals('sinc.8000.1s.10-30-50cm', './data/')
%%   analyze_intervals('sinc.8000.1s.100-1500cm.exp1', './data/')
%% ====================================

%% analyze_intervals: function description
function analyze_intervals(filename, input_dir)
    DEBUG0 = 0;
    DEBUG1 = 1;
    DEBUG2 = 0;  %% progress
    DEBUG3 = 0;  %% basic info
    DEBUG4 = 1;  %% process info
    DEBUG5 = 0;  %% final output
    DEBUG6 = 0;  %% show frequency shift

    % find_peak_method = 'first';
    % find_peak_method = 'max';
    find_peak_method = 'first2';
    % find_peak_method = 'max2';

    peak_thresh = 0.004;
    peak_thresh2 = 50;
    % sound_speed = 35280;
    sound_speed = 33100;
    half_band = 350;

    font_size = 16;

    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, filename = '1k.0.1-1s.10-30cm'; end
    if nargin < 2, input_dir = './rx_sound/'; end

    output_dir = './fig/';
    
    % expression = ['(?<freq>\d+).(?<sound_len>\d+\.*\d*)-(?<sound_itv>\d+\.*\d*)s.(?<dis1>\d+)-(?<dis2>\d+)-(?<dis3>\d+)cm'];
    % tokenNames = regexp(filename, expression, 'names');
    % if(length(tokenNames) > 0)
    %     % freq = str2num(tokenNames(1).freq) * 1000;
    %     freq = str2num(tokenNames(1).freq);
    %     sound_len = str2num(tokenNames(1).sound_len);
    %     sound_itv = str2num(tokenNames(1).sound_itv);
    %     dist1 = str2num(tokenNames(1).dis1);
    %     dist2 = str2num(tokenNames(1).dis2);
    %     dist3 = str2num(tokenNames(1).dis3);
    %     frame_len = sound_len + sound_itv;
    % else
    %     expression = ['(?<freq>\d+).(?<sound_len>\d+\.*\d*)-(?<sound_itv>\d+\.*\d*)s.(?<dis1>\d+)cm'];
    %     tokenNames = regexp(filename, expression, 'names');
    %     freq = str2num(tokenNames(1).freq);
    %     sound_len = str2num(tokenNames(1).sound_len);
    %     sound_itv = str2num(tokenNames(1).sound_itv);
    %     dist1 = str2num(tokenNames(1).dis1);
    %     frame_len = sound_len + sound_itv;
    % end
    expression = ['(snic.)*(?<freq>\d+).(?<sound_len>\d+\.*\d*)-*(?<sound_itv>\d+\.*\d*)*s.(?<dis1>\d+)-*(?<dis2>\d+)*-*(?<dis3>\d+)*cm.*'];
    tokenNames = regexp(filename, expression, 'names');
    if(length(tokenNames) > 0)
        % freq = str2num(tokenNames(1).freq) * 1000;
        freq = str2num(tokenNames(1).freq);
        sound_len = str2num(tokenNames(1).sound_len);
        sound_itv = str2num(tokenNames(1).sound_itv);
        dist1 = str2num(tokenNames(1).dis1);
        dist2 = str2num(tokenNames(1).dis2);
        dist3 = str2num(tokenNames(1).dis3);

        if(sound_itv)
            frame_len = sound_len + sound_itv;
        else
            frame_len = sound_len;
        end

        fprintf('  freq = %d\n', freq);
        fprintf('  sound len = %f\n', sound_len);
        fprintf('  sound itv = %f\n', sound_itv);
        fprintf('  dist1 = %d\n', dist1);
        fprintf('  frame length = %f\n', frame_len);

    else
        error(['wong format: ' filename]);
    end

    
    frame_len = calibrate_interval(filename, frame_len);


    %% ====================================
    %% Read the audio file
    %% ====================================
    if DEBUG2, fprintf('Read file\n'); end

    file_path_name = [input_dir filename '.aac'];
    fprintf('  file: %s\n', file_path_name);
    [wav_data, Fs] = audioread(file_path_name);
    Ts = 1/Fs;
    nbits = 16;

    
    %% ====================================
    %% Filter / Down-convertion
    %% ====================================
    if DEBUG2, fprintf('Down-convertion\n'); end

    f0 = freq;
    f_ctr = f0;
    f_min = f_ctr - half_band;
    f_max = f_ctr + half_band;
    
    if(1)
        wav_data = wav_data .* sin(2*pi*freq*Ts * [1:length(wav_data)]');
        wav_data = lowPassFilterByFFT(wav_data', Fs, 3000, 10)';
        % [b,a] = butter(6, 100/Fs);
        % wav_data = filter(b, a, wav_data);
        wav_data = abs(wav_data);
    else
        [B,A] = butter(2, [f_min/Fs f_max/Fs]);
        wav_data = filter(B, A, wav_data);
    end

    %% --------------
    %% skip the first second
    wav_data = wav_data(5*Fs:end);
    %% --------------

    wav_len = length(wav_data);
    wav_time = 0:1/Fs:(wav_len-1)/Fs;
    if DEBUG4, 
        fprintf('- wav_data %d x %d (Fs=%d, nbits=%d)\n', size(wav_data), Fs, nbits); 
        fprintf('  duration = %fs\n', wav_time(end));
    end
    

    %% ====================================
    %% plot wav in time domain
    %% ====================================
    fh = figure(1);
    clf;
    plot(wav_time, wav_data);
    ylabel('Amplitude');
    xlabel('Time (s)');
    % print(fh, '-dpsc', [output_dir filename '.time.ps']);
    


    %% ====================================
    %% Undersampling
    %% ====================================
    if DEBUG2, fprintf('Undersampling\n'); end
    
    scale = 1;
    f_bases = Fs / scale * [0:scale];
    tmp = sort(f_ctr - f_bases);
    idx = find(tmp > 0);
    f_ctr = tmp(idx(1));
    f_min = f_ctr - half_band;
    f_max = f_ctr + half_band;
    fprintf('  f_min=%f, f_max=%f, f_ctr=%f\n', f_min, f_max, f_ctr);
    
    wav_data = wav_data(1:scale:end);
    Fs = Fs / scale;


    %% ====================================
    %% Find intervals in time domain
    %% ====================================
    frame_samples = int32(floor(frame_len * Fs));
    % opt = ['method=''max'',thresh=' num2str(peak_thresh)];
    opt = ['method=''' find_peak_method ''',thresh=' num2str(peak_thresh)];
    peak_idx = find_peaks(wav_data(:,1), frame_samples, opt);

    %% --------
    %% skip the first k-1 peak
    k = 2;
    peak_idx = peak_idx(k:end);

    dlmwrite([output_dir filename '.td_peaks.txt'], [[1:length(peak_idx)]', wav_time(peak_idx)'], 'delimiter', '\t','precision', 10);
    

    peak_time = wav_time(peak_idx);
    T0 = peak_time(1) - dist1 / sound_speed;
    Ts = T0 + [0:length(peak_time)-1] * frame_len;
    dists = (peak_time - Ts) * sound_speed;
    fprintf('dist=%.2f\n', dists);
    % return;


    %% ====================================
    %% Calculate estimation error
    %% ====================================
    gd = [0:10:100 80:-20:0];
    gds = zeros(1, length(wav_data));
    for di = 1:length(gd)
        std_idx = (5+(di-1)*10) * Fs;
        end_idx = (15+(di-1)*10) * Fs - 1;
        gds(std_idx:min(end,end_idx)) = gd(di);
    end

    err = abs(dists - gds(peak_idx));
    idx = find(err < 5);
    if length(idx) > 0
        err = err(idx);
        % err = mean(err);
        fprintf('  > %d/%d samples: avg=%.2fcm, max=%.2f, std=%.5f\n', length(idx), length(peak_idx), mean(err), max(err), std(err));
    end

    %% ====================================
    %% plot wav and peaks in time domain
    %% ====================================
    fh = figure(1);
    clf;
    subplot(2,1,1);
    plot(wav_time, wav_data);
    hold on;
    plot(wav_time(peak_idx), wav_data(peak_idx), 'or');
    set(gca, 'XLim', [0 wav_time(end)]);
    ylabel('Amplitude', 'FontSize', font_size);
    xlabel('Time (s)', 'FontSize', font_size);
    set(gca, 'FontSize', font_size);
    grid on;
    % print(fh, '-dpsc', [output_dir filename '.time.ps']);

    subplot(2,1,2);
    lh = plot(peak_time, dists, '-b.');
    set(lh, 'LineWidth', 1);
    % set(lh, 'MarkerSize', 10);
    hold on;
    % plot([0  15], [10 10], '-r');
    % plot([20 35], [30 30], '-r');
    % plot([40 60], [50 50], '-r');
    % plot([0  15], [100 100], '-r');
    % plot([20 35], [300 300], '-r');
    % plot([40 60], [500 500], '-r');
    % gd = [0:10:100 80:-20:0];
    % for di = 1:length(gd)
    %     plot([5 15]+(di-1)*10, [gd(di) gd(di)], '-r');
    % end
    % plot(wav_time, gds, '-r');
    set(gca, 'XLim', [0 wav_time(end)]);
    set(gca, 'YLim', [0 3000]);
    ylabel('distance (cm)', 'FontSize', font_size);
    % set(gca, 'YTick', [-10:10:max(dists)+10]);
    xlabel('Time (s)', 'FontSize', font_size);
    set(gca, 'FontSize', font_size);
    grid on;
    
    print(fh, '-dpsc', [output_dir filename '.time_dists.ps']);
    return
    

    % %% ====================================
    % %% Short-Time Fourier Transform
    % %% ====================================
    % if DEBUG2, fprintf('Short-time Fourier transform\n'); end

    % % window = Fs/2; % Should be minimum twice the maximum frequency we want to analyze
    % % window = floor(Fs/32);
    % window = floor(Fs/128);
    % noverlap = floor(window/4); % 75% overlap
    % Nfft = Fs;
    % % Nfft = window;
    
    % % Spectrogram takes the STFT of the signal
    % % P matrix contains the power spectral density of each segment of the STFT
    % [S,F,T,P] = spectrogram(wav_data, window, noverlap, Nfft, Fs);

    % %% ====================================
    % %% find mag of target freq
    % %% ====================================
    % fh = figure(5); clf;
    % f_ind = freq2ind(freq, Fs/2, length(F));
    % power_f = 10*log10(P(f_ind, :));


    % %% ====================================
    % %% Find intervals in freq domain
    % %% ====================================
    % si = 1;
    % peak_idx = [];
    % frame_samples = int32(floor(length(power_f) / T(end) * frame_len));

    % while(si < length(power_f)) 
    %     wav_seg = power_f(si:min(end,si+frame_samples-1));

    %     % max(wav_seg) - min(wav_seg)
    %     if max(wav_seg) - min(wav_seg) < peak_thresh2
    %         si = si + frame_samples/2;
    %         continue;
    %     end
        
    %     % [val, idx] = max(wav_seg);
    %     idx = find(wav_seg > -120);
    %     idx = idx(1);
    %     peak_idx = [peak_idx, si + idx - 1];

    %     si = si + idx - 1 + frame_samples/2;
    % end

    % dlmwrite([output_dir filename '.fd_peaks.txt'], [[1:length(peak_idx)]', T(peak_idx)'], 'delimiter', '\t','precision', 10);

    % peak_time = T(peak_idx);
    % T0 = peak_time(1) - dist1 / sound_speed;
    % Ts = T0 + [0:length(peak_time)-1] * frame_len;
    % dists = (peak_time - Ts) * sound_speed;
    % fprintf('--------------\n');
    % fprintf('dist=%.2f\n', dists);

    % plot(T, power_f);
    % hold on;
    % plot(T(peak_idx), power_f(peak_idx), 'or');
    % ylabel('Amplitude');
    % xlabel('Time (s)');
    % print(fh, '-dpsc', [output_dir filename '.time_freq.ps']);
    
    % size(P)
    % return


    %% ====================================
    %% Plotting frequency-time Plot
    %% ====================================
    % if DEBUG4, fprintf('- plot S: %d x %d\n', size(S)); end;
    % fh = figure(2);
    % clf;
    % imagesc(T, F, real(S));
    % colorbar;
    % ylim([17500 18500]);
    % xlabel('Time (s)');
    % ylabel('Frequency (Hz)');
    % title('Time-Frequency plot of a Audio signal');
    % print(fh, '-dpsc', [output_dir filename '.freq-time.ps']);

    % %% ====================================
    % %% Plotting Power Spectral Density
    % %% ====================================
    % if DEBUG4, fprintf('- plot P: %d x %d\n', size(P)); end;
    % fh = figure(3);
    % clf;
    % imagesc(T, F, 10*log10(P)); % frequency-time Plot of the signal
    % colorbar;
    % % f_min = 17500; %16000;  %
    % % f_max = 18500; %22000;  %
    % % f_min = 3200;
    % % f_max = 3400;
    % ylim([f_min f_max]);
    % xlabel('Time (s)');
    % ylabel('Power/Frequency (dB/Hz)');
    % print(fh, '-dpsc', [output_dir filename '.f' num2str(f_ctr) '.psd.ps']);

    % return;

    % %% ====================================
    % %% Calculate freq shift
    % %% ====================================
    % if DEBUG2, fprintf('Highest peak:\n'); end;
    % % f_min = 17900;  f_min_ind = freq2ind(f_min, Fs/2, length(F));
    % % f_max = 18100;  f_max_ind = freq2ind(f_max, Fs/2, length(F));
    % % f_min = 17000;  f_min_ind = freq2ind(f_min, Fs/2, length(F));
    % % f_max = 21000;  f_max_ind = freq2ind(f_max, Fs/2, length(F));
    % f_min_ind = freq2ind(f_min, Fs/2, length(F));
    % f_max_ind = freq2ind(f_max, Fs/2, length(F));    
    % target_band = 10*log10(P(f_min_ind:f_max_ind, :));
    % % fh = figure(4); clf; imagesc(tmp); print(fh, '-dpsc', ['tmp/tmp.ps']);
    % target_freq = F(f_min_ind:f_max_ind);
    % target_f     = f_ctr;
    % target_f_min = f_ctr-0;
    % target_f_max = f_ctr+0;
    % % target_f_min = 18000;
    % % target_f_max = 20000;
    % % target_f     = 19000;
    % % target_f_min = 900;
    % % target_f_max = 900;
    % % target_f     = 900;
    % % plot_t = [6, 15, 18];
    % % plot_t = [7, 9, 11];
    % % plot_t = [5, 7, 9];
    % % plot_t = [56];
    % % plot_t = [11, 12, 14, 16, 17, 18];
    % plot_t = [1];

    % peak_freq = zeros(1, size(target_band, 2));
    % peak_freq_power = zeros(1, size(target_band, 2));
    % velocities = zeros(1, size(target_band, 2));
    
    
    % for t = 1:size(target_band, 2)
    %     if DEBUG6, fprintf('- t%d(%f):\n', t, T(t)); end;

    %     method_idx = 0;
    %     ts = target_band(:,t)';

    %     if isinf(ts)
    %         % idx = find(ts, inf)
    %         continue;
    %     end
        

    %     %% ====================================
    %     %% plot Power over frequency of this time period
    %     %% ====================================
    %     if ismember(t, plot_t), 
    %         % cps = detect_change_points(ts, num_bootstrap, conf_threshold, rank_method, filter_method, 'no_plot');
    %         fh = figure(4); clf; 
    %         plot(target_freq, ts); % hold on;
    %         % plot(target_freq(cps), ts(cps), 'or');
    %         xlim([target_freq(1), target_freq(end)]);
    %         xlabel('Frequency (Hz)'); ylabel('Power/Frequency (db/Hz)');
    %         print(fh, '-dpsc', [output_dir  filename '.t' int2str(t) '.ps']); 
    %     end;

        
    %     [spike_idx] = find_spikes(ts);
    %     peak_val = ts(spike_idx);
    %     [sort_val, sort_ind] = sort(peak_val, 'descend');
        
    %     first_val = sort_val(1);
    %     first_idx = spike_idx(sort_ind(1));
    %     % [vr, vs] = doppler_velocity(target_freq(first_idx), target_f, -1, 0);
    %     [vr, vs] = doppler_velocity(f0, target_freq(first_idx)-target_f, -1, 0);

    %     method_idx = method_idx + 1;
    %     peak_freq(method_idx, t) = target_freq(first_idx);
    %     peak_freq_power(method_idx, t) = first_val;
    %     velocities(method_idx, t) = vr;

    %     if DEBUG6, fprintf('  1st peak: freq=%f, pwer=%f, v=%f\n', target_freq(first_idx), first_val, vr); end;

    % end


    % %% ----------------------
    % %% Filter: anomalies in velocity
    % %% ----------------------
    % % v_difs = abs(velocities(:,2:end) - velocities(:,1:end-1));
    % % avg_v_dif = mean2(v_difs(:));
    % % std_v_dif = std(v_difs(:));
    % % method_idx = 1;
    % % fprintf('  v diff avg = %f, std = %f\n', avg_v_dif, std_v_dif);
    % % for t = 2:size(velocities, 2)
    % %     fprintf('  t%d: v=%f, v diff=%f\n', t, velocities(method_idx, t), abs(velocities(method_idx, t) - velocities(method_idx, t-1)));

    % %     if abs(velocities(method_idx, t) - velocities(method_idx, t-1)) > (avg_v_dif + 2*std_v_dif)
    % %         fprintf('  >> too large\n');
    % %         velocities(method_idx, t) = velocities(method_idx, t-1);
    % %     end
    % % end


    % %% velocity
    % avg_v = cumsum(velocities, 2) ./ (repmat([1:size(velocities,2)], size(velocities,1), 1));
    % traces = avg_v .* repmat(T, size(velocities,1), 1);


    % %% plot the center frequency over time
    % for mi = 1:size(peak_freq, 1)
    %     fprintf('md%d: move distance = %f\n', mi, traces(mi, end));

    %     this_filename = ['./tmp/' filename '.freq.method' int2str(mi) '.f' num2str(f_ctr)];
    %     plot_center_freq(T, peak_freq(mi, :), peak_freq_power(mi, :), traces(mi, :), this_filename);
    % end

end


%% freq2ind
function [ind] = freq2ind(freq, max_freq, F_len)
    ind = floor(freq / max_freq * F_len);
end


%% find_spikes: function description
function [spike_idx] = find_spikes(ts)
    phase = 0;
    spike_idx = [-1];
    
    for ti = 1:length(ts)
        val = ts(ti);
        if phase == 0,
            phase = 1;
        elseif phase == 1,
            if prev_val > val,
                if spike_idx(1) == -1,
                    spike_idx = [ti-1];
                else
                    spike_idx = [spike_idx, ti-1];
                end
                
                phase = 2;
            end
        elseif phase == 2,
            if prev_val < val,
                phase = 1;
            end
        end

        prev_val = val;
    end

    if spike_idx(1) == -1,
        spike_idx = [1];
    end
end


function [itvl] = calibrate_interval(filename, itvl)
    method = 'MANUAL';
    % method = 'AUTO';

    if strcmp(method, 'MANUAL')
        itvl = calibrate_interval_manual(filename, itvl);
    elseif strcmp(method, 'AUTO')
        itvl = calibrate_interval_auto(filename, itvl);
    end     
end

function [itvl] = calibrate_interval_manual(filename, itvl)
    % frame_len = frame_len + 1.43/44100;  %% 8000.0.01-0.5s.10-30-50cm
    % frame_len = frame_len + 2.8/44100; %% sinc.8000.1s.10-30-50cm
    % frame_len = frame_len + 0.27/44100; %% sinc.8000.0.1s.10-30-50cm.aac
    % frame_len = frame_len + 2.6/44100; %% sinc.8000.1s.100-300-500cm.aac
    % frame_len = frame_len + 2/44100; %% sinc.8000.1s.100-1500cm.exp1.aac
    % frame_len = frame_len + 3/44100; %% sinc.8000.1s.100-1500cm.exp3.aac
    if strcmp(filename, '8000.0.01-0.5s.10-30-50cm')
        itvl = itvl + 1.43/44100;
    elseif strcmp(filename, 'sinc.8000.1s.10-30-50cm')
        itvl = itvl + 2.8/44100;
    elseif strcmp(filename, 'sinc.8000.0.1s.10-30-50cm')
        itvl = itvl + 0.27/44100;
    elseif strcmp(filename, 'sinc.8000.1s.100-300-500cm')
        itvl = itvl + 2.6/44100;
    elseif strcmp(filename, 'sinc.8000.1s.100-1500cm.exp1')
        itvl = itvl + 2.6/44100;
    elseif strcmp(filename, 'sinc.8000.1s.100-1500cm.exp2')
        itvl = itvl + 2.7/44100;
    elseif strcmp(filename, 'sinc.8000.1s.100-1500cm.exp3')
        itvl = itvl + 2.8/44100;
    else
        itvl = itvl + 2.6/44100;
    end
end


function [itvl] = calibrate_interval_auto(filename, itvl)

end