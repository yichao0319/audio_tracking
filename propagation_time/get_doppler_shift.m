%% ====================================
%% Yi-Chao@UT Austin
%%
%% e.g.
%%   get_doppler_shift('./rx_sound/', '07.03.16000.40-10', 16000)
%% ====================================

%% get_doppler_shift: function description
function [T, velocities, traces] = get_doppler_shift(input_dir, filename, f0)

    DEBUG0 = 0;
    DEBUG1 = 1;
    DEBUG2 = 0;  %% progress
    DEBUG3 = 0;  %% basic info
    DEBUG4 = 0;  %% process info
    DEBUG5 = 0;  %% final output
    DEBUG6 = 0;  %% show frequency shift

    PLOT_PSD = 0;
    PLOT_FREQ_SHIFT = 0;


    %% parameters for change point detection
    num_bootstrap = 100;   % number of iterations for Bootstrap Analysis
    conf_threshold = 0.99;  % the threshold for Signiﬁcance Testing
    half_band = 200;
    % half_band = 1000;
    rank_method = 0;
    filter_method = 0;


    %% ====================================
    %% Read the .wav file
    %% ====================================
    if DEBUG2, fprintf('Read file\n'); end

    file_path_name = [input_dir filename '.aac'];
    % [wav_data, Fs, nbits] = wavread(file_path_name);
    [wav_data, Fs] = audioread(file_path_name);
    nbits = 16;
    
    %% --------------
    %% skip the first second
    wav_data = wav_data(Fs:end);
    %% --------------

    wav_len = length(wav_data);
    wav_time = 0:1/Fs:(wav_len-1)/Fs;
    if DEBUG4, 
        fprintf('- wav_data %d x %d (Fs=%d, nbits=%d)\n', size(wav_data), Fs, nbits); 
        fprintf('  duration = %f\n', wav_time(end));
    end

    %% ====================================
    %% plot wav in time domain
    %% ====================================
    % fh = figure(1);
    % clf;
    % plot(wav_time, wav_data);
    % ylabel('Amplitude');
    % xlabel('Time (s)');
    % print(fh, '-dpsc', ['tmp/' filename '.time.ps']);

    
    %% ====================================
    %% Filtering
    %% ====================================
    % if DEBUG2, fprintf('Filter\n'); end
    
    % f0 = 18000;
    % f0 = 17999.5;
    % f0 = 17999;
    f_ctr = f0;
    f_min = f_ctr - half_band;
    f_max = f_ctr + half_band;
    [B,A] = butter(2, [f_min/Fs f_max/Fs]);
    wav_data = filter(B, A, wav_data);


    %% ====================================
    %% Undersampling
    %% ====================================
    % if DEBUG2, fprintf('Undersampling\n'); end
    
    % scale = 1;
    % f_bases = Fs / scale * [0:scale];
    % tmp = sort(f_ctr - f_bases);
    % idx = find(tmp > 0);
    % f_ctr = tmp(idx(1));
    % f_min = f_ctr - half_band;
    % f_max = f_ctr + half_band;
    % fprintf('  f_min=%f, f_max=%f, f_ctr=%f\n', f_min, f_max, f_ctr);
    
    % wav_data = wav_data(1:scale:end);
    % Fs = Fs / scale;

    

    %% ====================================
    %% Short-Time Fourier Transform
    %% ====================================
    if DEBUG2, fprintf('Short-time Fourier transform\n'); end

    % window = Fs/2; % Should be minimum twice the maximum frequency we want to analyze
    window = floor(Fs/32);
    noverlap = floor(window/4); % 75% overlap
    Nfft = Fs;
    % Nfft = window;
    
    % Spectrogram takes the STFT of the signal
    % P matrix contains the power spectral density of each segment of the STFT
    [S,F,T,P] = spectrogram(wav_data, window, noverlap, Nfft, Fs);


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
    % print(fh, '-dpsc', ['tmp/' filename '.freq-time.ps']);

    %% ====================================
    %% Plotting Power Spectral Density
    %% ====================================
    if PLOT_PSD
        if DEBUG4, fprintf('- plot P: %d x %d\n', size(P)); end;
    
        fh = figure(3);
        clf;
        imagesc(T, F, 10*log10(P)); % frequency-time Plot of the signal
        colorbar;
        % f_min = 17500; %16000;  %
        % f_max = 18500; %22000;  %
        % f_min = 3200;
        % f_max = 3400;
        ylim([f_min f_max]);
        xlabel('Time (s)');
        ylabel('Power/Frequency (dB/Hz)');
        print(fh, '-dpsc', ['fig/' filename '.f' num2str(f_ctr) '.psd.ps']);
    end


    %% ====================================
    %% Calculate freq shift
    %% ====================================
    if DEBUG2, fprintf('Highest peak:\n'); end;
    % f_min = 17900;  f_min_ind = freq2ind(f_min, Fs/2, length(F));
    % f_max = 18100;  f_max_ind = freq2ind(f_max, Fs/2, length(F));
    % f_min = 17000;  f_min_ind = freq2ind(f_min, Fs/2, length(F));
    % f_max = 21000;  f_max_ind = freq2ind(f_max, Fs/2, length(F));
    f_min_ind = freq2ind(f_min, Fs/2, length(F));
    f_max_ind = freq2ind(f_max, Fs/2, length(F));    
    target_band = 10*log10(P(f_min_ind:f_max_ind, :));
    % fh = figure(4); clf; imagesc(tmp); print(fh, '-dpsc', ['tmp/tmp.ps']);
    target_freq = F(f_min_ind:f_max_ind);
    target_f     = f_ctr;
    target_f_min = f_ctr-0;
    target_f_max = f_ctr+0;
    % target_f_min = 18000;
    % target_f_max = 20000;
    % target_f     = 19000;
    % target_f_min = 900;
    % target_f_max = 900;
    % target_f     = 900;
    % plot_t = [6, 15, 18];
    % plot_t = [7, 9, 11];
    % plot_t = [5, 7, 9];
    % plot_t = [56];
    % plot_t = [11, 12, 14, 16, 17, 18];
    plot_t = [1];

    peak_freq = zeros(1, size(target_band, 2));
    peak_freq_power = zeros(1, size(target_band, 2));
    velocities = zeros(1, size(target_band, 2));
    
    
    for t = 1:size(target_band, 2)
        if DEBUG6, fprintf('- t%d(%f):\n', t, T(t)); end;

        method_idx = 0;
        ts = target_band(:,t)';

        if isinf(ts)
            % idx = find(ts, inf)
            continue;
        end
        

        %% ====================================
        %% plot Power over frequency of this time period
        %% ====================================
        % if ismember(t, plot_t), 
        %     % cps = detect_change_points(ts, num_bootstrap, conf_threshold, rank_method, filter_method, 'no_plot');
        %     fh = figure(4); clf; 
        %     plot(target_freq, ts); % hold on;
        %     % plot(target_freq(cps), ts(cps), 'or');
        %     xlim([target_freq(1), target_freq(end)]);
        %     xlabel('Frequency (Hz)'); ylabel('Power/Frequency (db/Hz)');
        %     print(fh, '-dpsc', ['tmp/'  filename '.t' int2str(t) '.ps']); 
        % end;

        
        %% ====================================
        %% 1) First spike
        %% ====================================
        [spike_idx] = find_spikes(ts);
        peak_val = ts(spike_idx);
        [sort_val, sort_ind] = sort(peak_val, 'descend');
        
        first_val = sort_val(1);
        first_idx = spike_idx(sort_ind(1));
        % [vr, vs] = doppler_velocity(target_freq(first_idx), target_f, -1, 0);
        [vr, vs] = doppler_velocity(f0, target_freq(first_idx)-target_f, -1, 0);

        method_idx = method_idx + 1;
        peak_freq(method_idx, t) = target_freq(first_idx);
        peak_freq_power(method_idx, t) = first_val;
        velocities(method_idx, t) = vr;

        if DEBUG6, fprintf('  1st peak: freq=%f, pwer=%f, v=%f\n', target_freq(first_idx), first_val, vr); end;


        %% ====================================
        %% 2) Second spike
        %% ====================================
        % second_val = sort_val(2);
        % second_idx = spike_idx(sort_ind(2));
        % % if t==6, 
        % %     fh = figure(4); clf; 
        % %     plot(target_freq, ts); hold on;
        % %     plot(target_freq(spike_idx), peak_val, 'or'); hold on;
        % %     plot(target_freq(second_idx), second_val, 'xg');
        % %     print(fh, '-dpsc', ['tmp/tmp.ps']); 
        % % end
        % [vr, vs] = doppler_velocity(target_freq(second_idx), target_f, -1, 0);
        % [vr, vs] = doppler_velocity(f0, target_freq(second_idx) - target_f, -1, 0);
        
        % method_idx = method_idx + 1;
        % peak_freq(method_idx, t) = target_freq(second_idx);
        % peak_freq_power(method_idx, t) = second_val;
        % velocities(method_idx, t) = vr;
        
        % if DEBUG6, fprintf('  2nd peak: freq=%f, pwer=%f, v=%f\n', target_freq(second_idx), second_val, vr); end;

        %% ====================================
        %% 3) Third spike
        %% ====================================
        % if(length(sort_val) >= 3)
        %     third_val = sort_val(3);
        %     third_idx = spike_idx(sort_ind(3));
        %     [vr, vs] = doppler_velocity(target_freq(third_idx), target_f, -1, 0);
        %     [vr, vs] = doppler_velocity(f0, target_freq(third_idx) - target_f, -1, 0);

        %     method_idx = method_idx + 1;
        %     peak_freq(method_idx, t) = target_freq(third_idx);
        %     peak_freq_power(method_idx, t) = third_idx;
        %     velocities(method_idx, t) = vr;

        %     if DEBUG6, fprintf('  3rd peak: freq=%f, pwer=%f, v=%f\n', target_freq(third_idx), third_val, vr); end;
        % end


        %% ====================================
        %% 4) Centor frequency of the first spike
        %% ====================================
        % [spk_min_ind, spk_max_ind] = find_spike_range(first_idx, ts);
        
        % first_spk_freq = target_freq([spk_min_ind:spk_max_ind])';
        % first_spk_pwer = ts([spk_min_ind:spk_max_ind]);
        % % [spk_freq, spk_rss] = weighted_freq(first_spk_freq, first_spk_pwer);
        
        % if ismember(t, plot_t), 
        %     fh = figure(5);
        %     clf;
        %     lh1 = plot(target_freq, ts); hold on;
        %     set(lh1, 'LineWidth', 2);
        %     plot(target_freq(first_idx), first_val, 'or'); hold on;
        %     lh2 = plot(first_spk_freq, first_spk_pwer, '.y-');
        %     set(lh2, 'LineWidth', 2);
        %     grid on;
        %     xlim([target_freq(1), target_freq(end)]);
        %     xlabel('Frequency (Hz)'); ylabel('Power/Frequency (db/Hz)');
        %     print(fh, '-dpsc', ['tmp/' filename '.spk' int2str(t) '.ps']); 
        % end
        
        % % [vr, vs] = doppler_velocity(mean(first_spk_freq), target_f, -1, 0);
        % [vr, vs] = doppler_velocity(f0, mean(first_spk_freq) - target_f, -1, 0);

        % method_idx = method_idx + 1;
        % peak_freq(method_idx, t) = mean(first_spk_freq);
        % peak_freq_power(method_idx, t) = mean(first_spk_pwer);
        % velocities(method_idx, t) = vr;

        % if DEBUG6, fprintf('  1st spike center: freq=%f, pwer=%f, v=%f\n', mean(first_spk_freq), mean(first_spk_pwer), vr); end;


        %% ====================================
        %% 5) Centor frequency of the second spike
        %% ====================================
        % [spk_min_ind, spk_max_ind] = find_spike_range(second_idx, ts);
        
        % second_spk_freq = target_freq([spk_min_ind:spk_max_ind])';
        % second_spk_pwer = ts([spk_min_ind:spk_max_ind]);
        
        % [vr, vs] = doppler_velocity(mean(second_spk_freq), target_f, -1, 0);
        % [vr, vs] = doppler_velocity(f0, mean(second_spk_freq) - target_f, -1, 0);

        % method_idx = method_idx + 1;
        % peak_freq(method_idx, t) = mean(second_spk_freq);
        % peak_freq_power(method_idx, t) = mean(second_spk_pwer);
        % velocities(method_idx, t) = vr;

        % if DEBUG6, fprintf('  2nd spike center: freq=%f, pwer=%f, v=%f\n', mean(second_spk_freq), mean(second_spk_pwer), vr); end;


        %% ====================================
        %% 6) Centor frequency of the third spike
        %% ====================================
        % if(length(sort_val) >= 3)
        %     [spk_min_ind, spk_max_ind] = find_spike_range(third_idx, ts);
            
        %     third_spk_freq = target_freq([spk_min_ind:spk_max_ind])';
        %     third_spk_pwer = ts([spk_min_ind:spk_max_ind]);
            
        %     [vr, vs] = doppler_velocity(mean(third_spk_freq), target_f, -1, 0);
        %     [vr, vs] = doppler_velocity(f0, mean(third_spk_freq) - target_f, -1, 0);

        %     method_idx = method_idx + 1;
        %     peak_freq(method_idx, t) = mean(third_spk_freq);
        %     peak_freq_power(method_idx, t) = mean(third_spk_pwer);
        %     velocities(method_idx, t) = vr;

        %     if DEBUG6, fprintf('  3rd spike center: freq=%f, pwer=%f, v=%f\n', mean(third_spk_freq), mean(third_spk_pwer), vr); end;
        % end


        %% ====================================
        %% 7) Centor frequency
        %% ====================================
        % [center_freq, avg_rss] = weighted_freq(target_freq', ts);
        % % [vr, vs] = doppler_velocity(center_freq, target_f, -1, 0);
        % [vr, vs] = doppler_velocity(f0, center_freq - target_f, -1, 0);

        % method_idx = method_idx + 1;
        % peak_freq(method_idx, t) = center_freq;
        % peak_freq_power(method_idx, t) = avg_rss;
        % velocities(method_idx, t) = vr;

        % if DEBUG6, fprintf('  Center: freq=%f, pwer=%f, v=%f\n', center_freq, avg_rss, vr); end;


        %% ====================================
        %% 8) Change Point Detection
        %% ====================================
        % cps = detect_change_points(ts, num_bootstrap, conf_threshold, rank_method, filter_method, 'no_plot');
        % cps_freq = target_freq(cps);
        % freq_shift_thresh = 0;
        % idx = find(cps_freq <= (target_f_min-freq_shift_thresh));
        % if length(idx) > 0
        %     cps_min_ind = cps(idx(end));
        % else
        %     cps_min_ind = cps(1);
        % end
            
        % idx = find(cps_freq >= (target_f_max+freq_shift_thresh));
        % if(length(idx) > 0)
        %     cps_max_ind = cps(idx(1));
        % else
        %     cps_max_ind = cps(end);
        % end
        
        % cps_target_freq = target_freq([cps_min_ind:cps_max_ind])';
        % cps_target_pwer = ts([cps_min_ind:cps_max_ind]);

        % % [center_freq, avg_rss] = weighted_freq(cps_target_freq, cps_target_pwer);
        % [vr, vs] = doppler_velocity(mean(cps_target_freq), target_f, -1, 0);
        % [vr, vs] = doppler_velocity(f0, mean(cps_target_freq) - target_f, -1, 0);

        % method_idx = method_idx + 1;
        % peak_freq(method_idx, t) = mean(cps_target_freq);
        % peak_freq_power(method_idx, t) = mean(cps_target_pwer);
        % velocities(method_idx, t) = vr;

        % if DEBUG6, fprintf('  CPs: freq=%f, pwer=%f, v=%f\n', mean(cps_target_freq), mean(cps_target_pwer), vr); end;
        
        % if ismember(t, plot_t), 
        %     fh = figure(4); clf; 
        %     lh1 = plot(target_freq, ts); hold on;
        %     set(lh1, 'LineWidth', 2);
        %     lh2 = plot(target_freq(cps), ts(cps), 'or'); hold on;
        %     lh3 = plot(cps_target_freq, cps_target_pwer, '.y-.'); hold on;
        %     set(lh3, 'LineWidth', 2);
        %     plot(target_freq([cps_min_ind, cps_max_ind]), ts([cps_min_ind, cps_max_ind]), 'xg');
        %     xlim([target_freq(1), target_freq(end)]);
        %     xlabel('Frequency (Hz)'); ylabel('Power/Frequency (db/Hz)');
        %     grid on;
        %     legend([lh2, lh3], {'change points', 'sound'});
        %     print(fh, '-dpsc', ['tmp/' filename '.cps' int2str(t) '.ps']); 
        % end

    end


    %% ----------------------
    %% Filter: anomalies in velocity
    %% ----------------------
    % v_difs = abs(velocities(:,2:end) - velocities(:,1:end-1));
    % avg_v_dif = mean2(v_difs(:));
    % std_v_dif = std(v_difs(:));
    % method_idx = 1;
    % fprintf('  v diff avg = %f, std = %f\n', avg_v_dif, std_v_dif);
    % for t = 2:size(velocities, 2)
    %     fprintf('  t%d: v=%f, v diff=%f\n', t, velocities(method_idx, t), abs(velocities(method_idx, t) - velocities(method_idx, t-1)));

    %     if abs(velocities(method_idx, t) - velocities(method_idx, t-1)) > (avg_v_dif + 2*std_v_dif)
    %         fprintf('  >> too large\n');
    %         velocities(method_idx, t) = velocities(method_idx, t-1);
    %     end
    % end


    %% velocity
    avg_v = cumsum(velocities, 2) ./ (repmat([1:size(velocities,2)], size(velocities,1), 1));
    traces = avg_v .* repmat(T, size(velocities,1), 1);
    traces = -traces;


    %% plot the center frequency over time
    if PLOT_FREQ_SHIFT
        for mi = 1:size(peak_freq, 1)
            fprintf('md%d: move distance = %f\n', mi, traces(mi, end));

            this_filename = ['./tmp/' filename '.freq.method' int2str(mi) '.f' num2str(f_ctr)];
            plot_center_freq(T, peak_freq(mi, :), peak_freq_power(mi, :), traces(mi, :), this_filename);
        end
    end

end


%% weighted_freq: weighted center
function [center_freq, avg_rss] = weighted_freq(freq, rss)
  avg_rss = mean(rss);
  center_freq = sum(freq .* (rss/sum(rss)));
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

                
            
%% find_spike_range: function description
function [min_ind, max_ind] = find_spike_range(ind, rss)
    spike_thresh = 3; %% ignore small spike
    min_ind = ind;
    max_ind = ind;

    next_val = rss(ind);
    for t = [ind-1:-1:1]
        val = rss(t) - spike_thresh;
        if val <= next_val
            next_val = rss(t);
        else
            min_ind = t + 1;
            break;
        end
    end

    prev_val = rss(ind);
    for t = [ind+1:1:length(rss)]
        val = rss(t) - spike_thresh;
        if val <= prev_val
            prev_val = rss(t);
        else
            max_ind = t - 1;
            break;
        end
    end

end


%% doppler_velocity: function description
% function [vr, vs] = doppler_velocity(f, f0, vr, vs)
function [vr, vs] = doppler_velocity(f0, f_offset, vr, vs)
    %% f = (c + vr) / (c + vs) * f0
    %%    c: sound speed
    %%    vr: velocity of the receiver (positive if moving toward, negative if otherwise)
    %%    vs: velocity of the sender (positive if moving away, negative if otherwise)
    c = 331 + 0.6 * 26;

    %% ------------------
    %% XXX: filtering
    % idx = find(abs(f_offset) <= 1);
    % f_offset(idx) = 0;
    %% ------------------

    f = f0 + f_offset;

    if(vr < 0)
        vr = f / f0 * (c + vs) - c;
    end

    if(vs < 0)
        vs = (c + vr) * f0 / f - c;
    end
end

    
%% plot_center_freq: function description
function plot_center_freq(time, freqs, powers, dists, filename)
    clf;
    fh = figure;
    font_size = 18;

    %% frequency
    subplot(3,1,1)
    lh1 = plot(time, freqs);
    set(lh1, 'Color', 'r');      %% color : r|g|b|c|m|y|k|w|[.49 1 .63]
    set(lh1, 'LineStyle', '-');  %% line  : -|--|:|-.
    set(lh1, 'LineWidth', 2);
    set(lh1, 'marker', '.');     %% marker: +|o|*|.|x|s|d|^|>|<|p|h
    set(lh1, 'MarkerSize', 5);
    % xlabel('Time (s)', 'FontSize', font_size);
    ylabel('Freq. (Hz)', 'FontSize', font_size);
    % ylim([f0-20 f0+20])
    
    %% Power
    subplot(3,1,2)
    bh1 = bar(time, powers);
    set(bh1, 'BarWidth', 0.6);
    set(bh1, 'EdgeColor', 'none');  %% color : r|g|b|c|m|y|k|w|[.49 1 .63]
    set(bh1, 'FaceColor', 'g');
    set(bh1, 'LineStyle', '-');  %% line  : -|--|:|-.
    set(bh1, 'LineWidth', 2);
    % xlabel('Time (s)', 'FontSize', font_size);
    ylabel('Power (dB)', 'FontSize', font_size);

    %% frequency
    subplot(3,1,3)
    lh3 = plot(time, dists);
    set(lh3, 'Color', 'r');      %% color : r|g|b|c|m|y|k|w|[.49 1 .63]
    set(lh3, 'LineStyle', '-');  %% line  : -|--|:|-.
    set(lh3, 'LineWidth', 2);
    set(lh3, 'marker', '.');     %% marker: +|o|*|.|x|s|d|^|>|<|p|h
    set(lh3, 'MarkerSize', 5);
    xlabel('Time (s)', 'FontSize', font_size);
    ylabel('Dist. (m)', 'FontSize', font_size);
    

    % set(gca, 'OuterPosition', [0 0 1 0.5]);

    print(fh, '-dpsc', [filename '.ps']);

end




%% remove_anomaly: function description
function [ts] = remove_anomaly(ts)
    avg_difs = mean2(abs(ts(:, 2:end) - ts(:, 1:end-1)));
    for ti = 2:size(ts, 2)
        idx = find(abs(ts(:, ti) - ts(:, ti-1)) > 2 * avg_difs);
        if length(idx) > 0
            ts(idx, ti) = ts(idx, ti-1);
        end
    end

end

