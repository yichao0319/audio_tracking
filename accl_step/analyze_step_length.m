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

function analyze_step_length()
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
    output_dir = './fig/';
    % filename = '0820.exp1.accl.walk.50m';
    % filename = '0820.exp2.accl.walk.55m';
    % filename = '0820.exp3.accl.walk.45m';
    % filename = '0820.exp4.accl.walk.50m';
    % filename = '0820.exp5.accl.walk.square';
    filenames = {'0820.exp1.accl.walk.50m', ...
                 '0820.exp4.accl.walk.50m', ...
                };
    avg_k = 0;
    fig_idx = 0;


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, arg = 1; end
    if nargin < 1, arg = 1; end


    %% --------------------
    %% Main starts
    %% --------------------
    
    for filei = 1:length(filenames)
        filename = char(filenames(filei));
        [dist] = get_file_info(filename);


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

        % fig_idx = fig_idx + 1;
        % fh = figure(fig_idx); clf;
        % subplot(3, 1, 1);
        % plot(data(12,:), data(13,:));
        % subplot(3, 1, 2);
        % plot(data(12,:), data(14,:));
        % subplot(3, 1, 3);
        % plot(data(12,:), data(15,:));

        % fig_idx = fig_idx + 1;
        % fh = figure(fig_idx); clf;
        % subplot(5, 1, 1);
        % plot(data(1,:), data(2,:));
        % subplot(5, 1, 2);
        % plot(data(1,:), data(3,:));
        % subplot(5, 1, 3);
        % plot(data(1,:), data(4,:));
        % subplot(5, 1, 4);
        % plot(data(1,:), data(5,:));
        % subplot(5, 1, 5);
        % plot(data(1,:), data(6,:));


        %% --------------------
        %% Remove DC (zero frequency component)
        %% --------------------
        if DEBUG2, fprintf('Remove DC\n'); end

        [B,A] = butter(10, 0.15/(Fs/2), 'high');
        for fi = [2:7, 9:11 13:15]
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
            
            [S,F,T,P] = spectrogram(data(9,:), window, noverlap, Nfft, Fs);

            fig_idx = fig_idx + 1;
            fh = figure(fig_idx); clf;
            imagesc(T, F, 10*log10(P)); % frequency-time Plot of the signal
            colorbar;
            % ylim([f_min f_max]);
            xlabel('Time (s)');
            ylabel('Power/Frequency (dB/Hz)');
            % print(fh, '-dpsc', [output_dir filename '.psd.eps']);
        end
        

        %% --------------------
        %% Find Periodic Portion and Steps
        %% --------------------
        % if DEBUG2, fprintf('Find Periodic Portion and Steps\n'); end

        ts = data(9,:);
        ts_time = data(8,:);
        ts_len = length(ts);
        [step_idx, periodic_idx, autocorr, fig_idx] = get_step_idx(ts, [], fig_idx);
        ts_duration = ts_time(periodic_idx(2)) - ts_time(periodic_idx(1));
        fprintf('  periodic portion duration = %f\n', ts_duration);
        
        
        %% plot ts and corrcoef
        fig_idx = fig_idx + 1;
        fh = figure(fig_idx); clf;
        subplot(2,1,1);
        plot(ts_time, ts);
        hold on;
        plot(ts_time(periodic_idx), ts(periodic_idx), 'ro');
        hold on;
        plot(ts_time(step_idx), ts(step_idx), 'g.');
        set(gca, 'XLim', [0 ts_time(end)]);
        subplot(2,1,2);
        plot(ts_time(1:length(autocorr)), autocorr);
        hold on;
        plot(ts_time(periodic_idx), autocorr(periodic_idx), 'ro');
        hold on;
        plot(step_idx, autocorr(step_idx), 'g.');
        set(gca, 'XLim', [0 ts_time(end)]);


        %% --------------------
        %% Step Length
        %% --------------------
        if DEBUG2, fprintf('Step Length\n'); end
            
        step_count = length(step_idx)-1;
        k = [];
        est_dist = 0;
        est_err = [];
        for si = 1:step_count
            std_idx = step_idx(si);
            end_idx = step_idx(si+1)-1;
            
            std_time = ts_time(std_idx);
            end_time = ts_time(end_idx);
            Ts = end_time - std_time;
            freq = 1 / Ts;

            step_len = dist * Ts / ts_duration;
            
            a_max = max(ts(std_idx:end_idx));
            a_min = min(ts(std_idx:end_idx));

            v = var(ts(std_idx:end_idx));

            fprintf('  step%d: Ts=%f, freq=%f, dist=%f, accl=[%f-%f], var=%f\n', si, Ts, freq, step_len, a_min, a_max, v);
            % input('')

            %% Step Length = K * (a_max - a_min )^(1/4)
            k(si) = step_len / ((a_max - a_min) ^ (1/4));

            if(avg_k ~= 0)
                est_step_len = avg_k * ((a_max - a_min) ^ (1/4));
                est_err(si) = abs(est_step_len - step_len) / step_len;
                est_dist = est_dist + est_step_len;
            end
        end
        avg_k = mean(k);

        if est_dist ~= 0
            err = abs(est_dist - dist) / dist;
            fprintf('  error = %f (per step error = %f)\n', err, mean(est_err));
        end
    end

    
end


function [dist] = get_file_info(filename)
    % filename = '0820.exp1.accl.walk.50m';
    % filename = '0820.exp2.accl.walk.55m';
    % filename = '0820.exp3.accl.walk.45m';
    % filename = '0820.exp4.accl.walk.50m';
    % filename = '0820.exp5.accl.walk.square';
    dist = 0;
    if strcmp(filename, '0820.exp1.accl.walk.50m')
        dist = 50;
    elseif strcmp(filename, '0820.exp2.accl.walk.55m')
        dist = 55;
    elseif strcmp(filename, '0820.exp3.accl.walk.45m')
        dist = 45;
    elseif strcmp(filename, '0820.exp4.accl.walk.50m')
        dist = 50;
    elseif strcmp(filename, '0820.exp5.accl.walk.square')
        dist = 150; %% not sure..
    end
end




