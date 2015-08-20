%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%% MAIN
%% - Input:
%%
%%
%% - Output:
%%
%%
%% example:
%%     analyze_camera_sound('07.04.8000.40-10');
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function analyze_camera_sound(filename, std_x, std_y)
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
    ctr_freq = get_cheat_freq(filename);
    font_size = 16;


    %% --------------------
    %% Variable
    %% --------------------
    input_dir_cam = '../../data/camera_tracking/traces/';
    input_dir_sound  = './rx_sound/';
    output_dir = './fig/';


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 2, std_x = 0; end
    if nargin < 3, std_y = 0; end
    std_pos = [std_x, std_y];


    %% --------------------
    %% Main starts
    %% --------------------

    %% ===========================
    %% parse filename to get basic infomation
    %% ===========================
    if DEBUG2, fprintf('parse filename to get basic infomation\n'); end

    [freq, spk] = parse_info(filename);

    fprintf('  freq = %d\n', freq);
    fprintf('  spk position = (%d,%d)\n', spk(1), spk(2));

    
    %% ===========================
    %% Read Camera Trace
    %% ===========================
    if DEBUG2, fprintf('Read Camera Trace\n'); end

    [cam_time, cam_trace_world, cam_trace_img] = read_camera_trace([input_dir_cam filename '.txt'], std_pos);


    %% ===========================
    %% Get distnace and speed for simulation
    %% ===========================
    if DEBUG2, fprintf('Get distnace and speed for simulation\n'); end

    [cam_dist_to_spk, cam_shift, cam_speed, cam_direction] = retrieve_dist_speed(cam_time, cam_trace_world, spk);

    
    % fh = figure(1); clf;
    % plot(cam_trace_world(:,1), cam_trace_world(:,2));
    % hold on;
    % lspk = plot(spk(1), spk(2), 'r*');
    % set(lspk, 'MarkerSize', 10);

    % set(gca, 'XLim', [-0.05 0.4]);
    % set(gca, 'YLim', [-0.05 0.25]);
    % ylabel('meters', 'FontSize', font_size);
    % xlabel('meters', 'FontSize', font_size);
    % set(gca, 'FontSize', font_size);

    
    %% ===========================
    %% get position shift calculated by Doppler Effect
    %% ===========================
    if DEBUG2, fprintf('get position shift calculated by Doppler Effect\n'); end

    [sound_time, sound_v, sound_shift] = get_doppler_shift(input_dir_sound, filename, ctr_freq);
    sound_time = sound_time';
    sound_v = sound_v';
    sound_shift = sound_shift';


    %% ===========================
    %% Sync camera and sound data
    %% ===========================
    if DEBUG2, fprintf('Sync camera and sound data\n'); end

    [cam_std_time, sound_std_time] = sync_camera_sound(cam_time, cam_shift, sound_time, sound_shift);
    
    cam_time_sync = cam_time - cam_std_time;
    idx = find(cam_time_sync >= 0);
    cam_time_sync = cam_time_sync(idx);
    cam_shift_sync = cam_shift(idx);
    cam_dist_to_spk_sync = cam_dist_to_spk(idx);
    cam_speed_sync = cam_speed(idx);
    cam_trace_world_sync = cam_trace_world(idx, :);

    sound_time_sync = sound_time - sound_std_time;
    idx = find(sound_time_sync >= 0);
    sound_time_sync = sound_time_sync(idx);
    sound_shift_sync = sound_shift(idx);
    sound_v_sync = sound_v(idx);


    fh = figure(2); clf;
    subplot(2,1,1);
    plot(cam_time-cam_time(1), cam_dist_to_spk);

    subplot(2,1,2);
    plot(cam_time_sync, cam_shift_sync, 'b');
    hold on;
    plot(sound_time_sync, sound_shift_sync, 'r');


    %% ===========================
    %% Infer 2D location
    %% ===========================
    if DEBUG2, fprintf('Infer 2D location\n'); end
    
    time = sound_time_sync;
    speed = sound_v_sync;
    real_trace = ...
        [interp1(cam_time_sync, cam_trace_world_sync(:,1), time, 'linear', 'extrap'), ...
         interp1(cam_time_sync, cam_trace_world_sync(:,2), time, 'linear', 'extrap')];
    % dist = interp1(cam_time_sync, cam_dist_to_spk_sync, time, 'linear', 'extrap');
    [cam_dist_to_spk2, cam_shift2, cam_speed2, cam_direction2] = retrieve_dist_speed(time, real_trace, spk);
    dist = cam_dist_to_spk2;
    direction = cam_direction2;


    % [est_trace] = cal_trace_dist_speed(spk, time, dist, speed, direction, real_trace);
    % [est_traces] = cal_trace_dist_speed_particle(spk, time, dist, speed, direction, real_trace);
    [est_traces] = cal_trace_dist_speed_opt(spk, time, dist, speed, direction, real_trace);

    for pidx = 1:size(est_traces,2)/2
        fprintf('  fig: %d/%d\n', pidx, size(est_traces,2)/2);
        % size(est_traces(:,(pidx-1)*2+1:pidx*2))
        % size(real_trace)
        avg_err = mean(cal_dist(est_traces(:,(pidx-1)*2+1:pidx*2), real_trace));
        fprintf('  avg error = %f\n', avg_err);

        fh = figure(1); clf;
        plot(real_trace(:,1), real_trace(:,2), '-bo');
        hold on;
        lspk = plot(spk(1), spk(2), 'r*');
        set(lspk, 'MarkerSize', 10);
        hold on;
        plot(est_traces(:,(pidx-1)*2+1), est_traces(:,pidx*2), '-r.');

        % set(gca, 'XLim', [-0.05 0.4]);
        % set(gca, 'YLim', [-0.05 0.25]);
        ylabel('meters', 'FontSize', font_size);
        xlabel('meters', 'FontSize', font_size);
        set(gca, 'FontSize', font_size);
        
        % print(fh, '-dpsc', [output_dir filename '.trace.eps']);

        waitforbuttonpress
    end



    %% ----------------------
    %% DEBUG
    % cam_speed_sync2 = interp1(cam_time_sync, cam_speed_sync, time);
    % [cam_dist_to_spk2, cam_shift2, cam_speed2, cam_direction2] = retrieve_dist_speed(time, real_trace, spk);

    % [est_trace] = cal_trace_dist_speed(spk, time, cam_dist_to_spk2, cam_speed2, cam_direction2, real_trace);
    % avg_err = mean(sqrt(sum((est_trace - real_trace) .^ 2, 2)));
    % fprintf('  avg error = %f\n', avg_err);


    % fh = figure(3); clf;
    % plot(real_trace(:,1), real_trace(:,2), '-bo');
    % hold on;
    % lspk = plot(spk(1), spk(2), 'r*');
    % set(lspk, 'MarkerSize', 10);
    % hold on;
    % plot(est_trace(:,1), est_trace(:,2), '-r.');

    % set(gca, 'XLim', [-0.05 0.4]);
    % set(gca, 'YLim', [-0.05 0.25]);
    % ylabel('meters', 'FontSize', font_size);
    % xlabel('meters', 'FontSize', font_size);
    % set(gca, 'FontSize', font_size);
    
    % print(fh, '-dpsc', [output_dir filename '.trace.eps']);

    
end



function [cam_time, cam_trace_world, cam_trace_img] = ...
    read_camera_trace(filename, std_pos)

    thresh = 0.01;

    data = load(filename);

    ok_idx = find(~isnan(data(:,5)));
    data = data(ok_idx(1):ok_idx(end), :);

    dist = sqrt((data(:,4)-std_pos(1)).^2 + (data(:,5)-std_pos(2)).^2);
    std_idx = find(dist < thresh);
    data = data(std_idx(1):end, :);


    nan_idx = find(~isnan(data(:,5)));
    for di = nan_idx
        if isnan(data(di,5))
            cur_time = data(di,1);

            prev_time = data(di-1,1);
            prev_loc = data(di-1, 2:5);

            next_ok_idx = find(isnan(data(di+1:end,5)));
            next_ok_idx = next_ok_idx(1);

            next_time = data(next_ok_idx,1);
            next_loc = data(next_ok_idx, 2:5);

            data(di, 2:5) = (cur_time - prev_time) * (next_loc - prev_loc) / (next_time - prev_time) + prev_loc;
        end
    end
    
    cam_time = data(:, 1);
    cam_trace_img = data(:, 2:3);
    cam_trace_world = data(:, 4:5);
end



function [freq, spk] = parse_info(filename)
    expression = ['\d+\d+.(?<freq>\d+).(?<cam_x>\d+)-(?<cam_y>\d+)'];
    tokenNames = regexp(filename, expression, 'names');
    if(length(tokenNames) > 0)
        % freq = str2num(tokenNames(1).freq) * 1000;
        freq = str2num(tokenNames(1).freq);
        spk = [str2num(tokenNames(1).cam_x), str2num(tokenNames(1).cam_y)] / 100;
    else
        error(['wong format: ' filename]);
    end
end





function ctr_freq = get_cheat_freq(filename)
    if strcmp(filename, '07.04.8000.40-10')
        ctr_freq = 7999.4;
    else
        error('wrong filename');
    end
end


%% dist: distance to the speaker
%% shift: relative shift to the speaker
%% speed: relative speed to the speaker
%% direction: moving direction
function [dist, shift, speed, direction] = retrieve_dist_speed(cam_time, cam_trace_world, spk_pos)
    dist = cal_dist(cam_trace_world, spk_pos);
    shift = dist - dist(1);
    
    velocity = (cam_trace_world(2:end, :) - cam_trace_world(1:end-1, :)) ./ [(cam_time(2:end) - cam_time(1:end-1)), (cam_time(2:end) - cam_time(1:end-1))];
    vec2spk = repmat(spk_pos, size(cam_trace_world,1)-1, 1) - cam_trace_world(1:end-1, :);
    speed = [(velocity(:,1) .* vec2spk(:,1) + velocity(:,2) .* vec2spk(:,2)) ./ (sqrt(vec2spk(:,1).^2 + vec2spk(:,2).^2)); 0];

    % direction = [cam_trace_world(2:end,:) - cam_trace_world(1:end-1,:); 0, 0];
    direction = velocity;
    % norm_dir = sqrt(direction(:,1) .^ 2 + direction(:,2) .^ 2);
    % direction = direction ./ repmat(norm_dir, 1, 2);
end


function [cam_std_time, sound_std_time] = sync_camera_sound(cam_time, cam_shift, sound_time, sound_shift)

    tic;

    %% undersampling
    itv = 0.01;
    cam_time_under = [cam_time(1):itv:cam_time(end)];
    cam_shift_under = interp1(cam_time, cam_shift, cam_time_under);

    sound_time_under = [sound_time(1):itv:sound_time(end)];
    sound_shift_under = interp1(sound_time, sound_shift, sound_time_under);


    %% sync by finding the shift with max corrcoef
    cam_len = length(cam_time_under);
    sound_len = length(sound_time_under);

    max_r = -1;
    max_t = 1;
    for ti = 1:cam_len - sound_len + 1
        r = corrcoef(cam_shift_under(ti:ti+sound_len-1), sound_shift_under);
        if r(1,2) > max_r
            max_r = r(1,2);
            max_t = ti;
        end
    end


    %% sync
    % cam_time_sync = cam_time - cam_time_under(max_t);
    % idx = find(cam_time_sync >= 0);
    % cam_time_sync = cam_time_sync(idx);
    % cam_shift_sync = cam_shift(idx);
    
    % sound_time_sync = sound_time - sound_time_under(1);
    % sound_shift_sync = sound_shift;

    cam_std_time = cam_time_under(max_t);
    sound_std_time = sound_time_under(1);


    elapse = toc;
    fprintf('  running time: %f\n', elapse);
    fprintf('  max corrcoef = %f\n', max_r);

end

