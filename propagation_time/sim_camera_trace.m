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
%%     sim_camera_trace('07.04.8000.40-10');
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sim_camera_trace(filename, std_x, std_y)
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
    
    fh = figure(2); clf;
    subplot(3,1,1);
    plot(cam_time, cam_dist_to_spk);
    ylabel('dist2spk', 'FontSize', font_size);
    set(gca, 'FontSize', font_size);
    
    subplot(3,1,2);
    plot(cam_time, cam_shift);
    ylabel('movement', 'FontSize', font_size);
    set(gca, 'FontSize', font_size);

    subplot(3,1,3);
    plot(cam_time, cam_speed);
    ylabel('speed', 'FontSize', font_size);
    set(gca, 'FontSize', font_size);


    %% ===========================
    %% Calculate trace
    %% ===========================
    if DEBUG2, fprintf('Calculate trace\n'); end

    % [est_traces] = cal_trace_dist_speed(spk, cam_time, cam_dist_to_spk, cam_speed, cam_direction, cam_trace_world);
    % [est_traces] = cal_trace_dist_speed_particle_prob(spk, cam_time, cam_dist_to_spk, cam_speed, cam_direction, cam_trace_world);
    % [est_traces] = cal_trace_dist_speed_particle(spk, cam_time, cam_dist_to_spk, cam_speed, cam_direction, cam_trace_world);
    % [est_traces] = cal_trace_dist_speed_opt(spk, cam_time, cam_dist_to_spk, cam_speed, cam_direction, cam_trace_world);
    [est_traces] = cal_trace_dist_speed_opt_multi_snapshot(spk, cam_time, cam_dist_to_spk, cam_speed, cam_direction, cam_trace_world);
    
    fh = figure(1); 

    for pidx = 1:size(est_traces,2)/2
        fprintf('  fig: %d/%d\n', pidx, size(est_traces,2)/2);
        avg_err = mean(cal_dist(est_traces(:,(pidx-1)*2+1:pidx*2), cam_trace_world));
        med_err = median(cal_dist(est_traces(:,(pidx-1)*2+1:pidx*2), cam_trace_world));
        max_err = max(cal_dist(est_traces(:,(pidx-1)*2+1:pidx*2), cam_trace_world));
        min_err = min(cal_dist(est_traces(:,(pidx-1)*2+1:pidx*2), cam_trace_world));
        fprintf('  Error: avg=%f, median=%f, max=%f, min=%f\n', avg_err, med_err, max_err, min_err);

        avg_err = mean(cal_dist(est_traces(floor(end/2):end,(pidx-1)*2+1:pidx*2), cam_trace_world(floor(end/2):end,:)));
        med_err = median(cal_dist(est_traces(floor(end/2):end,(pidx-1)*2+1:pidx*2), cam_trace_world(floor(end/2):end,:)));
        max_err = max(cal_dist(est_traces(floor(end/2):end,(pidx-1)*2+1:pidx*2), cam_trace_world(floor(end/2):end,:)));
        min_err = min(cal_dist(est_traces(floor(end/2):end,(pidx-1)*2+1:pidx*2), cam_trace_world(floor(end/2):end,:)));
        fprintf('  Error: avg=%f, median=%f, max=%f, min=%f\n', avg_err, med_err, max_err, min_err);

        clf;
        plot(cam_trace_world(:,1), cam_trace_world(:,2), '-bo');
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

        % waitforbuttonpress
    end
    
    % print(fh, '-dpsc', [output_dir filename '.trace.eps']);
    print(fh, '-dpng', [output_dir filename '.trace.png']);

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
    expression = ['\d+\d+.(?<freq>\d+).(?<cam_x>\d+)-(?<cam_y>\d+)\.*.*'];
    tokenNames = regexp(filename, expression, 'names');
    if(length(tokenNames) > 0)
        % freq = str2num(tokenNames(1).freq) * 1000;
        freq = str2num(tokenNames(1).freq);
        spk = [str2num(tokenNames(1).cam_x), str2num(tokenNames(1).cam_y)] / 100;
    else
        error(['wong format: ' filename]);
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

    direction = [cam_trace_world(2:end,:) - cam_trace_world(1:end-1,:); 0, 0];
    % direction = velocity / norm(velocity);
    norm_dir = sqrt(direction(:,1) .^ 2 + direction(:,2) .^ 2);
    direction = direction ./ repmat(norm_dir, 1, 2);
end


