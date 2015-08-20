%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%% MAIN
%% Given
%%   a) distance to the speaker over time
%%   b) acceleration
%%
%% Input:
%%   dist_err: meters
%%   v_err:    m/s
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sim_shake_positioning(dist_err, v_err)
    addpath('../propagation_time');

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
    spk_pos = [5, 5];
    dist_itv = 0.1;
    dist_f = 1 / dist_itv;

    
    %% --------------------
    %% Variable
    %% --------------------
    % input_dir  = '';
    output_dir = './fig/';


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, dist_err = 0; end
    if nargin < 2, v_err = 0; end


    %% --------------------
    %% Main starts
    %% --------------------
    
    %% --------------------
    %% real trace
    %% --------------------
    if DEBUG2, fprintf('Generate real trace\n'); end

    [time, real_trace] = gen_real_trace();

    move_dist = cal_dist(real_trace(end,:), real_trace(1,:));
    fprintf('  moving distance = %fcm\n', move_dist*100);

    
    %% --------------------
    %% infer distance to speaker, acceleration, etc info
    %% --------------------
    if DEBUG2, fprintf('Infer distance to speaker, acceleration, etc info\n'); end

    [real_dist2spk, real_velocity] = retrieve_trace_info(time, real_trace, spk_pos);


    %% --------------------
    %% Inject error
    %% --------------------
    if DEBUG2, fprintf('Inject error\n'); end

    dist2spk = real_dist2spk + dist_err * (rand(length(real_dist2spk),1) * 2 - 1);
    velocity = real_velocity + v_err * (rand(length(real_dist2spk),2) * 2 - 1);


    %% --------------------
    %% Synchronization samples from different sensors
    %% --------------------
    % if DEBUG2, fprintf('Synchronization samples from different sensors\n'); end


    %% --------------------
    %% Calculate traces
    %% --------------------
    if DEBUG2, fprintf('Calculate traces\n'); end
    
    est_traces = cal_trace_dist_acc(time, dist2spk, velocity, real_trace, spk_pos);


    fh = figure(1); 
    for pidx = 1:size(est_traces,2)/2
        this_est_trace = est_traces(:,(pidx-1)*2+1:pidx*2);
        fprintf('  fig %d/%d\n', pidx, size(est_traces,2)/2);
        
        err = mean(cal_dist(real_trace, this_est_trace));
        fprintf('    error = %f\n', err);

        avg_err = cal_dist(real_trace(1,:), mean(this_est_trace,1));
        fprintf('    avg error = %f\n', avg_err);

        clf;
        plot(real_trace(:,1), real_trace(:,2), '-bo');
        hold on;
        plot(this_est_trace(:,1), this_est_trace(:,2), '-r.');
        hold on;
        lh = plot(mean(this_est_trace(:,1)), mean(this_est_trace(:,2)), 'y.');
        set(lh, 'MarkerSize', 20);
        hold on;
        lh = plot(spk_pos(1), spk_pos(2), 'g*');
        set(lh, 'MarkerSize', 10);


        waitforbuttonpress
    end


end


function [time, real_trace] = gen_real_trace(spk_pos)
    f = 1/0.01;
    len = 0.4;  %% seconds
    time = [0:1/f:len]';
    ns = length(time);
    
    real_trace = [0, 0];
    speed = 0.2;  %% 20cm/s
    velocity = [3, 4];
    velocity = speed * velocity / norm(velocity);
    dir_update = [1, 1];

    for ti = 2:ns
        new_pos = real_trace(end, :) + velocity * (time(ti)-time(ti-1));
        real_trace = [real_trace; new_pos];

        velocity = velocity .* dir_update;
    end

end


function [dist2spk, velocity] = retrieve_trace_info(time, real_trace, spk_pos)
    dist2spk = cal_dist(real_trace, spk_pos);
    direction = [real_trace(2:end, :)-real_trace(1:end-1, :); 0, 0];
    velocity = [direction(1:end-1, :) ./ repmat(time(2:end)-time(1:end-1), 1, 2); 0, 0];
end


function [est_traces] = cal_trace_dist_acc(time, dist2spk, velocity, real_trace, spk_pos)
    epsilon = 10^-5; 
    est_traces = real_trace(1, :);
    first = 1;

    for ti = 2:length(time)
        a = spk_pos(1);
        b = spk_pos(2);
        r1 = dist2spk(ti-1);
        r2 = dist2spk(ti);
        d = norm(velocity(ti-1,:) * (time(ti)-time(ti-1)));
        dir_x = velocity(ti-1,1);
        dir_y = velocity(ti-1,2);

        if(abs(r1-r2) < epsilon)
            % fprintf('same location\n');
            x1 = est_traces(end, 1);
            y1 = est_traces(end, 2);
            x2 = x1;
            y2 = y1;
        else
            % fprintf('new location: %.10f\n', r1-r2);
            [x1, y1, x2, y2] = cal_intersections(a, b, r1, r2, d, dir_x, dir_y);
        end

        if first & isreal(x1(1)) & isreal(y1(1)) & isreal(x2(1)) & isreal(y2(1)) & ...
           isreal(x1(2)) & isreal(y1(2)) & isreal(x2(2)) & isreal(y2(2))
            % first = 0;
            fh = figure(5); clf;
            plot(real_trace(:,1), real_trace(:,2), '-bo');
            hold on;
            lh = plot(spk_pos(1), spk_pos(2), 'g*');
            set(lh, 'MarkerSize', 10);
            hold on;
            fprintf('  # points: %d\n', length(x1));
            for li = 1:length(x1)
                plot([x1(li) x2(li)], [y1(li) y2(li)], '-r.');
                hold on;
            end
        end


        li = 1;
        min_idx = li;
        [x1(li), y1(li)];
        min_dist = cal_dist(est_traces(ti-1,:), [x1(li), y1(li)]);
        for li = 2:length(x1)
            if cal_dist(est_traces(ti-1,:), [x1(li), y1(li)]) < min_dist
                min_idx = li;
                min_dist = cal_dist(est_traces(ti-1,:), [x1(li), y1(li)]);
            end
        end

        est_traces(ti, :) = [x2(min_idx), y2(min_idx)];
        % return;
    end

    est_traces = real(est_traces);
end


