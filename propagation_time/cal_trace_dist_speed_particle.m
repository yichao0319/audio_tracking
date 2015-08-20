%% dist: sound
%% speed: doppler
%% direction: accelerometer
function [est_trace, est_traces] = cal_trace_dist_speed_particle(spk_pos, time, dist, speed, direction, real_trace)

    REAL_TIME_FIGURE = 1;

    vel_err = 0.00;
    dist_err = 0.00;
    dist_err_thresh = 0.1; %10;
    merge_thresh = 0.00;
    wait_time = 0.05;

    %% GIVEN: init_pos, dist, speed

    %% get relative velocity
    vec2spk = spk_pos - real_trace(1,:);
    vec2spk = vec2spk / norm(vec2spk);
    rel_velocity = speed(1) * vec2spk;
    
    est_trace = [real_trace(1,:)];
    particles = est_trace;
    p_speed = [0];
    p_dist_err = [0];

    %% -------------
    %% real time updated figure
    if REAL_TIME_FIGURE
        callstr = ['set(gcbf,''Userdata'',double(get(gcbf,''Currentcharacter''))); uiresume'];
        run_fig = 1;
        fh = figure(5);
        clf(fh);
        set(fh, 'keypressfcn', callstr);
        set(fh, 'windowstyle', 'modal');
        set(fh, 'userdata', -1);
        tt = timer;
        tt.timerfcn = 'uiresume';
        tt.startdelay = wait_time;
    end
    %% -------------

    for ti = 2:length(time)
        %% known: real_trace(:, ti-1), speed(:, ti-1), dist(ti)
        new_particles = [];
        new_p_speed = [];
        new_p_dist_err = [];

        for pidx = 1:size(particles,2)/2
            % prev_loc = real_trace(ti-1,:);  %% XXX: for now...
            % prev_loc = est_trace(ti-1,:);
            prev_loc = particles(end, (pidx-1)*2+1:pidx*2);


            %% get relative velocity
            if ti > 2
                vec2spk(ti-1,:) = spk_pos - prev_loc;
                vec2spk(ti-1,:) = vec2spk(ti-1,:) / norm(vec2spk(ti-1,:));
                rel_velocity(ti-1,:) = speed(ti-1) * vec2spk(ti-1,:);
            end

            %% get new dist to speaker
            dist2spk = dist(ti);

            %% get moving direction
            this_direction = direction(ti-1, :);

            
            %% ---------------------
            %% inject error
            % rel_velocity(ti-1,:) = rel_velocity(ti-1,:) + vel_err * (rand(1,2)*2-1);
            % dist2spk = dist2spk + dist_err * (rand(1)*2-1);
            % rel_velocity(ti-1,:) = rel_velocity(ti-1,:) .* (1 + vel_err * (rand(1,2)*2-1));
            % dist2spk = dist2spk * (1 + dist_err * (rand(1)*2-1));
            rel_velocity(ti-1,:) = rel_velocity(ti-1,:) .* (1 + vel_err);
            dist2spk = dist2spk * (1 + dist_err);
            %% ---------------------


            %% algorithm :)
            proj_move_dist = rel_velocity(ti-1,:) * (time(ti) - time(ti-1));
            proj_pt = prev_loc + proj_move_dist;
            possible_line_vec = [-proj_move_dist(2), proj_move_dist(1)];
            slope = possible_line_vec(2) / possible_line_vec(1);

            %% proj_pt: (x0, y0)
            %% slope: m
            %% line: (x, m(x-x0)+y0)
            %% (x - spk_x)^2 + (m(x-x0)+y0 - spk_y)^2 = dist^2
            %% let w = -mx0 + y0 - spk_y
            %% (x - spk_x)^2 + (mx + w)^2 = dist^2
            %% (1+m^2)x^2 + 2(mw - spk_x)x + (spk_x^2 + w^2 - dist^2) = 0
            %% x = -b +- sqrt(b^2 - 4ac) / 2a
            a = (1 + slope^2);
            w = - slope * proj_pt(1) + proj_pt(2) - spk_pos(2);
            b = 2 * (slope * w - spk_pos(1));
            c = (spk_pos(1)^2 + w^2 - dist2spk^2);
            
            x1 = (- b + sqrt(b^2 - 4*a*c)) / 2 / a;
            y1 = slope * (x1 - proj_pt(1)) + proj_pt(2);
            x2 = (- b - sqrt(b^2 - 4*a*c)) / 2 / a;
            y2 = slope * (x2 - proj_pt(1)) + proj_pt(2);
            
            % break;

            %% ====================================
            %% assume we know the best solution
            % if norm(real_trace(ti,:) - [x1,y1]) < norm(real_trace(ti,:) - [x2,y2])
            %     est_trace(ti,:) = [x1,y1]; 
            % else
            %     est_trace(ti,:) = [x2,y2]; 
            % end
            %% ------------------------------------
            %% assume we have direction from accelerometer
            % dir1 = [x1,y1] - prev_loc;
            % dir2 = [x2,y2] - prev_loc;
            % if sum(dir1 .* this_direction) > sum(dir2 .* this_direction)
            %     est_trace(ti,:) = [x1,y1];
            % else
            %     est_trace(ti,:) = [x2,y2]; 
            % end
            %% ------------------------------------
            %% generate new particles
            dir1 = [x1,y1] - prev_loc;
            dir2 = [x2,y2] - prev_loc;
            speed1 = norm(dir1) / (time(ti)-time(ti-1));
            speed2 = norm(dir2) / (time(ti)-time(ti-1));
            dist_err1 = abs(norm([x1,y1] - spk_pos) - dist2spk);
            dist_err2 = abs(norm([x2,y2] - spk_pos) - dist2spk);

            %%    new particles
            this_new_particle1   = [particles(:, (pidx-1)*2+1:pidx*2); x1, y1];
            this_new_p_speed1    = [p_speed(:,pidx); speed1];
            this_new_p_dist_err1 = [p_dist_err(:,pidx); dist_err1];
            if dist_err1 < dist_err_thresh
                %% merge particles
                merged = 0;
                % merged = if_merged(new_particles, this_new_particle1, merge_thresh);

                if ~merged
                    new_particles  = [new_particles, this_new_particle1];
                    new_p_speed    = [new_p_speed, this_new_p_speed1];
                    new_p_dist_err = [new_p_dist_err, this_new_p_dist_err1]; 
                end
            end

            this_new_particle2   = [particles(:, (pidx-1)*2+1:pidx*2); x2, y2];
            this_new_p_speed2    = [p_speed(:,pidx); speed2];
            this_new_p_dist_err2 = [p_dist_err(:,pidx); dist_err2];
            if dist_err2 < dist_err_thresh
                %% merge particles
                merged = 0;
                % merged = if_merged(new_particles, this_new_particle2, merge_thresh);

                if ~merged
                    new_particles  = [new_particles, this_new_particle2];
                    new_p_speed    = [new_p_speed, this_new_p_speed2];
                    new_p_dist_err = [new_p_dist_err, this_new_p_dist_err2]; 
                end
            end

            % fprintf('  [%.0f,%.0f] (d=%.2f): sol1=[%.0f,%.0f]=%.2f (d=%.2f), sol2=[%.0f,%.0f]=%.2f (d=%.2f)\n', real_trace(ti,:) * 100, dist2spk * 100, x1 * 100, y1 * 100, norm(real_trace(ti,:) - [x1,y1]) * 100, norm([x1,y1]-spk_pos) * 100, x2 * 100, y2 * 100, norm(real_trace(ti,:) - [x2,y2]) * 100, norm([x2,y2]-spk_pos) * 100);
            
        end

        %% update particles and probability
        particles  = new_particles;
        p_speed    = new_p_speed;
        p_dist_err = new_p_dist_err;

        
        fprintf('  %d/%d: #particles=%d\n', ti, length(time), size(particles,2)/2);

        %% remove some particles
        n_p_thresh = 1000;
        n_p = floor(n_p_thresh / 10);
        if size(particles,2)/2 > n_p_thresh;
            fprintf('--- remove some particles ------------\n');
            [particles, p_dist_err, p_speed] = reduce_particles(particles, p_dist_err, p_speed, real_trace, n_p);
            fprintf('  %d/%d: #particles=%d\n', ti, length(time), size(particles,2)/2);
            fprintf('--- done ------------\n');
        end


        %% -------------
        %% real time updated figure
        if REAL_TIME_FIGURE
            %% choose the trajectory with highest probability
            [v, idx] = min(p_dist_err(end, :));
            est_trace = [particles(:, (idx-1)*2+1:idx*2)];

            clf(fh);
            plot(spk_pos(1), spk_pos(2), 'bo'); 
            hold on;
            plot(real_trace(1:ti,1), real_trace(1:ti,2), '-bo');
            hold on
            plot(real_trace(ti,1), real_trace(ti,2), 'b*')
            hold on;
            plot(est_trace(1:ti,1), est_trace(1:ti,2), '-g.');
            hold on
            plot(est_trace(ti,1), est_trace(ti,2), 'g*')
            % hold on;
            % set(gca, 'XLim', [-0.1 0.5]);
            % set(gca, 'YLim', [-0.1 0.5]);

            % waitforbuttonpress
            if run_fig
                try
                    start(tt);
                    uiwait;
                    ch = [];
                    ch = get(fh,'Userdata');
                catch
                    uiwait;
                    ch = [];
                end
            else
                uiwait;
                ch = get(fh, 'Userdata');
            end
                
            if length(ch) > 0 & ch(1) == 97, run_fig = 0;
            else, run_fig = 1;
            end
        end
        %% -------------

    end

    %% -------------
    %% real time updated figure
    if REAL_TIME_FIGURE
        stop(tt);
        delete(tt);
        close(fh);
    end
    %% -------------

    %% choose the trajectory with highest probability
    [v, idx] = min(p_dist_err(end, :));
    est_trace = particles(:, (idx-1)*2+1:idx*2);
    
    all_idx = find(p_dist_err(end, :) == v);
    est_traces = particles(:, sort([(all_idx-1)*2+1, all_idx*2]));

end


% function [merged] = if_merged(particles, new_pos, merge_thresh)
%     %% merge particles
%     merged = 0;
%     for pidx = 1:size(particles,2)/2
%         if norm(new_pos - particles(end,(pidx-1)*2+1:pidx*2)) <= merge_thresh
%             merged = 1;
%             break;
%         end
%     end
% end
function [merged] = if_merged(particles, new_particle, merge_thresh)
    %% merge particles
    merged = 0;
    for pidx = 1:size(particles,2)/2
        if norm(new_particle - particles(:,(pidx-1)*2+1:pidx*2)) <= merge_thresh
            merged = 1;
            break;
        end
    end
end


function [particles, p_dist_err, p_speed] = reduce_particles(particles, p_dist_err, p_speed, real_trace, num)
    % method = 'min_end_error';
    % method = 'min_avg_error';
    method = 'gt';
    % method = 'gt2';
    %% --------------------------------
    %% select particles w/ min final error
    if strcmp(method, 'min_end_error')
        [v, idx] = sort(p_dist_err(end, :), 'ascend');
        idx = idx(1:num);
        particles  = particles(:, sort([(idx-1)*2+1, idx*2]));
        p_speed    = p_speed(:, sort(idx));
        p_dist_err = p_dist_err(:, sort(idx));
    %% --------------------------------
    %% select particles w/ min avg error
    elseif strcmp(method, 'min_avg_error')
        [v, idx] = sort(sum(p_dist_err(:,:),1), 'ascend');
        idx = idx(1:num);
        particles  = particles(:, sort([(idx-1)*2+1, idx*2]));
        p_speed    = p_speed(:, sort(idx));
        p_dist_err = p_dist_err(:, sort(idx));
    %% --------------------------------
    %% select best particles by ground truth
    elseif strcmp(method, 'gt')
        len = size(particles,1);

        for pidx = 1:size(particles,2)/2
            % err(pidx) = norm(particles(:,(pidx-1)*2+1:pidx*2) - real_trace(1:len,:));
            err(pidx) = mean(cal_dist(particles(:,(pidx-1)*2+1:pidx*2), real_trace(1:len,:)));
        end

        [v, idx] = sort(err, 'ascend');
        idx = idx(1:num);
        particles  = particles(:, sort([(idx-1)*2+1, idx*2]));
        p_speed    = p_speed(:, sort(idx));
        p_dist_err = p_dist_err(:, sort(idx));
    %% --------------------------------
    %% select best particles by min end error and ground truth
    elseif strcmp(method, 'gt2')
        [v, idx] = min(sum(p_dist_err,1));
        idx = find(sum(p_dist_err,1) == v);
        particles  = particles(:, sort([(idx-1)*2+1, idx*2]));
        p_speed    = p_speed(:, sort(idx));
        p_dist_err = p_dist_err(:, sort(idx));

        if size(particles,2)/2 > num
            len = size(particles,1);

            for pidx = 1:size(particles,2)/2
                err(pidx) = norm(particles(:,(pidx-1)*2+1:pidx*2) - real_trace(1:len,:));
            end

            [v, idx] = sort(err, 'ascend');
            idx = idx(1:num);
            particles  = particles(:, sort([(idx-1)*2+1, idx*2]));
            p_speed    = p_speed(:, sort(idx));
            p_dist_err = p_dist_err(:, sort(idx));
        end
    end
end


