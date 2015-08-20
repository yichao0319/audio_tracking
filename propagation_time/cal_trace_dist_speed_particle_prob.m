%% dist: sound
%% speed: doppler
%% direction: accelerometer
function [est_trace] = cal_trace_dist_speed_particle_prob(spk_pos, time, dist, speed, direction, real_trace)
    vel_err = 0.00;
    dist_err = 0.00;

    %% GIVEN: init_pos, dist, speed

    %% get relative velocity
    vec2spk = spk_pos - real_trace(1,:);
    vec2spk = vec2spk / norm(vec2spk);
    rel_velocity = speed(1) * vec2spk;
    
    est_trace = [real_trace(1,:)];
    particles = est_trace;
    p_prob = [1];

    % fh = figure(5); clf(fh);

    for ti = 2:length(time)
        %% known: real_trace(:, ti-1), speed(:, ti-1), dist(ti)
        new_particles = [];
        new_p_prob = [];

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
            dir1 = [x1,y1] - prev_loc; norm_dir1 = dir1 / norm(dir1);
            dir2 = [x2,y2] - prev_loc; norm_dir2 = dir2 / norm(dir2);
            likelyhood_dir1 = sum(norm_dir1 .* this_direction);
            likelyhood_dir2 = sum(norm_dir2 .* this_direction);
            prob_dir1 = likelyhood_dir1 / (likelyhood_dir1+likelyhood_dir2);
            prob_dir2 = likelyhood_dir2 / (likelyhood_dir1+likelyhood_dir2);
            if prob_dir1 < 0 & prob_dir1 < prob_dir2
                prob_dir1 = 0.001;
                prob_dir2 = 0.999;
            elseif prob_dir2 < 0 & prob_dir2 < prob_dir1
                prob_dir1 = 0.999;
                prob_dir2 = 0.001;
            end

            speed1 = norm(dir1) / (time(ti)-time(ti-1));
            speed2 = norm(dir2) / (time(ti)-time(ti-1));
            likelyhood_speed1 = abs(speed1 - norm(this_direction));
            likelyhood_speed2 = abs(speed2 - norm(this_direction));
            prob_speed1 = likelyhood_speed1 / (likelyhood_speed1+likelyhood_speed2);
            prob_speed2 = likelyhood_speed2 / (likelyhood_speed1+likelyhood_speed2);


            this_prob1 = prob_dir1 / (prob_dir1+prob_dir2);
            this_prob2 = prob_dir2 / (prob_dir1+prob_dir2);
            % this_prob1 = prob_speed1 / (prob_speed1+prob_speed2);
            % this_prob2 = prob_speed2 / (prob_speed1+prob_speed2);

            %%    new particles
            new_p1 = [particles(:, (pidx-1)*2+1:pidx*2); x1, y1];
            new_p2 = [particles(:, (pidx-1)*2+1:pidx*2); x2, y2];
            new_particles = [new_particles, new_p1, new_p2];
            
            new_p1 = [p_prob(:,pidx); p_prob(end,pidx)*this_prob1];
            new_p2 = [p_prob(:,pidx); p_prob(end,pidx)*this_prob2];
            % new_p1 = [p_prob(:,pidx); p_prob(end,pidx)+this_prob1];
            % new_p2 = [p_prob(:,pidx); p_prob(end,pidx)+this_prob2];
            new_p_prob = [new_p_prob, new_p1, new_p2];


            fprintf('  [%f,%f] p%d/%d: sol1=[%f,%f]=%f (prob=%f), sol2=[%f,%f]=%f (prob=%f)\n', real_trace(ti,:), pidx, size(p_prob,2), x1, y1, norm(real_trace(ti,:) - [x1,y1]), this_prob1, x2, y2, norm(real_trace(ti,:) - [x2,y2]), this_prob2);
            

            % clf(fh);
            % plot(spk_pos(1), spk_pos(2), 'bo'); 
            % hold on;
            % plot(real_trace(1:ti,1), real_trace(1:ti,2), '-b.');
            % hold on;
            % plot(prev_loc(1), prev_loc(2), 'go');
            % hold on;
            % plot([prev_loc(1) prev_loc(1)+vec2spk(ti-1,1)], [prev_loc(2) prev_loc(2)+vec2spk(ti-1,2)], '-y');
            % hold on;
            
            % hold on;

            % set(gca, 'XLim', [-0.1 0.5]);
            % set(gca, 'YLim', [-0.1 0.5]);

            % waitforbuttonpress
        end

        %% normalize probability
        % new_p_prob(end, :) = new_p_prob(end, :) / sum(new_p_prob(end, :));
        new_p_prob(end, :)

        %% update particles and probability
        particles = new_particles;
        p_prob = new_p_prob;

        %% remove some particles
        idx = find(p_prob(end, :) > 0);
        particles = particles(:, sort([(idx-1)*2+1, idx*2]));
        p_prob = p_prob(:, sort(idx));

        if size(p_prob, 2) > 1000;
            fprintf('--- remove some particles ------------\n');
            [v, idx] = sort(p_prob(end, :), 'descend');
            idx = idx(1:500);
            % sort(idx)
            particles = particles(:, sort([(idx-1)*2+1, idx*2]));
            p_prob = p_prob(:, sort(idx));

            %% normalize probability
            % p_prob(end, :) = p_prob(end, :) / sum(p_prob(end, :));
            
            % particles(end, :)
            % p_prob(end, :)
            fprintf('--- done ------------\n');
        end

        %% choose the trajectory with highest probability
        % [v, idx] = max(p_prob(end, :));
        % est_trace = [est_trace; particles(end, (idx-1)*2+1:idx*2)];
    end

    %% choose the trajectory with highest probability
    [v, idx] = max(p_prob(end, :));
    est_trace = particles(:, (idx-1)*2+1:idx*2);

end
