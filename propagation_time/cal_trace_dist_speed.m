%% dist: sound
%% speed: doppler
%% direction: accelerometer
function [est_trace] = cal_trace_dist_speed(spk_pos, time, dist, speed, direction, real_trace)
    vel_err = 0.000;
    dist_err = 0.001;

    %% GIVEN: init_pos, dist, speed

    %% get relative velocity
    vec2spk = spk_pos - real_trace(1,:);
    vec2spk = vec2spk / norm(vec2spk);
    rel_velocity = speed(1) * vec2spk;
    
    est_trace = [real_trace(1,:)];


    % fh = figure(5); clf(fh);

    for ti = 2:length(time)
        %% known: real_trace(:, ti-1), speed(:, ti-1), dist(ti)

        % prev_loc = real_trace(ti-1,:);  %% XXX: for now...
        prev_loc = est_trace(ti-1,:);


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
        dist2spk = dist2spk + dist_err;
        rel_velocity(ti-1,:) = rel_velocity(ti-1,:) .* (1 + vel_err * (rand(1,2)*2-1));
        % dist2spk = dist2spk * (1 + dist_err * (rand(1)*2-1));
        % rel_velocity(ti-1,:) = rel_velocity(ti-1,:) .* (1 + vel_err);
        % dist2spk = dist2spk * (1 + dist_err);
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
        
        if b^2 - 4*a*c < 0
            % error('no intersection');
            fprintf('no intersection\n');
            % x1 = prev_loc(1)+this_direction(1)*(time(ti)-time(ti-1));
            % x2 = prev_loc(1)+this_direction(1)*(time(ti)-time(ti-1));
            % y1 = prev_loc(2)+this_direction(2)*(time(ti)-time(ti-1));
            % y2 = prev_loc(2)+this_direction(2)*(time(ti)-time(ti-1));
            x1 = (- b + sqrt(b^2 - 4*a*c)) / 2 / a;
            y1 = slope * (x1 - proj_pt(1)) + proj_pt(2);
            x2 = (- b - sqrt(b^2 - 4*a*c)) / 2 / a;
            y2 = slope * (x2 - proj_pt(1)) + proj_pt(2);
            % waitforbuttonpress
        else
            x1 = (- b + sqrt(b^2 - 4*a*c)) / 2 / a;
            y1 = slope * (x1 - proj_pt(1)) + proj_pt(2);
            x2 = (- b - sqrt(b^2 - 4*a*c)) / 2 / a;
            y2 = slope * (x2 - proj_pt(1)) + proj_pt(2);
        end
        
        % fprintf('  [%f,%f]: sol1=[%f,%f]=%f, sol2=[%f,%f]=%f\n', real_trace(ti,:), x1, y1, norm(real_trace(ti,:) - [x1,y1]), x2, y2, norm(real_trace(ti,:) - [x2,y2]));
        fprintf('  [%.0f,%.0f] (d=%.2f): sol1=[%.0f,%.0f]=%.2f (d=%.2f), sol2=[%.0f,%.0f]=%.2f (d=%.2f)\n', real_trace(ti,:) * 100, dist2spk * 100, x1 * 100, y1 * 100, norm(real_trace(ti,:) - [x1,y1]) * 100, norm([x1,y1]-spk_pos) * 100, x2 * 100, y2 * 100, norm(real_trace(ti,:) - [x2,y2]) * 100, norm([x2,y2]-spk_pos) * 100);
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
        dir1 = [x1,y1] - prev_loc; dir1 = dir1 / norm(dir1);
        dir2 = [x2,y2] - prev_loc; dir2 = dir2 / norm(dir2);
        if sum(dir1 .* this_direction) > sum(dir2 .* this_direction)
            est_trace(ti,:) = [x1,y1];
        else
            est_trace(ti,:) = [x2,y2]; 
        end


        % clf(fh);
        % plot(spk_pos(1), spk_pos(2), 'bo'); 
        % hold on;
        % plot(real_trace(1:ti,1), real_trace(1:ti,2), '-b.');
        % hold on;
        % plot(prev_loc(1), prev_loc(2), 'go');
        % hold on;
        % plot([prev_loc(1) prev_loc(1)+vec2spk(ti-1,1)], [prev_loc(2) prev_loc(2)+vec2spk(ti-1,2)], '-y');
        % hold on;
        % set(gca, 'XLim', [-0.1 0.5]);
        % set(gca, 'YLim', [-0.1 0.5]);

        % waitforbuttonpress
    end
end
