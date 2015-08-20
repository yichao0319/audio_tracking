%% dist [sound]: distance to speaker
%% speed [doppler]: relative speed w.r.t. speaker
%% direction [accel]: moving direction
function [est_trace] = cal_trace_dist_speed_opt(spk_pos, time, dist, speed, direction, real_trace)
    
    DEBUG1 = 0;
    DEBUG2 = 1;
    DEBUG3 = 0;

    v_err = 0.000;
    dist_err = 0.00;

    alpha1 = 10;
    alpha2 = 1;
    alpha3 = 0;
    alpha4 = 1;
    alpha5 = 5;

    epsilon = 0.00000;

    function f = one_snapshot_obj(xx)

        %% distance to speaker
        term1 = abs(norm(spk_pos-xx) - d2);

        %% velocity to speaker
        term2 = 0;
        if(abs(d1_err) > epsilon)
            term2 = abs(norm(dot(xx-l1, D1/d1_err) * (D1/d1_err)) - norm(Ms1_err));
        else
            term2 = norm(Ms1_err);
        end

        %% moving direction
        term3 = 0;
        if(norm(xx-l1) > epsilon)
            term3 = norm((xx-l1)/norm(xx-l1) - Dir1);
        else
            term3 = norm(Dir1);
        end

        %% moving speed
        speed_thresh = 0.5;  %% m/s
        term4 = 0;
        this_speed = norm((xx-l1) / t);
        if this_speed > speed_thresh
            term4 = 1;
        else
            term4 = 0;
        end

        %% moving distance
        term5 = abs(norm(xx-l1) - norm(M1));


        f = alpha1 * term1 + ...
            alpha2 * term2 + ...
            alpha3 * term3 + ...
            alpha4 * term4 + ...
            alpha5 * term5;
    end


    est_trace = [real_trace(1,:)];
    for ti = 2:length(time)
        %% known: real_trace(:, ti-1), speed(:, ti-1), dist(ti)

        % l1 = real_trace(ti-1,:);       %% XXX: for now...
        l1 = est_trace(ti-1,:);        %% location at ti-1
        D1 = spk_pos - l1;             %% direction to speaker at ti-1
        d1 = dist(ti-1);               %% [Audio] distance to speaker at ti-1
        d2 = dist(ti);                 %% [Audio] distance to speaker at ti
        Dir1 = direction(ti-1, :);     %% [Accel] moving direction at ti-1
        Ms1 = speed(ti-1)*D1/norm(D1); %% [Audio] relative velocity to speaker at ti-1
        t = time(ti) - time(ti-1);
        % M1 = (real_trace(ti,:)-real_trace(ti-1,:)) / t; %% ground truth -- shouldn't be used
        M1 = norm(real_trace(ti,:)-real_trace(ti-1,:));  %% [Accel] movement at ti-1

        
        %% ---------------------
        %% inject error
        Ms1_err = Ms1 + v_err * (rand(1,2)*2-1);
        d1_err = d1 + dist_err;
        % Ms1_err = Ms1 .* (1 + v_err * (rand(1,2)*2-1));
        % d1_err = d1 * (1 + dist_err * (rand(1)*2-1));
        % Ms1_err = Ms1_err .* (1 + v_err);
        % d1_err = d1 * (1 + dist_err);
        %% ---------------------


        if 1
            % init_l2 = real_trace(ti, :);
            % init_l2 = real_trace(ti, :) + [0.1, 0.1];
            init_l2 = l1;
            % init_l2 = l1 + speed(ti-1)*Dir1*t;
            % init_l2 = l1 + norm(M1)*Dir1*t;
            options = optimoptions(@fminunc,'MaxFunEvals', 3000);
            [l2, fval] = fminunc(@one_snapshot_obj, init_l2, options);
            x1 = l2(1);
            y1 = l2(2);
        else
            gran = 0.01;
            best_x1 = -0.05;
            best_y1 = -0.05;
            min_obj = Inf;
            for x1 = -0.01:gran:0.27
                for y1 = -0.01:gran:0.21
                    obj = one_snapshot_obj([x1,y1]);
                    if obj < min_obj
                        best_x1 = x1;
                        best_y1 = y1;
                        min_obj = obj;
                    end
                end
            end
            x1 = best_x1;
            y1 = best_y1;
        end
        est_trace(ti,:) = [x1,y1];
        

        fprintf('  [%.0f,%.0f] (d=%.2f): sol1=[%.0f,%.0f]=%.2f (d=%.2f)\n', real_trace(ti,:) * 100, d2 * 100, x1 * 100, y1 * 100, norm(real_trace(ti,:) - [x1,y1]) * 100, norm([x1,y1]-spk_pos) * 100);
        
        if 0
            fprintf('===================================\n');
            ti
            l2
            fval
            abs(norm(spk_pos-l2) - d2)
            abs(norm(dot(l2-l1, D1/d1_err) * (D1/d1_err)) - norm(Ms1_err))
            norm((l2-l1)/norm(l2-l1) - Dir1)

            fprintf('----------------------\n');

            l2_gt = real_trace(ti, :)
            one_snapshot_obj(l2_gt)
            abs(norm(spk_pos-l2_gt) - d2)
            abs(norm(dot(l2_gt-l1, D1/d1_err) * (D1/d1_err)) - norm(Ms1_err))
            norm((l2_gt-l1)/norm(l2_gt-l1) - Dir1)
            
            fprintf('[%f,%f] - [%f,%f] = %f\n', l2_gt, l2, norm(l2_gt-l2));
            fprintf('===================================\n');
            waitforbuttonpress
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
