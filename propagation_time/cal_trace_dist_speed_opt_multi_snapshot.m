%% dist [sound]: distance to speaker
%% speed [doppler]: relative speed w.r.t. speaker
%% direction [accel]: moving direction
function [est_traces] = cal_trace_dist_speed_opt_multi_snapshot(spk_pos, time, dist, speed, direction, real_trace)
    
    DEBUG1 = 0;
    DEBUG2 = 1;
    DEBUG3 = 0;

    v_err = 0.000;
    dist_err = 0.000;

    alpha1 = 1;
    alpha2 = 1;
    % alpha2 = 1;
    alpha3 = 1;
    % alpha4 = 0.000001;
    alpha4 = 0.000001;
    ns     = 1;   %% number of snapshots
    term1 = 0;
    term2 = 0;
    term3 = 0;
    term4 = 0;

    epsilon = 0.00000;

    function f = multi_snapshots_obj(xx)

        %% distance to the speaker
        est_dist = xx - repmat(spk_pos, this_ns, 1);
        est_dist = sqrt(est_dist(:,1).^2 + est_dist(:,2).^2);
        term1 = norm(est_dist - di) / this_ns;

        %% displacement
        est_displace = xx - [l1; xx(1:this_ns-1, :)];
        est_displace = sqrt(est_displace(:,1).^2 + est_displace(:,2).^2);
        term2 = norm(est_displace - mi) / this_ns;

        %% displacement from known
        est_displace_known = xx - repmat(l1, this_ns, 1);
        est_displace_known = sqrt(est_displace_known(:,1).^2 + est_displace_known(:,2).^2);
        term3 = norm(est_displace_known - mk) / this_ns;

        %% moving direction
        est_dir = xx - [l1; xx(1:this_ns-1, :)];
        est_dir = est_dir ./ repmat(est_displace, 1, 2);
        idx = find(est_displace == 0);
        est_dir(idx, :) = 0;
        term4 = -mean(sum(est_dir .* Mi, 2));

        f = alpha1 * term1 + ...
            alpha2 * term2 + ...
            alpha3 * term3 + ...
            alpha4 * term4;
    end


    est_trace = [real_trace(1,:)];
    for ti = 2:ns:length(time)
    % for ti = 2:1:length(time)
        this_ns = min(ns, length(dist)-ti+1);
        %% known: real_trace(:, ti-1), speed(:, ti-1), dist(ti)

        % l1 = real_trace(ti-1,:);       %% XXX: for now...
        l1 = est_trace(ti-1,:);        %% location at ti-1
        di = dist(ti-1:ti+this_ns-2);         %% [Audio] distance to speaker at ti-1:ti+ns-2
        mi = real_trace(ti:ti+this_ns-1,:)-real_trace(ti-1:ti+this_ns-2,:);
        mi = sqrt(mi(:,1).^2 + mi(:,2).^2);   %% [Accel] displacement from prev to next
        mk = real_trace(ti:ti+this_ns-1,:)-repmat(real_trace(ti-1,:), this_ns, 1);
        mk = sqrt(mk(:,1).^2 + mk(:,2).^2);   %% [Accel] displacement from first known to all
        Mi = real_trace(ti:ti+this_ns-1,:)-real_trace(ti-1:ti+this_ns-2,:);
        Mi = Mi ./ repmat(mi, 1, 2);          %% [Accel] moving direction
        idx = find(mi == 0);
        Mi(idx, :) = 0;

        % figure(3)
        % subplot(3,1,1);
        % plot(di)
        % subplot(3,1,2);
        % plot(mi)
        % subplot(3,1,3);
        % plot(mk)

        %% ------------
        %% DEBUG
        tmpdist = sqrt((real_trace(ti-1:ti+this_ns-2,1)-spk_pos(1)).^2 + (real_trace(ti-1:ti+this_ns-2,2)-spk_pos(2)).^2);
        if di ~= tmpdist, error('wrong dist to spk'); end
        %% ------------

        

        
        %% ---------------------
        %% inject error
        % mk = mk .* (1 + 0.5*rand(size(mk)));
        % Mi = Mi .* (1 + 0.5*rand(size(Mi)));
        %% ---------------------

        init_l2 = repmat(l1, this_ns, 1);
        % init_l2 = l1 + speed(ti-1)*Dir1*t;
        % init_l2 = l1 + norm(M1)*Dir1*t;
        options = optimoptions(@fminunc,'MaxFunEvals', 3000);
        [l2, fval] = fminunc(@multi_snapshots_obj, init_l2, options);
        est_trace(ti:ti+this_ns-1,:) = l2;

        multi_snapshots_obj(l2);
        fprintf('  %d: %f, %f, %f, %f\n', ti, term1, term2, term3, term4);

        % if term2 ~= term3
        %     waitforbuttonpress
        % end

        % fprintf('  [%.0f,%.0f] (d=%.2f): sol1=[%.0f,%.0f]=%.2f (d=%.2f)\n', real_trace(ti,:) * 100, d2 * 100, x1 * 100, y1 * 100, norm(real_trace(ti,:) - [x1,y1]) * 100, norm([x1,y1]-spk_pos) * 100);
        
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

    est_traces = [est_trace];


    %% second iterations
    ns = 1;
    % for itr_i = 1:2
    %     est_trace = [est_trace est_trace];
    %     for ti = 2:ns:length(time)
    %     % for ti = 2:1:length(time)
    %         this_ns = min(ns, length(dist)-ti+1);
    %         %% known: real_trace(:, ti-1), speed(:, ti-1), dist(ti)

    %         % l1 = real_trace(ti-1,:);       %% XXX: for now...
    %         l1 = est_trace(ti-1,3:4);        %% location at ti-1
    %         di = dist(ti-1:ti+this_ns-2);         %% [Audio] distance to speaker at ti-1:ti+ns-2
    %         mi = real_trace(ti:ti+this_ns-1,:)-real_trace(ti-1:ti+this_ns-2,:);
    %         mi = sqrt(mi(:,1).^2 + mi(:,2).^2);   %% [Accel] displacement from prev to next
    %         mk = real_trace(ti:ti+this_ns-1,:)-repmat(real_trace(ti-1,:), this_ns, 1);
    %         mk = sqrt(mk(:,1).^2 + mk(:,2).^2);   %% [Accel] displacement from first known to all
    %         Mi = real_trace(ti:ti+this_ns-1,:)-real_trace(ti-1:ti+this_ns-2,:);
    %         Mi = Mi ./ repmat(mi, 1, 2);          %% [Accel] moving direction
    %         idx = find(mi == 0);
    %         Mi(idx, :) = 0;

    %         %% ------------
    %         %% DEBUG
    %         tmpdist = sqrt((real_trace(ti-1:ti+this_ns-2,1)-spk_pos(1)).^2 + (real_trace(ti-1:ti+this_ns-2,2)-spk_pos(2)).^2);
    %         if di ~= tmpdist, error('wrong dist to spk'); end
    %         %% ------------

            

            
    %         %% ---------------------
    %         %% inject error
    %         % mk = mk .* (1 + 0.5*rand(size(mk)));
    %         % Mi = Mi .* (1 + 0.5*rand(size(Mi)));
    %         %% ---------------------

    %         init_l2 = est_trace(ti:ti+this_ns-1, 1:2);
    %         % init_l2 = l1 + speed(ti-1)*Dir1*t;
    %         % init_l2 = l1 + norm(M1)*Dir1*t;
    %         options = optimoptions(@fminunc,'MaxFunEvals', 3000);
    %         [l2, fval] = fminunc(@multi_snapshots_obj, init_l2, options);
    %         est_trace(ti:ti+this_ns-1,3:4) = l2;

    %         multi_snapshots_obj(l2);
    %         fprintf('  %d: %f, %f, %f, %f\n', ti, term1, term2, term3, term4);


    %     end


    %     est_trace = est_trace(:, 3:4);
    %     est_traces = [est_traces est_trace];
    % end
end
