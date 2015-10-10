%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%%
%% - Input:
%%   - filename
%%     x. loggingTime
%%     x. loggingSample
%%     x. identifierForVendor
%%     5. locationHeadingTimestamp_since1970
%%     6. locationHeadingX
%%     7. locationHeadingY
%%     8. locationHeadingZ
%%     9. locationTrueHeading
%%     10. locationMagneticHeading
%%     11. locationHeadingAccuracy
%%     12. accelerometerTimestamp_sinceReboot
%%     13. accelerometerAccelerationX
%%     14. accelerometerAccelerationY
%%     15. accelerometerAccelerationZ
%%     16. gyroTimestamp_sinceReboot
%%     17. gyroRotationX
%%     18. gyroRotationY
%%     19. gyroRotationZ
%%     20. motionTimestamp_sinceReboot
%%     21. motionYaw
%%     22. motionRoll
%%     23. motionPitch
%%     24. motionRotationRateX
%%     25. motionRotationRateY
%%     26. motionRotationRateZ
%%     27. motionUserAccelerationX
%%     28. motionUserAccelerationY
%%     29. motionUserAccelerationZ
%%     30. motionQuaternionX
%%     31. motionQuaternionY
%%     32. motionQuaternionZ
%%     33. motionQuaternionW
%%     34. motionGravityX
%%     35. motionGravityY
%%     36. motionGravityZ
%%     37. motionMagneticFieldX
%%     38. motionMagneticFieldY
%%     39. motionMagneticFieldZ
%%     40. motionMagneticFieldCalibrationAccuracy
%%
%% - Output:
%%
%%
%% example:
%%   visualize_accel_walk('0904.exp1.square')
%%   visualize_accel_walk('0906.exp1.mao.walk.square.hand')
%%   visualize_accel_walk('0906.exp2.mao.walk.square.pocket')
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visualize_accel_walk(filename)
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
    input_dir  = './accel/';
    output_dir = '';
    fig_dir = './fig/';


    %% --------------------
    %% Variable
    %% --------------------
    fig_idx = 0;


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, filename = '0904.exp1.square'; end


    %% --------------------
    %% Main starts
    %% --------------------
    
    %% --------------------
    %% Read data
    %% --------------------
    data = load([input_dir filename '.txt'])';

    time_base_mag = data(5,1);
    data(5,:) = data(5,:) - time_base_mag;
    time_base_accl = data(12,1);
    data(12,:) = data(12,:) - time_base_accl;
    time_base_gyro = data(16,1);
    data(16,:) = data(16,:) - time_base_gyro;
    time_base_motion = data(20,1);
    data(20,:) = data(20,:) - time_base_motion;
    
    Fs_mag = size(data,2) / data(5,end);
    Fs_accl = size(data,2) / data(12,end);
    Fs_gyro = size(data,2) / data(16,end);
    Fs_motion = size(data,2) / data(20,end);
    Fs = floor(min([Fs_mag, Fs_accl, Fs_gyro, Fs_motion]));
    
    fprintf('  size = %dx%d\n', size(data));
    fprintf('  Fs = %d (%f, %f, %f, %f)\n', Fs, Fs_mag, Fs_accl, Fs_gyro, Fs_motion);


    % %% --------------------
    % %% Sync data
    % %% --------------------
    % if DEBUG2, fprintf('Sync data\n'); end

    % common_time = 1:1/Fs:min([data(5,end),data(12,end),data(16,end),data(20,end)]);

    % data_sync(5,:) = common_time;
    % data_sync(12,:) = common_time;
    % data_sync(16,:) = common_time;
    % data_sync(20,:) = common_time;
    % [v, idx] = unique(data(5,:));
    % for fi = [6:11]
    %     data_sync(fi,:) = interp1(data(5,idx), data(fi,idx), common_time);
    % end
    % [v, idx] = unique(data(12,:));
    % for fi = [13:15]
    %     data_sync(fi,:) = interp1(data(12,idx), data(fi,idx), common_time);
    % end
    % [v, idx] = unique(data(16,:));
    % for fi = [17:19]
    %     data_sync(fi,:) = interp1(data(16,idx), data(fi,idx), common_time);
    % end
    % [v, idx] = unique(data(20,:));
    % for fi = [21:40]
    %     data_sync(fi,:) = interp1(data(20,idx), data(fi,idx), common_time);
    % end
    
    % data = data_sync;


    %% ----------------
    %% Accel and average
    %% ----------------
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(2,1,1);
    plot(data(12,:), mean(abs(data(13:15,:))), '-y.');
    hold on;
    plot(data(12,:), data(13,:), '-r');
    plot(data(12,:), data(14,:), '-g');
    plot(data(12,:), data(15,:), '-b.');
    grid();
    title('Accl raw')

    subplot(2,1,2);
    plot(data(20,:), mean(abs(data(27:29,:))), '-y.');
    hold on;
    plot(data(20,:), data(27,:), '-r');
    plot(data(20,:), data(28,:), '-g');
    plot(data(20,:), data(29,:), '-b');
    grid();
    title('Accl DC')


    % %% --------------------
    % %% Remove DC (zero frequency component)
    % %% --------------------
    % if DEBUG2, fprintf('Remove DC\n'); end

    % [B,A] = butter(10, 0.15/(Fs/2), 'high');
    % for fi = [27:29]
    %     tmp = filter(B,A,data(fi,:));
    %     data(fi,:) = tmp;
    %     % data(fi,:) = data(fi,:) - mean(data(fi,:));
    % end

    %% --------------------
    %% Find Periodic Portion and Steps
    %% --------------------
    % if DEBUG2, fprintf('Find Periodic Portion and Steps\n'); end

    ts = data(28,:);
    ts = smoothn(ts, 500);
    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % plot(ts);
    % return
    ts_time = data(20,:);
    ts_len = length(ts);
    [step_idx, periodic_idx, autocorr, fig_idx] = get_step_idx(ts, [0.5,1]*Fs, fig_idx);
    ts_duration = ts_time(periodic_idx(2)) - ts_time(periodic_idx(1));
    fprintf('  periodic portion duration = %f\n', ts_duration);
    
    % return
    
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
    print(fh, '-dpsc', [fig_dir filename '.xcorr.eps']);


    %% ----------------
    %% rotate accel to a reference point
    %% ----------------
    %% ZXY: yaw, pitch, roll
    %% 21. motionYaw
    %% 22. motionRoll
    %% 23. motionPitch
    dcm = angle2dcm(data(21,:)', data(23,:)', data(22,:)');
    direction_dif(1) = 0;
    for ti = 1:size(dcm, 3)
        accl(:,ti) = squeeze(dcm(:,:,ti)) * data(27:29,ti);

        gravity(:,ti) = squeeze(dcm(:,:,ti)) * data(34:36,ti);

        direction(:,ti) = squeeze(dcm(:,:,ti)) * [1; 0; 0];
        if ti > 1
            a = direction(:,1) / norm(direction(:,1));
            b = direction(:,ti) / norm(direction(:,ti));
            direction_dif(ti) = atan2(norm(cross(a,b)),dot(a,b));
            % direction_dif(ti) = mod(atan2(b(2)-a(2),b(1)-a(1)),2*pi);
        end
    end

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(2,1,1);
    plot(data(20,:), data(27,:), '-r');
    hold on;
    plot(data(20,:), data(28,:), '-g');
    plot(data(20,:), data(29,:), '-b');
    grid();
    title('Accl DM')

    subplot(2,1,2);
    plot(data(20,:), direction_dif, '-b');
    hold on;
    % plot(data(5,:), data(10,:)/180*pi, '-r');
    grid();
    title('Gyro direction')

    print(fh, '-dpsc', [fig_dir filename '.direction.eps']);

    %% get avg direction in each period
    % periods = [0 5; 5 10; 10 15; 17 24; 28 35; 35 40; 45 50; 50 55];
    % for ti = 1:size(periods,1)
    %     idx = find(data(20,:) >= periods(ti, 1) & data(20,:) < periods(ti, 2));
    %     avg_dir = sum(direction(:,idx), 2);
    %     avg_dir = avg_dir / norm(avg_dir)
    % end


end