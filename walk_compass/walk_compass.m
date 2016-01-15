%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%%
%% - Input:
%%   - filename [XArbitraryZVertical]
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
%%   walk_compass('0904.exp1.square')
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function walk_compass(filename)
    % addpath('../utils');
    
    %% --------------------
    %% DEBUG
    %% --------------------
    DEBUG0 = 0;
    DEBUG1 = 1;
    DEBUG2 = 1;  %% progress
    DEBUG3 = 1;  %% verbose
    DEBUG4 = 1;  %% results

    HEAL_STRIKE_THRESH = -0.1;
    


    %% --------------------
    %% Constant
    %% --------------------


    %% --------------------
    %% Variable
    %% --------------------
    % input_dir  = './raw/';
    input_dir  = './raw_combine_walk_sound/';
    output_dir = './mat/';

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
    if DEBUG2, fprintf('Read data: %s\n', [input_dir filename '.txt']); end

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
    fprintf('  Fs = %f, %f, %f, %f\n', Fs_mag, Fs_accl, Fs_gyro, Fs_motion);
    

    %% --------------------
    %% Sync data
    %% --------------------
    if DEBUG2, fprintf('Sync data\n'); end

    common_time = 1:1/Fs:min([data(5,end),data(12,end),data(16,end),data(20,end)]);

    data_sync(5,:) = common_time;
    data_sync(12,:) = common_time;
    data_sync(16,:) = common_time;
    data_sync(20,:) = common_time;
    [v, idx] = unique(data(5,:));
    for fi = [6:11]
        data_sync(fi,:) = interp1(data(5,idx), data(fi,idx), common_time);
    end
    [v, idx] = unique(data(12,:));
    for fi = [13:15]
        data_sync(fi,:) = interp1(data(12,idx), data(fi,idx), common_time);
    end
    [v, idx] = unique(data(16,:));
    for fi = [17:19]
        data_sync(fi,:) = interp1(data(16,idx), data(fi,idx), common_time);
    end
    [v, idx] = unique(data(20,:));
    for fi = [21:40]
        data_sync(fi,:) = interp1(data(20,idx), data(fi,idx), common_time);
    end
    
    data = data_sync;


    %% ----------------
    %% Accel, Gyro, Rotation
    %% ----------------
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(3,1,1);
    plot(data(12,:), data(13,:), '-r');
    hold on;
    plot(data(12,:), data(14,:), '-g');
    plot(data(12,:), data(15,:), '-b');
    grid();
    title('Accl')

    subplot(3,1,2);
    plot(data(16,:), data(17,:), '-r');
    hold on;
    plot(data(16,:), data(18,:), '-g');
    plot(data(16,:), data(19,:), '-b');
    grid();
    title('Gyro')

    subplot(3,1,3);
    plot(data(20,:), data(21,:), '-r');
    hold on;
    plot(data(20,:), data(22,:), '-g');
    plot(data(20,:), data(23,:), '-b');
    grid();
    title('Rotate')


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
        end
    end

    
    for ti = 1:size(dcm, 3)
        delta = - sum(accl(:,ti) .* gravity(:,ti)) / sum(gravity(:,ti).^2);
        accl_proj(:,ti) = accl(:,ti) + gravity(:,ti) * delta;
    end
    

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(3,1,1);
    plot(data(20,:), accl(1,:), '-r');
    hold on;
    plot(data(20,:), data(27,:), '-g');
    plot(data(12,:), data(13,:), '-b');
    plot(data(20,:), accl_proj(1,:), '-c');
    legend('rotate', 'DM', 'raw', 'proj');
    title('x')
    
    subplot(3,1,2);
    plot(data(20,:), accl(2,:), '-r');
    hold on;
    plot(data(20,:), data(28,:), '-g');
    plot(data(12,:), data(14,:), '-b');
    plot(data(20,:), accl_proj(2,:), '-c');
    legend('rotate', 'DM', 'raw', 'proj');
    title('y')

    subplot(3,1,3);
    plot(data(20,:), accl(3,:), '-r');
    hold on;
    plot(data(20,:), data(29,:), '-g');
    plot(data(12,:), data(15,:), '-b');
    plot(data(20,:), accl_proj(3,:), '-c');
    legend('rotate', 'DM', 'raw', 'proj');
    title('z')

    
    %% ----------------
    %% Direction using DM gyroscope
    %% ----------------
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    plot(data(20,:), direction_dif, '-r');
    title('direction change')


    %% ----------------
    %% Direction using Walk Compass
    %% ----------------
    heel_det = sum(accl_proj);
    idx = find(heel_det < HEAL_STRIKE_THRESH);
    
    accel_direction_dif(1) = 0;
    for ti = 2:size(accl_proj,2)
        a = accl_proj(:,1) / norm(accl_proj(:,1));
        b = accl_proj(:,ti) / norm(accl_proj(:,ti));
        accel_direction_dif(ti) = atan2(norm(cross(a,b)),dot(a,b));
    end


    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(3,1,1)
    plot(data(20,:), heel_det, '-r.');
    hold on;
    plot(data(20,:), direction_dif, '-b');
    plot(data(20,idx), direction_dif(idx), 'y.');
    % plot(data(20,:), accel_direction_dif, '-g');
    % plot(data(5,:), data(9,:)/180*pi, '-g');
    title('heel strike detection')
    grid()

    filtered_dir = direction_dif(idx(1));
    filtered_accel_dir = accel_direction_dif(idx(1));
    alpha = 0.8;
    for ti = 2:length(idx)
        filtered_dir(ti) = alpha*filtered_dir(ti-1) + (1-alpha)*direction_dif(idx(ti));
        filtered_accel_dir(ti) = alpha*filtered_accel_dir(ti-1) + (1-alpha)*accel_direction_dif(idx(ti));
    end

    subplot(3,1,2)
    plot(data(20,idx), direction_dif(idx), 'r.');
    hold on;
    plot(data(20,idx), filtered_dir, '-b.');
    % plot(data(5,idx), data(9,idx)/180*pi, '-g.');
    set(gca, 'XLim', [data(20,1) data(20,end)]);
    title('direction')
    grid()

    subplot(3,1,3)
    plot(data(20,idx), accel_direction_dif(idx), 'r.');
    hold on;
    plot(data(20,idx), filtered_accel_dir, '-b.');
    title('direction accel')
    grid()


end