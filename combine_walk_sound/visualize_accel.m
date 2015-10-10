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
    %% Accel and average
    %% ----------------
    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    subplot(2,1,1);
    plot(data(12,:), mean(abs(data(13:15,:))), '-y.');
    hold on;
    plot(data(12,:), data(13,:), '-r');
    plot(data(12,:), data(14,:), '-g');
    plot(data(12,:), data(15,:), '-b');
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
end