%% gen_pure_tone
function gen_pure_tone()
    
    output_dir = './tx_sound/';
    freq = 8000;
    total_len = 600;
    sound_len = 0.01;
    sound_itv = 0.5;
    amp = 1;
    Fs = 44100;


    %% generate base sound
    fprintf('generate base sound\n');
    signal = zeros(Fs*total_len, 2);
    ts = [linspace(0, sound_len, Fs*sound_len)];
    sound_sample = sin(2*pi*freq*ts)';

    %% concatenate base sound
    fprintf('concatenate base sound\n');
    n_frame = floor(total_len / (sound_len + sound_itv));
    for fi = 1:n_frame
        std_idx = int32((fi-1) * (sound_len + sound_itv) * Fs + 1);
        end_idx = int32((fi-1) * (sound_len + sound_itv) * Fs + sound_len * Fs);
        signal(std_idx:end_idx, 2) = sound_sample;
    end

    % signal = zeros(Fs*total_len, 2);
    % ts = [linspace(0, total_len, Fs*total_len)];
    % signal(:, 2) = amp * sin(2*pi*freq*ts)';

    % n_frame = floor(total_len / (sound_len + sound_itv));
    % for fi = 1:n_frame
    %     std_idx = int32((fi-1) * (sound_len + sound_itv) * Fs + sound_len * Fs + 1);
    %     end_idx = int32(fi * (sound_len + sound_itv) * Fs);
    %     fprintf('  frame %d: %d-%d\n', fi, std_idx, end_idx);
    %     signal(std_idx:end_idx, 2) = signal(std_idx:end_idx, 2) / 100;

    % end

    fprintf('write to the file\n');
    audiowrite([output_dir num2str(freq) '.' num2str(sound_len) '.' num2str(sound_itv) '.wav'], signal, Fs);
