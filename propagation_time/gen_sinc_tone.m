%% gen_pure_tone
function gen_sinc_tone()
    
    output_dir = './tx_sound/';
    freq = 8000;
    total_len = 60;
    frame_len = 0.1;
    amp = 2;
    Fs = 44100;


    %% generate base sound
    fprintf('generate base sound\n');
    signal = zeros(Fs*total_len, 2);

    Ts = 1/Fs;
    L = 1000;
    Tb = 1/20000;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate pulse signals for two channels
    % index=10*Fs;
    % Tb=1/10000;
    % interval=0.04;
    % K=Fs*interval;
    % round=1000;
    % T=round*interval;
    % N=Fs*T;
    % K=K*3;
    % signal_left=zeros(N,1);
    
    % for i=-L:L
    %     signal_left(index+i)=1/2*sinc(i*Ts/Tb)+1/4*sinc(i*Ts/Tb-1)+1/4*sinc(i*Ts/Tb+1);
    % end
    % index=10*Fs+2*K;
    % for i=-L:L
    %     signal_left(index+i)=1/2*sinc(i*Ts/Tb)+1/4*sinc(i*Ts/Tb-1)+1/4*sinc(i*Ts/Tb+1);
    % end
    % index=10*Fs+4*K;
    % for i=-L:L
    %     signal_left(index+i)=1/2*sinc(i*Ts/Tb)+1/4*sinc(i*Ts/Tb-1)+1/4*sinc(i*Ts/Tb+1);
    % end
    % % up-convert
    % freq=16000;
    % for i=1:N
    %     signal_left(i)=signal_left(i)*sin(2*pi*freq*i*Ts);
    % end

    % audiowrite([output_dir 'sinc.' num2str(freq) '.wav'], [zeros(length(signal_left),1), amp*signal_left], Fs);
    % return

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ts = [-L:L];
    sound_sample = 1/2*sinc(ts*Ts/Tb) + 1/4*sinc(ts*Ts/Tb-1) + 1/4*sinc(ts*Ts/Tb+1);
    % sound_sample = 1/4*sinc(ts*Ts/Tb-1);
    sound_sample = sound_sample';
    % up-convert
    sound_sample = sound_sample .* sin(2*pi*freq*Ts * [1:length(sound_sample)])';
    n_sound = length(sound_sample);


    %% concatenate base sound
    fprintf('embed base sound\n');
    n_frame = floor(total_len / frame_len);
    for fi = 1:n_frame
        fprintf('  frame %d\n', fi);
        std_idx = int32((fi-1) * frame_len * Fs) + 1;
        end_idx = int32((fi-1) * frame_len * Fs) + n_sound;
        signal(std_idx:end_idx, 2) = sound_sample;
    end

    fprintf('write to the file\n');
    audiowrite([output_dir 'sinc.' num2str(freq) '.' num2str(frame_len) '.wav'], amp*signal, Fs);
