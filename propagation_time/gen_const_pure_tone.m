%% gen_pure_tone
function gen_const_pure_tone()

    output_dir = './tx_sound/';
    freqs = [0 8000];
    t = 1*300;
    amp = 1;
    % t = 1*60;
    Fs = 44100;

    signals = zeros(Fs*t, length(freqs));

    for fi = 1:length(freqs)
        f = freqs(fi);
        ts = [linspace(0, t, Fs*t)];
        signals(:, fi) = sin(2*pi*f*ts)';
    end


    %% freq L[17k], R[18k]
    idxl = find(freqs == 0);
    idxr = find(freqs == 8000);
    out = amp*signals(:, [idxl, idxr]);
    audiowrite([output_dir '8000.wav'], out, Fs);

