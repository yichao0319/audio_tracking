function show_progress(t, T, t0)
    persistent status;

    if t == t0
        status = zeros(1,11);
    else
        idx = floor(t/T*10)+1;
        if status(idx) == 0
            status(idx) = 1;
            fprintf('  %d%%\n', (idx-1)*10);
        end
    end
end
