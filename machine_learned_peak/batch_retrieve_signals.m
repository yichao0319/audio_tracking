function batch_retrieve_signals

    for expi = 1:23
        exps{expi} = sprintf('rx.3.%d', expi);
    end
    confs = [ones(1,23)*3];


    for expi = 3:23
        retrieve_signals(exps{expi}, confs(expi));
    end

end
