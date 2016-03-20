function batch_process_rx_both
    font_size = 28;

    comb_methods = {'avg'};
    % exps = {'rx1', 'rx2', 'rx3', 'rx4', 'rx5', 'rx6', 'rx7', 'rx8', 'rx9', 'rx10', 'rx11', 'rx12', 'rx.2.1', 'rx.2.2', 'rx.2.3', 'rx.2.4', 'rx.2.5', 'rx.2.6', 'rx.2.7', 'rx.2.8', 'rx.2.9'};
    for expi = 1:22
        exps{expi} = sprintf('rx.3.%d', expi);
    end
    confs = [ones(1,22)*3];

    output_dir = './data/';

    for ci = 1:length(comb_methods)
        combine_method = char(comb_methods{ci});

        avg_errs = zeros(length(exps), 3);
        std_errs = zeros(length(exps), 3);
        for expi = 1:length(exps)
            exp_name = char(exps{expi});
            filename = sprintf('%s%s.%s.results.txt', output_dir, exp_name, combine_method);

            if exist(filename, 'file') == 2
                tmp = load(filename);
                avg_errs(expi, :) = tmp(1, :);
                std_errs(expi, :) = tmp(2, :);
            else
                [tmp_avg_errs, tmp_std_errs] = process_rx_both(exp_name, confs(expi), combine_method);
                avg_errs(expi, :) = tmp_avg_errs;
                std_errs(expi, :) = tmp_std_errs;

                dlmwrite(filename, [tmp_avg_errs; tmp_std_errs], 'delimiter', '\t');
            end
        end

        % avg_errs = [1 2 3; 4 5 6];
        % std_errs = [1 2 3; 4 5 6] * 0.1;

        fh = figure(21); clf;
        hb = bar(1:length(exps), avg_errs);
        hold on;

        for bi = 1:length(hb)
            xtik = hb(bi).XData;
            idx = bi - mean(1:length(hb));
            act_xtik = xtik + idx * 0.23;
            h = errorbar(act_xtik, avg_errs(:,bi), std_errs(:,bi), '-k', 'linestyle', 'none');
            set(h, 'LineWidth', 2);
        end
        ylim(max(0, get(gca, 'ylim')))

        set(gca, 'FontSize', font_size);
        legend('combine', 'PN', 'FMCW');
        xlabel('exp', 'FontSize', font_size);
        ylabel('error (cm)', 'FontSize', font_size);
        title(combine_method);
        hold on;
    end

end
