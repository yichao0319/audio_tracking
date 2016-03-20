function read_batch_results
    font_size = 28;

    comb_methods = {'avg', 'mrc_var', 'mrc_spk', 'mrc_spk2'};
    exps = {'rx1', 'rx2', 'rx3', 'rx4', 'rx5', 'rx6', 'rx7', 'rx8', 'rx9', 'rx10', 'rx11', 'rx12', 'rx.2.1', 'rx.2.2', 'rx.2.3', 'rx.2.4', 'rx.2.5', 'rx.2.6', 'rx.2.7', 'rx.2.8', 'rx.2.9'};
    % exps = {'rx1', 'rx2'};
    confs = [ones(1,12)*1, ones(1,9)*2];

    input_dir = './data/';


    n_methods = 2+length(comb_methods);
    ret = ones(length(exps), n_methods, 2) * -1;

    for ci = 1:length(comb_methods)
        combine_method = char(comb_methods{ci});

        for expi = 1:length(exps)
            exp_name = char(exps{expi});
            filename = sprintf('%s%s.%s.results.txt', input_dir, exp_name, combine_method);

            if exist(filename, 'file') == 2
                tmp = load(filename);
                ret(expi, 2+ci,  1) = tmp(1, 1);
                ret(expi, 2+ci,  2) = tmp(2, 1);

                if ci == 1
                    ret(expi, [1,2], 1) = tmp(1, [2,3]);
                    ret(expi, [1,2], 2) = tmp(2, [2,3]);
                end
            end
        end
    end

    dlmwrite(sprintf('%ssummary.avg.txt', input_dir), [squeeze(ret(:,:,1))], 'delimiter', '\t');
    dlmwrite(sprintf('%ssummary.std.txt', input_dir), [squeeze(ret(:,:,2))], 'delimiter', '\t');


end
