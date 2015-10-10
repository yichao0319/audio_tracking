fig_dir = './fig/';
fig_idx = 0;
ei = 0;
err = {[1]};

colors   = {'r', 'b', [0 0.8 0], 'm', [1 0.85 0], [0 0 0.47], [0.45 0.17 0.48], 'k'};
lines    = {'-', '--', '-.', ':'};
markers  = {'+', 'o', '*', '.', 'x', 's', 'd', '^', '>', '<', 'p', 'h'};
font_size = 28;


%%=============================================
% ei = ei + 1;
% err{ei} = [est_fix_dist_pn('0904.exp1.fix.pn2047.loc1', 8000, 2047, sqrt(2.8^2+1^2), 1)];
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0904.exp1.fix.pn2047.loc2', 8000, 2047, 1, 1);
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0904.exp1.fix.pn2047.loc3', 8000, 2047, sqrt(3.1^2+1^2), 1);
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0904.exp1.fix.pn2047.loc4', 8000, 2047, sqrt(3.1^2+2^2), 1);
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0904.exp1.fix.pn2047.loc5', 8000, 2047, 2, 1);
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0904.exp1.fix.pn2047.loc6', 8000, 2047, sqrt(2.8^2+2^2), 1);
% mean_err(ei) = mean(err{ei});

%%=============================================
% ei = ei + 1;
% err{ei} = est_fix_dist_sinc('0904.exp2.fix.sinc.loc1', 8000, sqrt(2.8^2+1^2), 1)
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_sinc('0904.exp2.fix.sinc.loc2', 8000, 1, 1)
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_sinc('0904.exp2.fix.sinc.loc3', 8000, sqrt(3.1^2+1^2), 1)
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_sinc('0904.exp2.fix.sinc.loc4', 8000, sqrt(3.1^2+2^2), 1)
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_sinc('0904.exp2.fix.sinc.loc5', 8000, 2, 1)
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_sinc('0904.exp2.fix.sinc.loc6', 8000, sqrt(2.8^2+2^2), 1)
% mean_err(ei) = mean(err{ei});

%%=============================================
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp4.fix.pn2047.loc1', 8000, 2047, sqrt(2.8^2+1^2), 1);
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp4.fix.pn2047.loc2', 8000, 2047, 1, 1);
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp4.fix.pn2047.loc3', 8000, 2047, sqrt(3.1^2+1^2), 1);
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp4.fix.pn2047.loc4', 8000, 2047, sqrt(3.1^2+2^2), 1);
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp4.fix.pn2047.loc5', 8000, 2047, 2, 1);
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp4.fix.pn2047.loc6', 8000, 2047, sqrt(2.8^2+2^2), 1);
% mean_err(ei) = mean(err{ei});


%%=============================================
ei = ei + 1;
err{ei} = est_fix_dist_pn('0906.exp2.3m.wall.degree0', 8000, 2047, 1, 1)
mean_err(ei) = mean(err{ei});
ei = ei + 1;
err{ei} = est_fix_dist_pn('0906.exp2.3m.wall.degree45', 8000, 2047, 1, 1)
mean_err(ei) = mean(err{ei});
ei = ei + 1;
err{ei} = est_fix_dist_pn('0906.exp2.3m.wall.degree90', 8000, 2047, 1, 1)
mean_err(ei) = mean(err{ei});
ei = ei + 1;
err{ei} = est_fix_dist_pn('0906.exp2.3m.wall.degree135', 8000, 2047, 1, 1)
mean_err(ei) = mean(err{ei});
ei = ei + 1;
err{ei} = est_fix_dist_pn('0906.exp2.3m.wall.degree180', 8000, 2047, 1, 1)
mean_err(ei) = mean(err{ei});
ei = ei + 1;
err{ei} = est_fix_dist_pn('0906.exp2.3m.wall.degree225', 8000, 2047, 1, 1)
mean_err(ei) = mean(err{ei});
ei = ei + 1;
err{ei} = est_fix_dist_pn('0906.exp2.3m.wall.degree270', 8000, 2047, 1, 1)
mean_err(ei) = mean(err{ei});
ei = ei + 1;
err{ei} = est_fix_dist_pn('0906.exp2.3m.wall.degree315', 8000, 2047, 1, 1)
mean_err(ei) = mean(err{ei});

%%=============================================
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp1.1m.degree0', 8000, 2047, 1, 1)
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp1.1m.degree45', 8000, 2047, 1, 1)
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp1.1m.degree90', 8000, 2047, 1, 1)
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp1.1m.degree135', 8000, 2047, 1, 1)
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp1.1m.degree180', 8000, 2047, 1, 1)
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp1.1m.degree225', 8000, 2047, 1, 1)
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp1.1m.degree270', 8000, 2047, 1, 1)
% mean_err(ei) = mean(err{ei});
% ei = ei + 1;
% err{ei} = est_fix_dist_pn('0906.exp1.1m.degree315', 8000, 2047, 1, 1)
% mean_err(ei) = mean(err{ei});


fh = figure(1); clf;
for ei = 1:length(err)
    [f,x] = ecdf(err{ei});

    lh = plot(x,f);
    set(lh, 'Color', colors{mod(ei-1,length(colors))+1});
    set(lh, 'LineStyle', char(lines{mod(ei-1,length(lines))+1}));
    set(lh, 'LineWidth', 4);
    set(lh, 'marker', markers{mod(ei-1,length(markers))+1});
    hold on;
    % return
    % legends{ei} = ['loc' num2str(ei)];
    legends{ei} = ['angle' num2str((ei-1)*45)];
end
% set(gca, 'XScale', 'log');
legend(legends, 'Location','SouthEast');
% set(gca, 'XLim', [0 10]);
xlabel('Error (m)', 'FontSize', font_size)
ylabel('CDF', 'FontSize', font_size)
set(gca, 'FontSize', font_size);

print(fh, '-dpsc', [fig_dir 'sound.eps']);
