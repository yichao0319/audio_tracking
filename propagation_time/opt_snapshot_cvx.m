%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%%
%% - Input:
%%
%%
%% - Output:
%%
%%
%% example:
%%
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function opt_snapshot_cvx()
    
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


    %% --------------------
    %% Variable
    %% --------------------
    input_dir  = '';
    output_dir = '';


    %% --------------------
    %% Check input
    %% --------------------
    % if nargin < 1, arg = 1; end
    % if nargin < 1, arg = 1; end


    %% --------------------
    %% Main starts
    %% --------------------
    alpha1 = 1;
    alpha2 = 1;
    alpha3 = 1;

    x1 = [0, 0];
    x2_gt = [0.1, 0.1];
    spk_pos = [0.4, 0.1];

    D1 = spk_pos - x1;
    D2 = spk_pos - x2_gt(2);
    d1 = norm(D1);
    d2 = norm(D2);
    M1 = x2_gt - x1;
    Mspk1 = (dot(M1, D1) / d1) * (D1 / d1);

    % Dacc1 = M1 / norm(M1);
    Dacc1 = M1;


    %% ==================================================================
    %% algorithm :)
    prev_loc = x1;
    dist2spk = d2;
    proj_move_dist = Mspk1;
    proj_pt = prev_loc + proj_move_dist;
    possible_line_vec = [-proj_move_dist(2), proj_move_dist(1)];
    slope = possible_line_vec(2) / possible_line_vec(1);

    a = (1 + slope^2);
    w = - slope * proj_pt(1) + proj_pt(2) - spk_pos(2);
    b = 2 * (slope * w - spk_pos(1));
    c = (spk_pos(1)^2 + w^2 - dist2spk^2);

    cal_x1 = (- b + sqrt(b^2 - 4*a*c)) / 2 / a;
    cal_y1 = slope * (cal_x1 - proj_pt(1)) + proj_pt(2);
    cal_x2 = (- b - sqrt(b^2 - 4*a*c)) / 2 / a;
    cal_y2 = slope * (cal_x2 - proj_pt(1)) + proj_pt(2);
    %% ==================================================================


    %% ==================================================================
    %% fmincon
    %%   inject error
    % x1 = x1 + 0.1;

    function f = one_snapshot_obj(xx)
        f = alpha1 * abs(norm(spk_pos - xx) - d2) + ...
            alpha2 * abs(norm(dot(xx-x1, D1/d1) * (D1/d1)) - norm(Mspk1)) + ...
            alpha3 * norm((xx-x1) - Dacc1);
    end

    x2_0 = x1;
    [x2,fval] = fminunc(@one_snapshot_obj,x2_0);
    x2
    x2_gt
    fval
    one_snapshot_obj(x2_gt)
    return
    %% ==================================================================

    cvx_begin
        cvx_quiet(true)
        % variable x2(2)
        % variable u2(2)
        % variable v2(2)
        variable x2(1,3)
        % size(x2)
        % minimize( alpha1 * (norm(spk_pos - x2(3:4)) - d2) + ...
        %           alpha2 * (norm(dot(x2(5:6)-x1, D1/d1) * (D1/d1)) - norm(Mspk1)) );
        %         % - alpha3 * dot((x2(1:2)-x1)/norm(x2(1:2)-x1), Dacc1) );
        % minimize( alpha1 * x2(3) + ...
        %           alpha2 * x2(4) );
        % minimize( alpha1 * norm(((spk_pos(1)-x2(1)) + (spk_pos(2)-x2(2))) - d2^2, 1));
        minimize( alpha1 * x2(3) )
        subject to
            % x2(1:2) <= x2(3:4)
            % -x2(1:2) <= x2(3:4)
            % x2(1:2) <= x2(5:6)
            % -x2(1:2) <= x2(5:6)
            % norm(spk_pos - x2(1:2), 2) - d2 - x2(3) <= 0
            % -norm(spk_pos - x2(1:2), 2) - d2 - x2(3) <= 0
            % norm(dot(x2(1:2)-x1, D1/d1) * (D1/d1), 2) - norm(Mspk1, 2) - x2(4) <= 0
            % -norm(dot(x2(1:2)-x1, D1/d1) * (D1/d1), 2) - norm(Mspk1, 2) - x2(4) <= 0
            norm(spk_pos(1)-x2(1),1) - norm(spk_pos(1)-x2_gt(1),1) < x2(3)
            % -(norm(spk_pos(1)-x2(1),1) - norm(spk_pos(1)-x2_gt(1),1)) < x2(3)
    cvx_end
    x2
    norm(spk_pos(1)-x2(1),1)

    (norm(spk_pos(1)-x2(1),1) - norm(spk_pos(1)-x2_gt(1),1))
    x2_gt
    [cal_x1, cal_y1]
    [cal_x2, cal_y2]
end

