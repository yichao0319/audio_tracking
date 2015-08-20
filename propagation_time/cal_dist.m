function [dists] = cal_dist(pos1, pos2)
    dists = sqrt(abs((pos1(:,1) - pos2(:,1)) .^ 2 + (pos1(:,2) - pos2(:,2)) .^ 2));
end
