% Given a point and obstacle cell array (of ellipsotopes), determine if the
% point lies inside any of the obstacles
function check = is_in_obs(point, obs)
    check = false;
    for i = 1:length(obs)
        if obs{i}.contains(point)
            check = true;
            return
        end
    end
end