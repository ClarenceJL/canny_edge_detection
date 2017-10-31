function linkedMap = edgeLink( M, J, theta )
%% edgeLink
%  Input:
%   M - original edge map
%   J - magnitute of derivative
%   theta - orientation of derivative
%  Output:
%   linkedMap - linked edge map
%
%  Jiani Li, 09/25/2016

[row,col] = size(J);   % size of image

%% define thresholds
edge_prop = 0.28;
edgeJ = M.*J;
sorted_J = sort(J(:));

threshold_h = sorted_J(ceil(numel(J)*(1-edge_prop)));
threshold_l = threshold_h * 0.4;


%% double thresholding
strong_map = (edgeJ >= threshold_h);
weak_map = (edgeJ < threshold_h) & (edgeJ >= threshold_l);
weakandstrong_map = (edgeJ >= threshold_l);

% [rstrong,cstrong] = find(strong_map);
% linkedMap = bwselect(weakandstrong_map, cstrong, rstrong, 8);

%% edgelink by hysteresis
gamma = mod(theta + pi/2,pi); % tangent direction map gamma, 0 <= gamma < pi

% This function provide 3 different strategies to link edges:
% method:
% 1 - Discretize tangent directon gamma to 4 directions; Link a pair of 
%     neighbor points (if they are weak edge) for each strong edge point.
% 2 - Divide [0,pi) to 4 intervals, when gamma falls in one interval,
%     consider neighbor points on both sides; Link two pairs of
%     neighbor points (if they are weak edge) for each strong point.
% 3 - For each strong point, consider all its 8 neighbor points regardless
%     of the edge direction. This method is theoretically the same as the 
%     MATLAB implementation for canny edge detection (or bwselect). 

method = 2; % method = 1,2,3; default method is 2
itr = 10000;
linkedMap  = strong_map; %initialization

if method == 1
    % gamma discretization
    gamma((gamma >= 0) & (gamma < pi/8)) = 1;
    gamma((gamma >= pi/8) & (gamma < 3*pi/8)) = 2;
    gamma((gamma >= 3*pi/8) & (gamma < 5*pi/8)) = 3;
    gamma((gamma >= 5*pi/8) & (gamma < 7*pi/8)) = 4;
    gamma((gamma >= 7*pi/8) & (gamma <= pi)) = 1;

    while itr > 0
        % gamma == 1
        strong1 = (gamma == 1) & strong_map;
        left_nb = false(row,col); left_nb(:,1:end-1) = strong1(:,2:end);
        right_nb = false(row,col); right_nb(:,2:end) = strong1(:,1:end-1);
        % gamma == 2
        strong2 = (gamma == 2) & strong_map;
        leftup_nb = false(row,col); leftup_nb(1:end-1,1:end-1) = strong2(2:end,2:end);
        rightdown_nb = false(row,col); rightdown_nb(2:end,2:end) = strong2(1:end-1,1:end-1);
        % gamma == 3
        strong3 = (gamma == 3) & strong_map;
        up_nb = false(row,col); up_nb(1:end-1,:) = strong3(2:end,:);
        down_nb = false(row,col); down_nb(2:end,:) = strong3(1:end-1,:);
        % gamma == 4
        strong4 = (gamma == 4) & strong_map;
        leftdown_nb = false(row,col); leftdown_nb(2:end,1:end-1) = strong4(1:end-1,2:end);
        rightup_nb = false(row,col); rightup_nb(1:end-1,2:end) = strong4(2:end,1:end-1);
        % 
        weak_nb = (left_nb|right_nb|up_nb|down_nb|leftup_nb|rightdown_nb|rightup_nb|leftdown_nb) & weak_map;

        linkedMap = linkedMap | weak_nb;    % add selected weak points to final edge map
        weak_nb_ind = find(weak_nb == 1);
        weak_map(weak_nb_ind) = 0;   % remove selected weak points from weak_map

        if isempty(weak_nb_ind)  % cannot link to other points, end iteration
            break;
        else
            strong_map = weak_nb;   % update strong_map to newly added points
        end

        itr = itr - 1;

    end

else if method == 2
         while itr > 0
            % divide strong edge set to four classes according to gradient
            % direction
            strong1 = (gamma >= 0) & (gamma < pi/4) & strong_map;
            strong2 = (gamma >= pi/4) & (gamma < pi/2) & strong_map;
            strong3 = (gamma >= pi/2) & (gamma < 3*pi/4) & strong_map;
            strong4 = (gamma >= 3*pi/4) & (gamma <= pi) & strong_map;
            % select neighbor points
            type1 = strong1 | strong4; 
            left_nb = false(row,col);left_nb(:,1:end-1) = type1(:,2:end);
            right_nb = false(row,col); right_nb(:,2:end) = type1(:,1:end-1);
            type2 = strong1 | strong2;
            leftup_nb = false(row,col); leftup_nb(1:end-1,1:end-1) = type2(2:end,2:end);
            rightdown_nb = false(row,col); rightdown_nb(2:end,2:end) = type2(1:end-1,1:end-1);
            type3 = strong2 | strong3;
            up_nb = false(row,col); up_nb(1:end-1,:) = type3(2:end,:);
            down_nb = false(row,col); down_nb(2:end,:) = type3(1:end-1,:);
            type4 = strong3 | strong4;
            leftdown_nb = false(row,col); leftdown_nb(2:end,1:end-1) = type4(1:end-1,2:end);
            rightup_nb = false(row,col); rightup_nb(1:end-1,2:end) = type4(2:end,1:end-1);
            % 
            weak_nb = (left_nb|right_nb|up_nb|down_nb|leftup_nb|rightdown_nb|rightup_nb|leftdown_nb) & weak_map;

            linkedMap = linkedMap | weak_nb;    % add selected weak points to final edge map
            weak_nb_ind = find(weak_nb == 1);
            weak_map(weak_nb_ind) = 0;   % remove selected weak points from weak_map

            if isempty(weak_nb_ind)  % cannot link to other points, end iteration
                break;
            else
                strong_map = weak_nb;   % update strong_map to newly added points
            end

            itr = itr - 1;

         end
    else % method == 3
        while itr > 0
            left_nb = false(row,col); left_nb(:,1:end-1) = strong_map(:,2:end);
            right_nb = false(row,col); right_nb(:,2:end) = strong_map(:,1:end-1);
            leftup_nb = false(row,col); leftup_nb(1:end-1,1:end-1) = strong_map(2:end,2:end);
            rightdown_nb = false(row,col); rightdown_nb(2:end,2:end) = strong_map(1:end-1,1:end-1);
            up_nb = false(row,col); up_nb(1:end-1,:) = strong_map(2:end,:);
            down_nb = false(row,col); down_nb(2:end,:) = strong_map(1:end-1,:);
            leftdown_nb = false(row,col); leftdown_nb(2:end,1:end-1) = strong_map(1:end-1,2:end);
            rightup_nb = false(row,col); rightup_nb(1:end-1,2:end) = strong_map(2:end,1:end-1);

            weak_nb = (left_nb|right_nb|up_nb|down_nb|leftup_nb|rightdown_nb|rightup_nb|leftdown_nb) & weak_map;

            linkedMap = linkedMap | weak_nb;    % add selected weak points to final edge map
            weak_nb_ind = find(weak_nb == 1);
            weak_map(weak_nb_ind) = 0;   % remove selected weak points from weak_map

            if isempty(weak_nb_ind)  % cannot link to other points, end iteration
                break;
            else
                strong_map = weak_nb;   % update strong_map to newly added points
            end

            itr = itr - 1;
        end
    end
end


end



